"""NCBI FTP assembly download pipeline (Phase 2).

Orchestrates parallel downloading of NCBI assemblies listed in a transfer
manifest.  Settings, batching, CLI entry point, and CTS integration live here;
domain-specific download logic is in :mod:`cdm_data_loaders.ncbi_ftp.assembly`.
"""

import json
import logging
import shutil
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from ftplib import error_temp
from logging import Logger, getLogger
from pathlib import Path, PurePosixPath
from typing import Any

import tqdm
from pydantic import AliasChoices, Field
from tenacity import before_sleep_log, retry, retry_if_exception_type, stop_after_attempt, wait_exponential

from cdm_data_loaders.core.fields import INPUT_MOUNT, OUTPUT_MOUNT
from cdm_data_loaders.core.settings import LoggerSettings
from cdm_data_loaders.ncbi_ftp.assembly import (
    FTP_HOST,
    build_accession_path,
    download_assembly_to_local,
    parse_assembly_path,
)
from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.utils.file_transfer.ftp_client import ThreadLocalFTP
from cdm_data_loaders.utils.file_transfer.s3.client import get_s3_client
from cdm_data_loaders.utils.file_transfer.s3.object_utils import upload_file

logger: Logger = getLogger(__name__)


DEFAULT_STAGING_KEY_PREFIX: PurePosixPath = PurePosixPath("staging")


class DownloadSettings(LoggerSettings):
    """Configuration for the NCBI FTP assembly download pipeline."""

    manifest: Path = Field(
        default=Path(INPUT_MOUNT) / "transfer_manifest.txt",
        description="Path to the transfer manifest file listing FTP paths to download",
        validation_alias=AliasChoices("m", "manifest"),
    )
    output_dir: Path = Field(
        default=Path(OUTPUT_MOUNT),
        description="Output directory for downloaded assembly files",
        validation_alias=AliasChoices("output-dir", "output_dir"),
    )
    threads: int = Field(
        default=4,
        ge=1,
        le=32,
        description="Number of parallel download threads",
        validation_alias=AliasChoices("t", "threads"),
    )
    ftp_host: str = Field(
        default=FTP_HOST,
        description="NCBI FTP hostname",
        validation_alias=AliasChoices("ftp-host", "ftp_host"),
    )
    limit: int | None = Field(
        default=None,
        ge=1,
        description="Limit to first N assemblies (for testing)",
        validation_alias=AliasChoices("l", "limit"),
    )


# Private helpers


def _upload_assembly_dir(
    assembly_dir: Path,
    tmp_root: Path,
    bucket: PurePosixPath,
    staging_key_prefix: PurePosixPath,
) -> int:
    """Upload all files under *assembly_dir* to S3, deleting each file immediately after upload.

    Empty directories are removed after all files are uploaded.  If the
    directory does not exist (e.g. the assembly had no files) the function
    returns zero without raising.

    :param assembly_dir: local directory for one assembly
    :param tmp_root: root of the temp directory (used to compute relative S3 paths)
    :param bucket: destination S3 bucket
    :param staging_key_prefix: S3 key prefix within *bucket*
    :return: number of files uploaded
    """
    if not assembly_dir.exists():
        return 0
    count = 0
    for f in sorted(assembly_dir.rglob("*")):
        if f.is_file():
            relative = f.relative_to(tmp_root)
            dest_prefix = bucket / staging_key_prefix / relative.parent
            if upload_file(f, str(dest_prefix), show_progress=False):
                count += 1
            else:
                logger.warning("Failed to upload %s to %s", f, dest_prefix)
            f.unlink()
    shutil.rmtree(assembly_dir, ignore_errors=True)
    return count


# Batch download


def download_batch(
    manifest_path: Path,
    output_dir: Path,
    threads: int = 4,
    ftp_host: str = FTP_HOST,
    limit: int | None = None,
) -> dict[str, Any]:
    """Download all assemblies listed in the manifest.

    :param manifest_path: path to the transfer manifest file
    :param output_dir: base output directory
    :param threads: number of parallel download threads
    :param ftp_host: FTP hostname
    :param limit: optional limit for testing
    :return: report dict with overall stats
    """
    with manifest_path.open() as f:
        assembly_paths = [PurePosixPath(line.strip()) for line in f if line.strip() and not line.startswith("#")]

    if limit:
        assembly_paths = assembly_paths[:limit]

    logger.info("Starting download of %d assemblies with %d threads", len(assembly_paths), threads)

    pool = ThreadLocalFTP(ftp_host)
    lock = threading.Lock()
    success_count = 0
    failed: list[dict[str, str]] = []
    all_stats: list[dict[str, Any]] = []

    def _download_one(path: PurePosixPath) -> tuple[str, Exception | None]:
        nonlocal success_count

        @retry(
            retry=retry_if_exception_type((error_temp, BrokenPipeError, ConnectionResetError, EOFError)),
            stop=stop_after_attempt(5),
            wait=wait_exponential(multiplier=1, min=5, max=60),
            reraise=True,
            before_sleep=before_sleep_log(logger, logging.WARNING),
        )
        def _attempt() -> dict[str, Any]:
            return download_assembly_to_local(path, output_dir, ftp_host=ftp_host, ftp=pool.get())

        try:
            stats = _attempt()
        except Exception as e:  # noqa: BLE001
            return str(path), e
        else:
            with lock:
                success_count += 1
                all_stats.append(stats)
            return str(path), None

    try:
        with (
            tqdm.tqdm(
                total=len(assembly_paths), unit="assembly", desc="Downloading from NCBI FTP", smoothing=0.01
            ) as pbar,
            ThreadPoolExecutor(max_workers=threads) as executor,
        ):
            futures = {executor.submit(_download_one, p): p for p in assembly_paths}
            for future in as_completed(futures):
                path, error = future.result()
                if error:
                    logger.error("FAILED: %s: %s", path, error)
                    with lock:
                        failed.append({"path": path, "error": str(error)})
                pbar.update(1)
    finally:
        pool.close_all()

    report: dict[str, Any] = {
        "timestamp": datetime.now(UTC).isoformat(),
        "total_attempted": len(assembly_paths),
        "succeeded": success_count,
        "failed": len(failed),
        "failures": failed,
        "assembly_stats": all_stats,
    }

    report_path = output_dir / "download_report.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as f:
        json.dump(report, f, indent=2)
    logger.info("Download report written to: %s", report_path)

    logger.info(
        "SUMMARY: %d attempted, %d succeeded, %d failed",
        len(assembly_paths),
        success_count,
        len(failed),
    )

    return report


# CTS entry point


def run_download(config: DownloadSettings) -> None:
    """Main CTS entry point for Phase 2 download.

    :param config: validated download settings
    """
    report = download_batch(
        manifest_path=config.manifest,
        output_dir=config.output_dir,
        threads=config.threads,
        ftp_host=config.ftp_host,
        limit=config.limit,
    )
    if report["failed"] > 0:
        msg = f"Download completed with {report['failed']} failures"
        raise RuntimeError(msg)


# Notebook / interactive entry point


def download_and_stage(  #  noqa: PLR0913, PLR0915
    *,
    bucket: PurePosixPath,
    staging_key_prefix: PurePosixPath,
    manifest_s3_key: PurePosixPath | None = None,
    manifest_local_path: Path | None = None,
    threads: int = 4,
    ftp_host: str = FTP_HOST,
    limit: int | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Download assemblies from NCBI FTP and stage them to S3 (Phase 2).

    Exactly one of *manifest_s3_key* or *manifest_local_path* must be given.

    Downloads and uploads are pipelined per assembly: each worker downloads one
    assembly, immediately uploads its files to S3, then deletes the local copies
    before picking up the next assembly.  At most *threads* assembly directories
    exist on disk simultaneously, preventing disk exhaustion on large batches.

    :param bucket: destination S3 bucket name
    :param staging_key_prefix: key prefix inside the bucket (e.g. ``"staging/run1/"``)
    :param manifest_s3_key: S3 object key of the transfer manifest within *bucket*
    :param manifest_local_path: local path to the transfer manifest file
    :param threads: number of parallel download-and-upload workers
    :param ftp_host: NCBI FTP hostname
    :param limit: optional limit for testing
    :param dry_run: when ``True``, download but skip all S3 uploads
    :return: download report extended with ``staged_objects``, ``staging_key_prefix``, ``dry_run``
    """
    if manifest_s3_key is not None and manifest_local_path is not None:
        msg = "Provide exactly one of manifest_s3_key or manifest_local_path, not both"
        raise ValueError(msg)
    if manifest_s3_key is None and manifest_local_path is None:
        msg = "One of manifest_s3_key or manifest_local_path must be provided"
        raise ValueError(msg)

    with tempfile.TemporaryDirectory() as _tmpdir:
        tmp = Path(_tmpdir)
        manifest_dest = tmp / "transfer_manifest.txt"

        if manifest_s3_key is not None:
            s3 = get_s3_client()
            response = s3.get_object(Bucket=str(bucket), Key=str(manifest_s3_key))
            manifest_dest.write_bytes(response["Body"].read())
            logger.info("Manifest read from S3: s3://%s/%s", bucket, manifest_s3_key)
        elif manifest_local_path:
            manifest_dest.write_bytes(manifest_local_path.read_bytes())
            logger.info("Manifest read from local path: %s", manifest_local_path)

        with manifest_dest.open() as f:
            assembly_paths = [PurePosixPath(line.strip()) for line in f if line.strip() and not line.startswith("#")]

        if limit:
            assembly_paths = assembly_paths[:limit]

        logger.info("Starting download & stage of %d assemblies with %d threads", len(assembly_paths), threads)

        pool = ThreadLocalFTP(ftp_host)
        lock = threading.Lock()
        success_count = 0
        staged_objects = 0
        failed: list[dict[str, str]] = []
        all_stats: list[dict[str, Any]] = []

        def _download_upload_one(path: PurePosixPath) -> tuple[str, Exception | None]:
            nonlocal success_count, staged_objects

            @retry(
                retry=retry_if_exception_type((error_temp, BrokenPipeError, ConnectionResetError, EOFError)),
                stop=stop_after_attempt(5),
                wait=wait_exponential(multiplier=1, min=5, max=60),
                reraise=True,
                before_sleep=before_sleep_log(logger, logging.WARNING),
            )
            def _attempt() -> dict[str, Any]:
                return download_assembly_to_local(path, tmp, ftp_host=ftp_host, ftp=pool.get())

            try:
                stats = _attempt()
            except Exception as e:  # noqa: BLE001
                return str(path), e

            if not dry_run:
                _db, assembly_dir_name, _accession = parse_assembly_path(path)
                assembly_local_dir = tmp / build_accession_path(assembly_dir_name)
                count = _upload_assembly_dir(assembly_local_dir, tmp, bucket, staging_key_prefix)
                with lock:
                    staged_objects += count

            with lock:
                success_count += 1
                all_stats.append(stats)
            return str(path), None

        try:
            with (
                tqdm.tqdm(
                    total=len(assembly_paths), unit="assembly", desc="Downloading & staging", smoothing=0.01
                ) as pbar,
                ThreadPoolExecutor(max_workers=threads) as executor,
            ):
                futures = {executor.submit(_download_upload_one, p): p for p in assembly_paths}
                for future in as_completed(futures):
                    path, error = future.result()
                    if error:
                        logger.error("FAILED: %s: %s", path, error)
                        with lock:
                            failed.append({"path": path, "error": str(error)})
                    pbar.update(1)
        finally:
            pool.close_all()

        report: dict[str, Any] = {
            "timestamp": datetime.now(UTC).isoformat(),
            "total_attempted": len(assembly_paths),
            "succeeded": success_count,
            "failed": len(failed),
            "failures": failed,
            "assembly_stats": all_stats,
        }

        if not dry_run:
            report_path = tmp / "download_report.json"
            with report_path.open("w") as f:
                json.dump(report, f, indent=2)
            if upload_file(report_path, str(bucket / staging_key_prefix), show_progress=False):
                staged_objects += 1
            else:
                logger.warning("Failed to upload download report to s3://%s/%s", bucket, staging_key_prefix)
            logger.info("Staged %d objects to s3://%s/%s", staged_objects, bucket, staging_key_prefix)

        logger.info(
            "SUMMARY: %d attempted, %d succeeded, %d failed",
            len(assembly_paths),
            success_count,
            len(failed),
        )

    return {
        **report,
        "staged_objects": staged_objects,
        "staging_key_prefix": staging_key_prefix,
        "dry_run": dry_run,
    }


cli = PipelineCli(DownloadSettings, run_download)  # pyright: ignore[reportArgumentType]


if __name__ == "__main__":
    cli()
