"""Phase 3: Promote staged files to final Lakehouse paths in S3.

Walks staged files in an S3 staging prefix (written by CTS after Phase 2),
uploads each to the final Lakehouse path with MD5 metadata from sidecar files,
archives replaced/suppressed and updated assemblies, and trims the transfer
manifest so that a re-run of Phase 2 only downloads remaining entries.
"""

import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import UTC, datetime
from logging import Logger, getLogger
from pathlib import Path, PurePosixPath
from typing import Any

import botocore.exceptions
import tqdm

from cdm_data_loaders.ncbi_ftp.constants import ACCESSION_PARTS_REGEX, ASSEMBLY_PATH_REGEX
from cdm_data_loaders.ncbi_ftp.metadata import (
    DescriptorResource,
    archive_descriptor,
    build_descriptor_key,
    create_descriptor,
    upload_descriptor,
)
from cdm_data_loaders.utils.file_transfer.s3 import client
from cdm_data_loaders.utils.file_transfer.s3.object_utils import (
    copy_object,
    delete_objects,
    list_objects,
    object_exists,
    upload_file,
)

logger: Logger = getLogger(__name__)

_MAX_DRY_RUN_LOGS = 10


# Promote from S3 staging prefix


def promote_from_s3(  # noqa: PLR0913
    *,
    staging_bucket: PurePosixPath,
    staging_key_prefix: PurePosixPath,
    lakehouse_bucket: PurePosixPath,
    lakehouse_key_prefix: PurePosixPath,
    removed_manifest_path: Path | None = None,
    updated_manifest_path: Path | None = None,
    manifest_s3_key: PurePosixPath | None = None,
    ncbi_release: str | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Promote files from an S3 staging prefix to the final Lakehouse path.

    Downloads each file to a temp location and re-uploads to the final path
    with MD5 metadata from ``.md5`` sidecar files.

    :param staging_bucket: S3 bucket containing the staged files (e.g. ``"cts"``)
    :param staging_key_prefix: S3 key prefix where CTS output was written
    :param lakehouse_bucket: S3 bucket for the final Lakehouse destination (e.g. ``"cdm-lake"``)
    :param lakehouse_key_prefix: S3 key prefix for final Lakehouse locations
    :param removed_manifest_path: local path to the removed_manifest file
    :param updated_manifest_path: local path to the updated_manifest file
    :param manifest_s3_key: S3 object key for transfer_manifest.txt (for trimming)
    :param ncbi_release: NCBI release version tag for archiving
    :param dry_run: if True, log actions without side effects
    :return: report dict with counts
    """
    # Get list of objects under the staging prefix
    staged_objects: list[dict[str, Any]] = list_objects(f"{staging_bucket / staging_key_prefix}/")

    # Separate data files from sidecars
    sidecars = {PurePosixPath(k["Key"]) for k in staged_objects if k["Key"].endswith((".crc64nvme", ".md5"))}
    data_files = [PurePosixPath(k["Key"]) for k in staged_objects if PurePosixPath(k["Key"]) not in sidecars]

    logger.info("Found %d data files and %d sidecars in staging", len(data_files), len(sidecars))

    # Archive all affected assemblies BEFORE promoting or deleting
    archived = 0
    for manifest_file, reason, delete in [
        (updated_manifest_path, "updated", False),
        (removed_manifest_path, "replaced_or_suppressed", True),
    ]:
        if manifest_file and Path(str(manifest_file)).is_file():
            archived += _archive_assemblies(
                manifest_file,
                lakehouse_bucket=lakehouse_bucket,
                lakehouse_key_prefix=lakehouse_key_prefix,
                ncbi_release=ncbi_release,
                archive_reason=reason,
                delete_source=delete,
                dry_run=dry_run,
            )

    promoted, failed, descriptors_written, promoted_accessions = _promote_data_files(
        data_files,
        sidecars,
        staging_key_prefix,
        lakehouse_key_prefix,
        staging_bucket,
        lakehouse_bucket,
        dry_run=dry_run,
    )

    # Trim manifest for resumability
    if manifest_s3_key and promoted_accessions and not dry_run:
        _trim_manifest(manifest_s3_key, staging_bucket, promoted_accessions)

    if descriptors_written:
        logger.info("Wrote %d frictionless descriptor(s)", descriptors_written)

    report: dict[str, Any] = {
        "timestamp": datetime.now(UTC).isoformat(),
        "promoted": promoted,
        "archived": archived,
        "failed": failed,
        "dry_run": dry_run,
    }

    logger.info(
        "PROMOTE SUMMARY: %d promoted, %d archived, %d failed%s",
        promoted,
        archived,
        failed,
        " (dry-run)" if dry_run else "",
    )
    return report


# Promote data files (per-file loop)


def _group_files_by_assembly(
    data_files: list[PurePosixPath], staging_prefix: PurePosixPath
) -> defaultdict[tuple[PurePosixPath, str], list[PurePosixPath]]:
    """Group files by assembly; skip download_report.json and non-raw_data paths.

    :param data_files: list of S3 keys for staged data files
    :param staging_prefix: staging prefix in S3
    :return: dict mapping (assembly_dir, accession) to list of staged keys for that assembly
    """
    assembly_files: defaultdict[tuple[PurePosixPath, str], list[PurePosixPath]] = defaultdict(list)
    for staged_key in data_files:
        if staged_key.match("**/download_report.json"):
            continue
        rel_path = staged_key.relative_to(staging_prefix)
        if not rel_path.is_relative_to("raw_data/"):
            continue
        m = ASSEMBLY_PATH_REGEX.search(str(staged_key))
        if m:
            assembly_files[(PurePosixPath(m.group(1)), m.group(2))].append(staged_key)
    return assembly_files


def _promote_file(  # noqa: PLR0913
    staged_key: PurePosixPath,
    staging_prefix: PurePosixPath,
    lakehouse_key_prefix: PurePosixPath,
    staging_bucket: PurePosixPath,
    lakehouse_bucket: PurePosixPath,
    sidecars: set[PurePosixPath],
) -> tuple[DescriptorResource, PurePosixPath]:
    """Download one staged file, re-upload to Lakehouse with MD5 metadata.

    :param staged_key: S3 key of the staged file
    :param staging_prefix: staging prefix in S3
    :param lakehouse_key_prefix: S3 key prefix for final Lakehouse locations
    :param staging_bucket: S3 bucket containing the staged file
    :param lakehouse_bucket: S3 bucket for the final Lakehouse destination
    :param sidecars: set of S3 keys for sidecar files (to check for MD5 metadata)
    :return: ``(resource_dict, staged_key)`` on success; raises on failure.
    """
    s3 = client.get_s3_client()
    rel_path = staged_key.relative_to(staging_prefix)
    final_key = lakehouse_key_prefix / rel_path
    final_key_path = PurePosixPath(final_key)

    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        tmp_path = tmp.name
    try:
        s3.download_file(Bucket=str(staging_bucket), Key=str(staged_key), Filename=tmp_path)

        metadata: dict[str, str] = {}
        md5_key = staged_key.with_name(f"{staged_key.name}.md5")
        if md5_key in sidecars:
            md5_obj = s3.get_object(Bucket=str(staging_bucket), Key=str(md5_key))
            metadata["md5"] = md5_obj["Body"].read().decode().strip()

        upload_succeeded = upload_file(
            tmp_path,
            str(lakehouse_bucket / final_key_path.parent),
            user_metadata=metadata,
            object_name=final_key_path.name,
            show_progress=False,
        )
        if not upload_succeeded:
            msg = f"upload_file returned False for {staged_key}"
            raise RuntimeError(msg)

        fname = final_key_path.name
        ext = fname.rsplit(".", 1)[-1] if "." in fname else ""
        resource: DescriptorResource = {
            "name": fname.lower(),
            "path": final_key,
            "format": ext,
            "bytes": Path(tmp_path).stat().st_size,
            "hash": metadata.get("md5"),
        }
        return resource, staged_key
    finally:
        Path(tmp_path).unlink()


def _write_descriptor_for_assembly(
    assembly_dir: PurePosixPath,
    accession: str,
    resources: list[DescriptorResource],
    lakehouse_bucket: PurePosixPath,
    lakehouse_key_prefix: PurePosixPath,
) -> bool:
    """Create and upload a frictionless descriptor for an assembly.

    :param assembly_dir: full assembly directory name
    :param accession: assembly accession (e.g. "GCF_000001405.39")
    :param resources: list of DescriptorResource dicts for the assembly's files
    :param lakehouse_bucket: S3 bucket for the final Lakehouse destination
    :param lakehouse_key_prefix: S3 key prefix for final Lakehouse locations
    :return: True if the descriptor was successfully written, False otherwise
    """
    try:
        descriptor_key = build_descriptor_key(assembly_dir, lakehouse_key_prefix)
        if object_exists(str(lakehouse_bucket / descriptor_key)):
            logger.debug("Descriptor already exists, skipping: %s", descriptor_key)
        else:
            descriptor = create_descriptor(assembly_dir, accession, resources)
            descriptor_key = upload_descriptor(
                descriptor, assembly_dir, lakehouse_bucket, lakehouse_key_prefix, dry_run=False
            )
            logger.debug("Uploaded descriptor: %s", descriptor_key)
            return True
    except Exception:
        logger.exception("Failed to write descriptor for %s", assembly_dir)
    return False


def _batch_delete(
    promoted_keys: list[PurePosixPath],
    sidecars: set[PurePosixPath],
    staging_bucket: PurePosixPath,
) -> None:
    """Batch-delete all staged data files and their sidecars in one API call.

    :param promoted_keys: list of staged keys that were successfully promoted
    :param sidecars: set of all sidecar keys in staging (to check for existence of sidecars for promoted files)
    :param staging_bucket: S3 bucket containing the staged files
    """
    keys_to_delete = list(promoted_keys)
    keys_to_delete.extend(
        key.with_name(f"{key.name}{sidecar_ext}")
        for key in promoted_keys
        for sidecar_ext in (".md5", ".crc64nvme")
        if key.with_name(f"{key.name}{sidecar_ext}") in sidecars
    )
    string_keys_to_delete = [str(key) for key in keys_to_delete]
    del_errors = delete_objects(str(staging_bucket), string_keys_to_delete)
    for err in del_errors:
        logger.warning("Failed to delete staged file %s: %s", err.get("Key"), err.get("Message"))


def _promote_data_files(  # noqa: PLR0913
    data_files: list[PurePosixPath],
    sidecars: set[PurePosixPath],
    staging_prefix: PurePosixPath,
    lakehouse_key_prefix: PurePosixPath,
    staging_bucket: PurePosixPath,
    lakehouse_bucket: PurePosixPath,
    *,
    dry_run: bool,
) -> tuple[int, int, int, set[str]]:
    """Promote each data file from staging to the final Lakehouse path.

    Files are grouped by assembly.  When all files for an assembly are promoted
    successfully, the frictionless descriptor is written immediately and the staged
    files (including sidecars) are deleted from staging.  This prevents staging
    accumulation across runs and ensures partial runs leave descriptors for all
    completed assemblies.

    :return: (promoted_count, failed_count, descriptors_written, promoted_accessions)
    """
    promoted = 0
    failed = 0
    descriptors_written = 0
    promoted_accessions: set[str] = set()
    assembly_files = _group_files_by_assembly(data_files, staging_prefix)

    def _promote_one(staged_key: PurePosixPath) -> tuple[DescriptorResource, PurePosixPath]:
        return _promote_file(
            staged_key,
            staging_prefix,
            lakehouse_key_prefix,
            staging_bucket,
            lakehouse_bucket,
            sidecars,
        )

    total_files = sum(len(v) for v in assembly_files.values())
    _dry_run_log_count = 0
    with tqdm.tqdm(total=total_files, unit="file", desc="Promoting") as pbar:
        for (adir, acc), files in assembly_files.items():
            assembly_failed = 0
            resources: list[DescriptorResource] = []
            promoted_keys: list[PurePosixPath] = []

            if dry_run:
                for staged_key in files:
                    rel_path = staged_key.relative_to(staging_prefix)
                    final_key = lakehouse_key_prefix / rel_path
                    if _dry_run_log_count < _MAX_DRY_RUN_LOGS:
                        logger.info("[dry-run] would promote: %s -> %s", staged_key, final_key)
                    else:
                        logger.debug("[dry-run] would promote: %s -> %s", staged_key, final_key)
                    _dry_run_log_count += 1
                    promoted += 1
                    pbar.update(1)
                continue

            # Download and re-upload all files for this assembly concurrently
            n_workers = min(32, len(files))
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                futures = {executor.submit(_promote_one, key): key for key in files}
                for future in as_completed(futures):
                    staged_key = futures[future]
                    try:
                        resource, _ = future.result()
                        resources.append(resource)
                        promoted_keys.append(staged_key)
                        promoted += 1
                        promoted_accessions.add(acc)
                    except Exception:
                        logger.exception("Failed to promote %s", staged_key)
                        assembly_failed += 1
                    pbar.update(1)

            failed += assembly_failed

            # Write descriptor and delete staged files immediately after a fully successful assembly
            if assembly_failed == 0 and promoted_keys:
                if _write_descriptor_for_assembly(adir, acc, resources, lakehouse_bucket, lakehouse_key_prefix):
                    descriptors_written += 1

                # delete staged files in batch with their sidecars (if any)
                _batch_delete(promoted_keys, sidecars, staging_bucket)

    return promoted, failed, descriptors_written, promoted_accessions


# Archive assemblies


def _get_accession_path_prefix(accession: str, lakehouse_key_prefix: PurePosixPath) -> PurePosixPath | None:
    """Get the S3 key prefix for all files related to an accession.

    :param accession: assembly accession (e.g. "GCF_000001405.39_Some_description")
    :param lakehouse_key_prefix: S3 key prefix for the Lakehouse dataset root
    :return: S3 key prefix under which all files for the accession are stored, or None if the accession format is invalid
    """
    m = ACCESSION_PARTS_REGEX.match(accession)  # returns, e.g., "GCF", "000", "001"', "405" captured groups
    if not m:
        logger.warning("Invalid accession format: %s", accession)
        return None
    return PurePosixPath(lakehouse_key_prefix, "raw_data", *m.groups(), accession)


def _get_source_dest_pairs_for_accession(
    accession: str,
    lakehouse_bucket: PurePosixPath,
    lakehouse_key_prefix: PurePosixPath,
    release_tag: str,
    archive_reason: str,
) -> list[tuple[PurePosixPath, PurePosixPath]]:
    """Get list of (source_key, archive_key) pairs for all objects related to an accession.

    :param accession: assembly accession (e.g. "GCF_000001405.39_Some_description")
    :param lakehouse_bucket: S3 bucket for the Lakehouse (source and archive destination)
    :param lakehouse_key_prefix: S3 key prefix for the Lakehouse dataset root
    :param release_tag: release tag for the archive
    :param archive_reason: reason for archiving
    :return: list of (source_key, archive_key) pairs for all objects related to the accession
    """
    source_prefix = _get_accession_path_prefix(accession, lakehouse_key_prefix)
    if not source_prefix:
        return []
    matching_objs: list[dict[str, Any]] = list_objects(f"{lakehouse_bucket / source_prefix}")
    return [
        (
            PurePosixPath(obj["Key"]),
            lakehouse_key_prefix
            / "archive"
            / release_tag
            / archive_reason
            / PurePosixPath(obj["Key"]).relative_to(lakehouse_key_prefix),
        )
        for obj in matching_objs
    ]


def _dry_run_output(key_pairs: list[tuple[PurePosixPath, PurePosixPath]], log_count: int) -> int:
    """Log source and archive key pairs for a dry run, with a limit on how many are logged at INFO level.

    :param key_pairs: list of (source_key, archive_key) pairs
    :param log_count: current count of logged entries
    :return: updated count of logged entries
    """
    _dry_run_log_count = log_count
    for source_key, archive_key in key_pairs:
        if _dry_run_log_count < _MAX_DRY_RUN_LOGS:
            logger.info("[dry-run] would archive: %s -> %s", source_key, archive_key)
        else:
            logger.debug("[dry-run] would archive: %s -> %s", source_key, archive_key)
        _dry_run_log_count += 1
    return _dry_run_log_count


def _archive_objects(
    key_pairs: list[tuple[PurePosixPath, PurePosixPath]], lakehouse_bucket: PurePosixPath, *, delete_source: bool
) -> int:
    """Copy objects from source keys to archive keys, optionally deleting the source objects.

    :param key_pairs: list of (source_key, archive_key) pairs
    :param lakehouse_bucket: S3 bucket for the Lakehouse (source and archive destination)
    :param delete_source: if True, delete the source object after copying
    :return: number of objects successfully archived
    """
    archived = 0
    if not key_pairs:
        return archived
    keys_to_delete: list[str] = []  # strings so they can be passed to delete_objects() in s3 module
    n_workers = min(32, len(key_pairs))
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(
                copy_object,
                str(lakehouse_bucket / src),
                str(lakehouse_bucket / arch),
            ): src
            for src, arch in key_pairs
        }
        for future in as_completed(futures):
            src = futures[future]
            try:
                future.result()
                archived += 1
                if delete_source:
                    keys_to_delete.append(str(src))
                logger.debug("  Archived: %s", src)
            except Exception:
                logger.exception("Failed to archive %s", src)

    if delete_source and keys_to_delete:
        del_errors = delete_objects(str(lakehouse_bucket), keys_to_delete)
        for err in del_errors:
            logger.warning("Failed to delete %s: %s", err.get("Key"), err.get("Message"))

    return archived


def _archive_assemblies(  # noqa: PLR0913
    manifest_local_path: Path,
    lakehouse_bucket: PurePosixPath,
    lakehouse_key_prefix: PurePosixPath,
    ncbi_release: str | None = None,
    archive_reason: str = "unknown",
    *,
    delete_source: bool = False,
    dry_run: bool = False,
) -> int:
    """Archive assembly objects to ``archive/{release_tag}/``.

    Copies S3 objects matching each accession to the archive prefix.
    When *delete_source* is True (replaced/suppressed), the original
    objects are deleted after copying.  When False (updated), the
    originals remain in place to be overwritten by the promote step.

    :param manifest_local_path: local path to a manifest file (one accession per line)
    :param lakehouse_bucket: S3 bucket for the Lakehouse (source and archive destination)
    :param ncbi_release: release tag used in the archive path
    :param lakehouse_key_prefix: S3 key prefix for the Lakehouse dataset root
    :param archive_reason: metadata value describing why the object was archived
    :param delete_source: if True, delete the source object after copying
    :param dry_run: if True, log without making changes
    :return: number of objects archived
    """
    release_tag = ncbi_release or "unknown"
    archived = 0

    with Path(manifest_local_path).open() as f:
        accessions = [line.strip() for line in f if line.strip()]

    _dry_run_log_count = 0
    for accession in tqdm.tqdm(accessions, unit="accession", desc="Archiving"):
        # get list of (source_key, archive_key) pairs for all objects related to this accession
        key_pairs: list[tuple[PurePosixPath, PurePosixPath]] = _get_source_dest_pairs_for_accession(
            accession,
            lakehouse_bucket,
            lakehouse_key_prefix,
            release_tag,
            archive_reason,
        )
        if not key_pairs:
            logger.debug("No objects found for %s, skipping archive", accession)
            continue

        # Archive all files for this accession
        if dry_run:
            _dry_run_log_count = _dry_run_output(key_pairs, _dry_run_log_count)
            archived += len(key_pairs)
            continue
        archived += _archive_objects(key_pairs, lakehouse_bucket, delete_source=delete_source)

        # Infer assembly_dir from key paths for descriptor archival
        assembly_dir: PurePosixPath | None = None
        for src, _ in key_pairs:
            if adir_match := ASSEMBLY_PATH_REGEX.search(f"{src}"):
                assembly_dir = PurePosixPath(adir_match.group(1))
                break

        # Archive the frictionless descriptor alongside raw data
        if assembly_dir:
            try:
                archived_desc = archive_descriptor(
                    assembly_dir,
                    lakehouse_bucket,
                    lakehouse_key_prefix,
                    release_tag,
                    archive_reason=archive_reason,
                    dry_run=dry_run,
                )
                if not archived_desc:
                    logger.debug("No descriptor found to archive for %s", assembly_dir)
            except Exception:
                logger.exception("Failed to archive descriptor for %s", assembly_dir)

    logger.info("Archived %d objects for %d accessions (%s)", archived, len(accessions), archive_reason)
    return archived


# Manifest trimming


def _trim_manifest(
    manifest_s3_key: PurePosixPath, staging_bucket: PurePosixPath, promoted_accessions: set[str]
) -> None:
    """Remove promoted accessions from the transfer manifest in S3.

    :param manifest_s3_key: S3 object key of the transfer_manifest.txt
    :param staging_bucket: S3 bucket containing the transfer manifest
    :param promoted_accessions: set of accessions that were successfully promoted
    """
    s3 = client.get_s3_client()

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".txt") as tmp:
        tmp_path = tmp.name

    try:
        try:
            s3.download_file(Bucket=str(staging_bucket), Key=str(manifest_s3_key), Filename=tmp_path)
        except s3.exceptions.NoSuchKey:
            logger.warning("Manifest not found in S3 (s3://%s/%s) — skipping trim", staging_bucket, manifest_s3_key)
            return
        except botocore.exceptions.ClientError as e:
            if e.response["Error"]["Code"] == "404":
                logger.warning("Manifest not found in S3 (s3://%s/%s) — skipping trim", staging_bucket, manifest_s3_key)
                return
            raise

        with Path(tmp_path).open() as f:
            lines = f.readlines()

        remaining = [line for line in lines if line.strip() and not any(acc in line for acc in promoted_accessions)]

        with Path(tmp_path).open("w") as f:
            f.writelines(remaining)

        s3.upload_file(Filename=tmp_path, Bucket=str(staging_bucket), Key=str(manifest_s3_key))
        logger.info(
            "Trimmed manifest: %d -> %d entries (%d promoted)",
            len(lines),
            len(remaining),
            len(lines) - len(remaining),
        )
    finally:
        Path(tmp_path).unlink()
