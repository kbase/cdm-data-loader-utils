"""Phase 1: Assembly summary diffing and manifest generation.

Downloads the current NCBI assembly summary from FTP, compares it against a
previous snapshot, and produces ``transfer_manifest.txt`` (assemblies to
download), ``removed_manifest.txt`` (assemblies to archive), and a JSON diff
summary.  All filtering logic (prefix range, limit) lives here so that
downstream phases receive a final, pre-filtered manifest.
"""

import contextlib
import csv
import json
import time
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from datetime import UTC, datetime
from http import HTTPStatus
from logging import Logger, getLogger
from pathlib import Path, PurePosixPath
from typing import Any
from urllib.parse import urlsplit

from botocore.exceptions import ClientError

from cdm_data_loaders.ncbi_ftp.assembly import (
    FILE_FILTERS,
    FTP_HOST,
    build_accession_path,
    parse_md5_checksums_file,
)
from cdm_data_loaders.ncbi_ftp.constants import ACCESSION_PARTS_REGEX, ASSEMBLY_PATH_REGEX
from cdm_data_loaders.utils.file_transfer.ftp_client import FTP, connect_ftp, ftp_noop_keepalive, ftp_retrieve_text
from cdm_data_loaders.utils.file_transfer.s3.object_utils import (
    head_object,
    list_objects,
)

logger: Logger = getLogger(__name__)

_DATABASE_ACC_PREFIX: dict[str, str] = {
    "refseq": "GCF_",
    "genbank": "GCA_",
}

SUMMARY_FTP_PATHS: dict[str, PurePosixPath] = {
    "refseq": PurePosixPath("/") / "genomes" / "ASSEMBLY_REPORTS" / "assembly_summary_refseq.txt",
    "genbank": PurePosixPath("/") / "genomes" / "ASSEMBLY_REPORTS" / "assembly_summary_genbank.txt",
}

# Assembly summary file columns of interest
_ACCESSION_COL = 0
_STATUS_COL = 10
_SEQ_REL_DATE_COL = 14
_FTP_URL_COL = 19
_MIN_COL = 20

# Data structures


@dataclass
class AssemblyRecord:
    """Parsed row from an NCBI assembly summary file.

    Attributes:
        accession: Assembly accession (e.g., "GCF_000001215.4").
        status: "latest", "replaced", or "suppressed".
        seq_rel_date: Release date.
        ftp_url: Full FTP URL.
        assembly_dir: Assembly directory name (final path segment from the FTP URL).

    """

    accession: str
    status: str
    seq_rel_date: str
    ftp_url: str
    assembly_dir: PurePosixPath


@dataclass
class DiffResult:
    """Result of comparing current and previous assembly summaries."""

    new: list[str] = field(default_factory=list)
    updated: list[str] = field(default_factory=list)
    replaced: list[str] = field(default_factory=list)
    suppressed: list[str] = field(default_factory=list)


# Assembly summary download & parsing


def download_assembly_summary(database: str = "refseq", ftp_host: str = FTP_HOST) -> str:
    """Download the assembly summary file from NCBI FTP.

    :param database: ``"refseq"`` or ``"genbank"``
    :param ftp_host: FTP hostname
    :return: raw text content of the summary file
    """
    ftp_path = SUMMARY_FTP_PATHS.get(database)
    if not ftp_path:
        msg = f"Unknown database: {database}"
        raise ValueError(msg)

    logger.info("Downloading assembly_summary_%s.txt from NCBI FTP ...", database)
    ftp = connect_ftp(ftp_host)
    try:
        content = ftp_retrieve_text(ftp, ftp_path)
    finally:
        with contextlib.suppress(Exception):
            ftp.quit()

    logger.info("Downloaded assembly summary (%d bytes)", len(content))
    return content


def parse_assembly_summary(source: str | Path | list[str]) -> dict[str, AssemblyRecord]:
    """Parse an NCBI assembly summary into a dict of assembly records.

    Accepts a file path, raw text string, or list of lines.

    Columns of interest (0-indexed):
        0: assembly_accession (e.g. GCF_000001215.4)
        10: version_status ("latest", "replaced", "suppressed")
        14: seq_rel_date
        19: ftp_path (full FTP URL or "na")

    :param source: file path, raw text, or list of lines
    :return: dict mapping accession to :class:`AssemblyRecord`
    """
    assemblies: dict[str, AssemblyRecord] = {}

    def _parse_lines(lines: Iterable[str]) -> None:
        reader = csv.reader(
            (line.rstrip("\n") for line in lines if not line.startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            if len(row) < _MIN_COL:
                continue
            accession = row[_ACCESSION_COL]
            ftp_url = row[_FTP_URL_COL]
            if "https://ftp" not in ftp_url:
                for col in row:
                    if "https://ftp" in col:
                        ftp_url = col
                        break
            if "https://ftp" not in ftp_url:
                msg = f"Missing ftp path for record {accession}."
                logger.warning(msg)
                continue
            assembly_dir = PurePosixPath(_ftp_dir_from_url(ftp_url).name)
            assemblies[accession] = AssemblyRecord(
                accession=accession,
                status=row[_STATUS_COL],
                seq_rel_date=row[_SEQ_REL_DATE_COL],
                ftp_url=ftp_url,
                assembly_dir=assembly_dir,
            )

    if isinstance(source, Path) or (isinstance(source, str) and "\n" not in source and Path(source).is_file()):
        with Path(source).open() as f:
            _parse_lines(f)
    elif isinstance(source, list):
        _parse_lines(source)
    else:
        _parse_lines(source.splitlines(keepends=True))

    logger.info("Parsed %d assemblies from summary", len(assemblies))
    return assemblies


def get_latest_assembly_paths(
    assemblies: dict[str, AssemblyRecord], ftp_host: str = FTP_HOST
) -> list[tuple[str, PurePosixPath]]:
    """Extract FTP directory paths for all assemblies with ``latest`` status.

    :param assemblies: parsed assembly records
    :param ftp_host: FTP hostname for URL stripping
    :return: list of ``(accession, ftp_dir_path)`` tuples
    """
    paths: list[tuple[str, PurePosixPath]] = []
    for accession, rec in assemblies.items():
        if rec.status != "latest":
            continue
        ftp_path = _ftp_dir_from_url(rec.ftp_url, ftp_host)
        paths.append((accession, ftp_path))
    return paths


# Prefix filtering


def accession_prefix(accession: str) -> str:
    """Extract the 3-digit prefix from an accession (e.g. ``GCF_000005845.2`` → ``"000"``)."""
    if m := ACCESSION_PARTS_REGEX.match(accession):
        return m.group(2)
    msg = f"Could not parse accession: {accession}"
    raise ValueError(msg)


def filter_by_prefix_range(
    assemblies: dict[str, AssemblyRecord],
    prefix_from: str | None = None,
    prefix_to: str | None = None,
) -> dict[str, AssemblyRecord]:
    """Filter assemblies to those whose 3-digit accession prefix is in range.

    Both bounds are inclusive. If neither is set, returns all assemblies.

    :param assemblies: dict of parsed assembly records
    :param prefix_from: lower bound (inclusive), e.g. ``"000"``
    :param prefix_to: upper bound (inclusive), e.g. ``"003"``
    :return: filtered dict
    """
    if prefix_from is None and prefix_to is None:
        return assemblies
    filtered: dict[str, AssemblyRecord] = {}
    for acc, rec in assemblies.items():
        pfx = accession_prefix(acc)
        if (
            pfx is None
            or (prefix_from is not None and pfx < prefix_from)
            or (prefix_to is not None and pfx > prefix_to)
        ):
            continue
        filtered[acc] = rec
    return filtered


# Diff computation


def compute_diff(
    current: dict[str, AssemblyRecord],
    previous_assemblies: dict[str, AssemblyRecord] | None = None,
    previous_accessions: set[str] | None = None,
) -> DiffResult:
    """Compute the diff between current and previous assembly state.

    :param current: the new NCBI summary (parsed)
    :param previous_assemblies: full parsed previous summary, or None if using fallback
    :param previous_accessions: set of known accessions (store-scan fallback)
    :return: diff result with new/updated/replaced/suppressed lists
    """
    diff = DiffResult()

    known = set(previous_assemblies) if previous_assemblies is not None else (previous_accessions or set())

    for acc, rec in current.items():
        if rec.status in ("replaced", "suppressed"):
            if acc in known:
                getattr(diff, rec.status).append(acc)
            continue
        if rec.status != "latest":
            continue

        if acc not in known:
            diff.new.append(acc)
        elif previous_assemblies is not None:
            prev = previous_assemblies.get(acc)
            if prev and (rec.seq_rel_date > prev.seq_rel_date or rec.assembly_dir != prev.assembly_dir):
                diff.updated.append(acc)

    # Accessions in previous but entirely absent from current (withdrawn)
    current_accs = set(current)
    diff.suppressed.extend(known - current_accs)

    diff.new.sort()
    diff.updated.sort()
    diff.replaced.sort()
    diff.suppressed.sort()
    return diff


# FTP URL helpers


def _ftp_dir_from_url(ftp_url: str, ftp_host: str = FTP_HOST) -> PurePosixPath:
    """Convert an FTP URL from the assembly summary to an FTP directory path."""
    parsed = urlsplit(ftp_url)

    # Full URL
    if parsed.scheme in {"ftp", "https"}:
        if parsed.netloc and parsed.netloc != ftp_host:
            msg = f"Unexpected FTP host: {parsed.netloc}"
            raise ValueError(msg)
        return PurePosixPath(parsed.path)

    # unparsable url
    msg = f"Could not parse FTP URL: {ftp_url}"
    raise ValueError(msg)


# Synthetic summary from S3 store scan


def _extract_accession_dir_and_id_from_s3_key(key: PurePosixPath) -> tuple[str, str]:
    """Extract both accession and assembly directory from an S3 object key.

    e.g. "some/prefix/GCF_000001215.4_Release_6_plus_ISO1_MT/file.gz"
         → ("GCF_000001215.4_Release_6_plus_ISO1_MT", "GCF_000001215.4")
    """
    if m := ASSEMBLY_PATH_REGEX.search(str(key)):
        return (m.group(1), m.group(2))
    msg = f"Could not parse S3 key for accession info: {key}"
    raise ValueError(msg)


def scan_store_to_synthetic_summary(
    bucket: PurePosixPath,
    key_prefix: PurePosixPath,
    release_date: str,
    database: str = "refseq",
    progress_callback: Callable[[int, str], None] | None = None,
) -> dict[str, AssemblyRecord]:
    """Scan S3 store and build a synthetic assembly summary from existing objects.

    This function is useful when bootstrapping a diffs against an existing,
    pre-populated S3 store that lacks a baseline assembly summary.

        For each assembly found in the store:
    - Extracts the accession and assembly directory name from S3 paths
        - Applies the provided ``release_date`` as synthetic ``seq_rel_date`` for
            all assemblies
    - Creates an ``AssemblyRecord`` with ``status="latest"``
    - Filters to accessions matching the expected prefix for ``database``
      (``GCF_`` for ``"refseq"``, ``GCA_`` for ``"genbank"``)

    The function paginates through S3 to handle large stores efficiently.

    :param bucket: S3 bucket name
    :param key_prefix: S3 key prefix (all objects under this prefix are scanned)
    :param release_date: release date string in ``YYYY/MM/DD`` format used for
        all synthetic records
    :param database: ``"refseq"`` or ``"genbank"`` — controls which accession
        prefix is included (``GCF_`` or ``GCA_`` respectively)
    :param progress_callback: optional callable invoked after each accession is
        processed with ``(count, accession)`` where count is the running total
        of unique accessions found
    :return: dict mapping accession to ``AssemblyRecord``
    """
    try:
        datetime.strptime(release_date, "%Y/%m/%d").astimezone(UTC)
    except ValueError as exc:
        msg = f"Invalid release_date '{release_date}'. Expected format YYYY/MM/DD."
        raise ValueError(msg) from exc

    acc_prefix = _DATABASE_ACC_PREFIX.get(database)
    if acc_prefix is None:
        msg = f"Unknown database: {database!r}. Expected 'refseq' or 'genbank'."
        raise ValueError(msg)

    assemblies: dict[str, AssemblyRecord] = {}
    processed_count = 0

    try:
        objs = list_objects(str(bucket / key_prefix))

        for obj in objs:
            try:
                assembly_dir, acc = _extract_accession_dir_and_id_from_s3_key(obj["Key"])
            except ValueError:
                continue

            if not acc.startswith(acc_prefix):
                continue

            if acc not in assemblies:
                # First object for this accession; store it.
                # Construct a fake FTP path that ends with assembly_dir so
                # that round-tripping through parse_assembly_summary (which
                # derives assembly_dir via ftp_path.rstrip("/").split("/")[-1])
                # yields the correct assembly_dir and therefore correct diffs.
                fake_ftp_path = f"https://ftp.ncbi.nlm.nih.gov/synthetic/{assembly_dir}"
                assemblies[acc] = AssemblyRecord(
                    accession=acc,
                    status="latest",
                    seq_rel_date=release_date,
                    ftp_url=fake_ftp_path,
                    assembly_dir=PurePosixPath(assembly_dir),
                )
                processed_count += 1
                if progress_callback is not None:
                    progress_callback(processed_count, acc)

    except Exception:
        logger.exception("Error scanning store")
        raise

    logger.info("Scanned S3 store: found %d unique assemblies", len(assemblies))
    return assemblies


# Checksum verification against S3 store


def _fetch_accession_checksums_from_ftp(
    ftp: FTP | None, ftp_host: str, last_activity: float, current_accession: str
) -> tuple[dict[PurePosixPath, str], float]:
    """Fetch and parse the md5checksums.txt file for a given accession from FTP.

    :param ftp: FTP connection object
    :param ftp_host: NCBI FTP hostname
    :param last_activity: timestamp of the last FTP activity
    :param current_accession: accession identifier
    :return: dictionary mapping file names to MD5 checksums and updated last_activity timestamp
    """
    if ftp is None:
        ftp = connect_ftp(ftp_host)
    last_activity = ftp_noop_keepalive(ftp, last_activity)
    ftp_dir = _ftp_dir_from_url(current_accession, ftp_host)
    try:
        md5_text = ftp_retrieve_text(ftp, ftp_dir / "md5checksums.txt")
        last_activity = time.monotonic()
        ftp_checksums = parse_md5_checksums_file(md5_text)
    except Exception:  # noqa: BLE001
        logger.warning("Cannot fetch md5checksums.txt for %s, keeping in transfer list", current_accession)
        return {}, last_activity
    return {
        fname: md5 for fname, md5 in ftp_checksums.items() if any(fname.match(f"**{suffix}") for suffix in FILE_FILTERS)
    }, last_activity


def _does_accession_need_update(
    target_checksums: dict[PurePosixPath, str],
    bucket: PurePosixPath,
    s3_prefix: PurePosixPath,
) -> bool:
    """Check if any file for an accession needs updating by comparing FTP checksums to S3 metadata.

    :param target_checksums: dict mapping file names to expected MD5 checksums from FTP
    :param bucket: S3 bucket name
    :param s3_prefix: S3 key prefix for the assembly in question
    :return: True if any file is missing or has a checksum mismatch, False otherwise
    """
    for fname, expected_md5 in target_checksums.items():
        s3_path = bucket / s3_prefix / fname
        s3_md5 = ""
        try:
            obj_info = head_object(str(s3_path))
            s3_md5 = obj_info.get("Metadata", {}).get("md5", "")
        except ClientError as e:
            code = e.response.get("Error", {}).get("Code")
            # boto3 error codes are not always ints, so check the conversion first
            try:
                code_int = int(code)
            except (TypeError, ValueError):
                code_int = None

            if code_int == HTTPStatus.NOT_FOUND:
                logger.debug("File missing from store: %s", s3_path)
                return True
            raise

        if s3_md5 != expected_md5:
            logger.debug("MD5 mismatch for %s: S3=%s FTP=%s", s3_path, s3_md5, expected_md5)
            return True

    return False


def verify_transfer_candidates(  # noqa: PLR0913
    accessions: list[str],
    current_assemblies: dict[str, AssemblyRecord],
    bucket: PurePosixPath,
    key_prefix: PurePosixPath,
    ftp_host: str = FTP_HOST,
    progress_callback: Callable[[int, int, str], None] | None = None,
) -> list[str]:
    """Verify which transfer candidates actually need downloading.

    For each accession, downloads ``md5checksums.txt`` from NCBI FTP and
    compares the checksums of filtered files against the ``md5`` user metadata
    on corresponding S3 objects.  Only accessions where at least one file
    differs or is missing from S3 are returned.

    This acts as a final gate before Phase 2: even if the summary diff flags an
    assembly, we skip it if every file in the store already matches.

    :param accessions: list of candidate accessions (new + updated from diff)
    :param current_assemblies: parsed current assembly summary
    :param bucket: S3 bucket name
    :param key_prefix: S3 key prefix for the Lakehouse dataset root
    :param ftp_host: NCBI FTP hostname
    :param progress_callback: optional callable invoked after each accession is
        processed with ``(done, total, accession)`` so callers can display a
        progress bar.  ``done`` is the 1-based count of completed accessions.
    :return: filtered list of accessions that actually need downloading
    """
    if not accessions:
        return []

    ftp: Any = None  # lazily connected only when needed
    confirmed: list[str] = []
    pruned = 0
    skipped_missing = 0
    last_activity = time.monotonic()

    def _progress(done: int, total: int, acc: str) -> None:
        if progress_callback is not None:
            progress_callback(done, total, acc)

    try:
        for done, acc in enumerate(accessions, start=1):
            rec = current_assemblies.get(acc)
            if not rec:
                confirmed.append(acc)
                _progress(done, len(accessions), acc)
                continue

            # Build S3 prefix for this assembly
            s3_rel = build_accession_path(rec.assembly_dir)
            s3_prefix = key_prefix / s3_rel

            # Quick check: does *anything* exist under this prefix?
            resp = list_objects(str(bucket / s3_prefix), max_keys=1)
            if not resp:
                # Nothing in the store — definitely needs downloading
                confirmed.append(acc)
                skipped_missing += 1
                _progress(done, len(accessions), acc)
                continue

            # Objects exist — need FTP md5 checksums to decide
            target_checksums, last_activity = _fetch_accession_checksums_from_ftp(
                ftp, ftp_host, last_activity, rec.ftp_url
            )

            if not target_checksums:
                confirmed.append(acc)
                _progress(done, len(accessions), acc)
                continue

            # Short-circuit: if any file differs or is missing, keep the assembly
            needs_update = _does_accession_need_update(target_checksums, bucket, s3_prefix)

            if needs_update:
                confirmed.append(acc)
            else:
                pruned += 1
                logger.debug("Pruned %s — all files match S3 checksums", acc)

            _progress(done, len(accessions), acc)
    finally:
        if ftp is not None:
            with contextlib.suppress(Exception):
                ftp.quit()

    logger.info(
        "Checksum verification: %d confirmed (%d missing from store), %d pruned (of %d candidates)",
        len(confirmed),
        skipped_missing,
        pruned,
        len(accessions),
    )
    return confirmed


# Manifest writing


def write_transfer_manifest(
    diff: DiffResult,
    current_assemblies: dict[str, AssemblyRecord],
    output_path: Path,
    ftp_host: str = FTP_HOST,
) -> list[PurePosixPath]:
    """Write the transfer manifest (new + updated assemblies).

    Each line is an FTP directory path suitable for Phase 2 download.

    :param diff: computed diff result
    :param current_assemblies: parsed current assembly summary
    :param output_path: path to write the manifest file
    :param ftp_host: FTP hostname for URL stripping
    :return: list of FTP paths written
    """
    to_transfer = diff.new + diff.updated
    paths: list[PurePosixPath] = []
    for acc in sorted(to_transfer):
        rec = current_assemblies.get(acc)
        if not rec:
            continue
        ftp_path = _ftp_dir_from_url(rec.ftp_url, ftp_host)
        paths.append(ftp_path)

    with Path(output_path).open("w") as f:
        f.writelines(f"{p}\n" for p in paths)

    logger.info("Wrote %d entries to transfer manifest: %s", len(paths), output_path)
    return paths


def write_removed_manifest(diff: DiffResult, output_path: Path) -> list[str]:
    """Write the removed manifest (replaced + suppressed accessions).

    :param diff: computed diff result
    :param output_path: path to write the manifest file
    :return: list of accessions written
    """
    removed = sorted(diff.replaced + diff.suppressed)
    with output_path.open("w") as f:
        f.writelines(acc + "\n" for acc in removed)
    logger.info("Wrote %d entries to removed manifest: %s", len(removed), output_path)
    return removed


def write_updated_manifest(diff: DiffResult, output_path: Path) -> list[str]:
    """Write the updated manifest (accessions whose content changed).

    This file is consumed by Phase 3 to archive existing S3 objects
    before they are overwritten by the new versions.

    :param diff: computed diff result
    :param output_path: path to write the manifest file
    :return: list of accessions written
    """
    updated = sorted(diff.updated)
    with output_path.open("w") as f:
        f.writelines(acc + "\n" for acc in updated)
    logger.info("Wrote %d entries to updated manifest: %s", len(updated), output_path)
    return updated


def write_diff_summary(
    diff: DiffResult,
    output_path: Path,
    database: str,
    prefix_from: str | None = None,
    prefix_to: str | None = None,
) -> dict[str, Any]:
    """Write a JSON diff summary file.

    :param diff: computed diff result
    :param output_path: path to write the JSON file
    :param database: database name (``"refseq"`` or ``"genbank"``)
    :param prefix_from: lower bound of prefix filter (if any)
    :param prefix_to: upper bound of prefix filter (if any)
    :return: the summary dict that was written
    """
    summary: dict[str, Any] = {
        "database": database,
        "timestamp": datetime.now(UTC).isoformat(),
        "prefix_range": {
            "from": prefix_from,
            "to": prefix_to,
        },
        "counts": {
            "new": len(diff.new),
            "updated": len(diff.updated),
            "replaced": len(diff.replaced),
            "suppressed": len(diff.suppressed),
            "total_to_transfer": len(diff.new) + len(diff.updated),
            "total_to_remove": len(diff.replaced) + len(diff.suppressed),
        },
        "accessions": {
            "new": diff.new,
            "updated": diff.updated,
            "replaced": diff.replaced,
            "suppressed": diff.suppressed,
        },
    }
    with output_path.open("w") as f:
        json.dump(summary, f, indent=2)
    logger.info("Wrote diff summary to: %s", output_path)
    return summary
