"""NCBI FTP assembly-specific domain logic.

Provides path helpers, file filters, MD5 checksum parsing, and single-assembly
download logic for NCBI GenBank/RefSeq assemblies.  Orchestration (batching,
threading, CLI) lives in :mod:`cdm_data_loaders.pipelines.ncbi_ftp_download`.
"""

import contextlib
import time
from ftplib import FTP
from logging import Logger, getLogger
from pathlib import Path, PurePosixPath
from typing import Any

from cdm_data_loaders.ncbi_ftp.constants import ACCESSION_PARTS_REGEX, ASSEMBLY_PATH_REGEX, FTP_HOST
from cdm_data_loaders.utils.checksums import compute_file_checksum
from cdm_data_loaders.utils.file_transfer.ftp_client import connect_ftp, ftp_noop_keepalive, ftp_retrieve_text

logger: Logger = getLogger(__name__)

FILE_FILTERS = [
    "_gene_ontology.gaf.gz",
    "_genomic.fna.gz",
    "_genomic.gff.gz",
    "_protein.faa.gz",
    "_ani_contam_ranges.tsv",
    "_assembly_regions.txt",
    "_assembly_report.txt",
    "_assembly_stats.txt",
    "_gene_expression_counts.txt.gz",
    "_normalized_gene_expression_counts.txt.gz",
]


def parse_md5_checksums_file(text: str) -> dict[PurePosixPath, str]:
    """Parse an NCBI ``md5checksums.txt`` file into a filename-to-hash mapping.

    Each line has the format ``<md5>  ./<filename>`` (two-space separator).

    :param text: raw text of the md5checksums.txt file
    :return: dict mapping filename to MD5 hex digest
    """
    checksums: dict[PurePosixPath, str] = {}
    for raw_line in text.strip().splitlines():
        parts = raw_line.strip().split("  ", maxsplit=1)
        if len(parts) == 2:  # noqa: PLR2004
            md5_hash, filename = parts
            checksums[PurePosixPath(filename.removeprefix("./"))] = md5_hash.strip()
    return checksums


# Path helpers


def build_accession_path(assembly_dir: PurePosixPath) -> PurePosixPath:
    """Build the relative output path for an assembly directory.

    Produces ``raw_data/{GCF|GCA}/{000}/{001}/{215}/{assembly_dir}/``.

    :param assembly_dir: full assembly directory name (e.g. ``GCF_000001215.4_Release_6...``)
    :return: relative path
    :raises ValueError: if the assembly directory name cannot be parsed
    """
    if m := ACCESSION_PARTS_REGEX.match(f"{assembly_dir}"):  # returns e.g., "GCF", "000", "001", "215" captured groups
        return PurePosixPath("raw_data", *m.groups(), assembly_dir)
    msg = f"Cannot parse accession: {assembly_dir}"
    raise ValueError(msg)


def parse_assembly_path(assembly_path: PurePosixPath) -> tuple[str, PurePosixPath, str]:
    """Extract database, assembly_dir, and accession from an FTP assembly path.

    :param assembly_path: FTP directory path (e.g. ``/genomes/all/GCF/000/.../GCF_000001215.4_Rel.../``)
    :return: tuple of ``(database, assembly_dir, accession)``
    :raises ValueError: if the path cannot be parsed
    """
    if m := ASSEMBLY_PATH_REGEX.search(str(assembly_path)):
        return m.group(3), PurePosixPath(m.group(1)), m.group(2)
    msg = f"Cannot parse assembly path: {assembly_path}"
    raise ValueError(msg)


# Single assembly download


def _download_and_verify(  # noqa: PLR0913
    ftp: FTP,
    filename: PurePosixPath,
    dest_dir: Path,
    md5_checksums: dict[PurePosixPath, str],
    stats: dict[str, Any],
    last_activity: float,
) -> float:
    """Download one file, verify its MD5, and write a sidecar if valid."""
    last_activity = ftp_noop_keepalive(ftp, last_activity)
    local_file = dest_dir / filename
    expected_md5 = md5_checksums.get(filename)

    # Define a helper function to validate the MD5 checksum of the downloaded file if one exists.
    # Returns False if file checksums mismatch, True otherwise
    def validate_file(local_file: Path) -> bool:
        if expected_md5:
            if compute_file_checksum(local_file, "md5") == expected_md5:
                (dest_dir / f"{filename}.md5").write_text(expected_md5)
            else:
                logger.warning(
                    "  MD5 mismatch for %s: expected %s",
                    filename,
                    expected_md5,
                )
                return False
            logger.debug("  MD5 verified: %s", filename)
        else:
            stats["files_without_checksum"] += 1
        return True

    # Try 3 times to download and validate the file; if it fails, delete any partial file and return
    for attempt in range(1, 4):
        logger.debug("  Downloading %s (attempt %d/3)", filename, attempt)
        with local_file.open("wb") as f:
            ftp.retrbinary(f"RETR {filename}", f.write)
        last_activity = time.monotonic()

        if validate_file(local_file):
            stats["files_downloaded"] += 1
            return last_activity

    stats["files_skipped_checksum_mismatch"] += 1
    local_file.unlink(missing_ok=True)
    return last_activity


def download_assembly_to_local(
    assembly_path: PurePosixPath,
    output_dir: Path,
    ftp_host: str = FTP_HOST,
    ftp: FTP | None = None,
) -> dict[str, Any]:
    """Download one assembly from NCBI FTP to a local directory.

    Creates a directory structure under *output_dir* matching the S3 layout,
    downloads filtered files, verifies MD5 checksums, and writes ``.md5``
    sidecar files for downstream metadata.

    :param assembly_path: FTP directory path for the assembly
    :param output_dir: base output directory
    :param ftp_host: FTP hostname
    :param ftp: optional existing FTP connection (caller manages lifecycle)
    :return: json parsable dict with download statistics
    """
    _database, assembly_dir, accession = parse_assembly_path(assembly_path)
    rel_path = build_accession_path(assembly_dir)
    dest_dir: Path = output_dir / rel_path
    dest_dir.mkdir(parents=True, exist_ok=True)

    logger.debug("Downloading %s -> %s", accession, dest_dir)

    owns_ftp = ftp is None
    if owns_ftp:
        ftp = connect_ftp(ftp_host)
    stats: dict[str, Any] = {
        "accession": accession,
        "assembly_dir": str(assembly_dir),
        "files_downloaded": 0,
        "files_skipped_checksum_mismatch": 0,
        "files_without_checksum": 0,
    }

    try:
        ftp.cwd(str(assembly_path))

        file_strings: list[str] = []
        ftp.retrlines("NLST", file_strings.append)
        files: list[PurePosixPath] = [PurePosixPath(f) for f in file_strings]

        # Download and parse md5checksums.txt
        md5_checksums: dict[PurePosixPath, str] = {}
        if PurePosixPath("md5checksums.txt") in files:
            md5_text = ftp_retrieve_text(ftp, PurePosixPath("md5checksums.txt"))
            md5_checksums = parse_md5_checksums_file(md5_text)
            (dest_dir / "md5checksums.txt").write_text(md5_text)
            stats["files_downloaded"] += 1

        target_files = [f for f in files if any(f.match(f"**{s}") for s in FILE_FILTERS)]
        last_activity = time.monotonic()

        for filename in target_files:
            last_activity = _download_and_verify(ftp, filename, dest_dir, md5_checksums, stats, last_activity)

        logger.debug("  %s: %d files downloaded", accession, stats["files_downloaded"])

    finally:
        if owns_ftp:
            with contextlib.suppress(Exception):
                ftp.quit()

    return stats
