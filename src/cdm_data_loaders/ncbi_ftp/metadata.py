"""Frictionless data package descriptor creation for NCBI FTP assemblies.

Creates KBase credit metadata descriptors for each promoted assembly,
matching the schema produced by ``kbase-transfers/scripts/ncbi/download_genomes.py``.

Each descriptor is a frictionless ``Package``-compatible JSON document
describing the assembly's data files, stored at::

    {key_prefix}metadata/{assembly_dir}_datapackage.json

and archived alongside raw data at::

    {key_prefix}archive/{release_tag}/metadata/{assembly_dir}_datapackage.json

The descriptor ``resources`` list records the final Lakehouse S3 key, byte
size, file format, and MD5 hash of each promoted data file.
"""

import json
import tempfile
from datetime import UTC, datetime
from logging import Logger, getLogger
from pathlib import Path, PurePosixPath
from typing import Any, TypedDict

from frictionless import Package

from cdm_data_loaders.utils.file_transfer.s3 import client
from cdm_data_loaders.utils.file_transfer.s3.object_utils import copy_object

logger: Logger = getLogger(__name__)

_NCBI_CONTRIBUTOR = {
    "contributor_type": "Organization",
    "name": "National Center for Biotechnology Information",
    "contributor_id": "ROR:02meqm098",
    "contributor_roles": "DataCurator",
}
_NCBI_PUBLISHER = {
    "organization_name": "National Center for Biotechnology Information",
    "organization_id": "ROR:02meqm098",
}
_SAVED_BY = "cdm-data-loaders-ncbi-ftp"
_SCHEMA_VERSION = "1.0"


class DescriptorResource(TypedDict, total=False):
    """A single resource entry in the frictionless descriptor ``resources`` list."""

    name: str
    path: PurePosixPath
    format: str
    bytes: int | None
    hash: str | None


# Public helpers


def build_descriptor_key(assembly_dir: PurePosixPath, key_prefix: PurePosixPath) -> PurePosixPath:
    """Return the S3 key for the live descriptor of *assembly_dir*.

    :param assembly_dir: full assembly directory name, e.g. ``GCF_000001215.4_Release_6_plus_ISO1_MT``
    :param key_prefix: Lakehouse key prefix (trailing slash optional)
    :return: S3 key, e.g. ``tenant-general-warehouse/.../ncbi/metadata/GCF_..._datapackage.json``
    """
    return key_prefix / "metadata" / assembly_dir.with_name(f"{assembly_dir.name}_datapackage.json")


def build_archive_descriptor_key(
    assembly_dir: PurePosixPath, release_tag: str, key_prefix: PurePosixPath, archive_reason: str = "unknown"
) -> PurePosixPath:
    """Return the S3 key for the archived descriptor of *assembly_dir*.

    :param assembly_dir: full assembly directory name
    :param release_tag: NCBI release tag used in the archive path, e.g. ``"2024-01"``
    :param key_prefix: Lakehouse key prefix
    :param archive_reason: reason for archival, encoded as a path segment
    :return: S3 key under ``archive/{release_tag}/{archive_reason}/metadata/``
    """
    return (
        key_prefix
        / "archive"
        / release_tag
        / archive_reason
        / "metadata"
        / assembly_dir.with_name(f"{assembly_dir.name}_datapackage.json")
    )


def create_descriptor(
    assembly_dir: PurePosixPath,
    accession_full: str,
    resources: list[DescriptorResource],
    *,
    timestamp: int | None = None,
) -> dict[str, Any]:
    """Build a KBase credit metadata descriptor for an NCBI assembly.

    Matches the schema produced by
    ``kbase-transfers/scripts/ncbi/download_genomes.py::create_frictionless_descriptor()``.

    Resource names are lowercased.  Resources whose ``hash`` value is ``None``
    have the ``hash`` key removed entirely (frictionless does not accept null
    hash values).

    :param assembly_dir: full assembly directory name (includes the accession
        suffix, e.g. ``GCF_000001215.4_Release_6_plus_ISO1_MT``)
    :param accession_full: accession without suffix, e.g. ``GCF_000001215.4``
    :param resources: list of :class:`DescriptorResource` dicts
    :param timestamp: Unix timestamp to embed; defaults to ``datetime.now(UTC)``
    :return: descriptor dict ready for serialisation and frictionless validation
    """
    ts = timestamp if timestamp is not None else int(datetime.now(UTC).timestamp())
    version = accession_full.rsplit(".", 1)[-1]  # e.g. "4" from "GCF_000001215.4"

    # Normalise resources: lowercase name, drop null hash
    normalised: list[dict[str, Any]] = []
    for res in resources:
        entry: dict[str, Any] = {
            "name": res["name"].lower() if "name" in res else None,
            "path": str(res["path"]) if "path" in res else None,
            "format": res.get("format", ""),
        }
        if entry["name"] is None:
            msg = f"Missing name for file descriptor in {accession_full}"
            raise ValueError(msg)
        if entry["path"] is None:
            msg = f"Missing path for file descriptor in {accession_full}"
            raise ValueError(msg)
        if res.get("bytes") is not None:
            entry["bytes"] = res.get("bytes")
        if res.get("hash") is not None:
            entry["hash"] = res.get("hash")
        normalised.append(entry)

    return {
        "identifier": f"NCBI:{accession_full}",
        "resource_type": "dataset",
        "version": version,
        "titles": [{"title": f"NCBI Genome Assembly {assembly_dir}"}],
        "descriptions": [
            {"description_text": (f"Genome assembly files for {accession_full} downloaded from NCBI Datasets")}
        ],
        "url": f"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession_full}/",
        "contributors": [_NCBI_CONTRIBUTOR],
        "publisher": _NCBI_PUBLISHER,
        "license": {},
        "meta": {
            "credit_metadata_schema_version": _SCHEMA_VERSION,
            "credit_metadata_source": [
                {
                    "source_name": "NCBI Genomes FTP",
                    "source_url": "ftp.ncbi.nlm.nih.gov/genomes/all/",
                    "access_timestamp": ts,
                }
            ],
            "saved_by": _SAVED_BY,
            "timestamp": ts,
        },
        "resources": normalised,
    }


def validate_descriptor(descriptor: dict[str, Any], accession_full: str) -> None:
    """Validate a descriptor with frictionless.

    :param descriptor: descriptor dict from :func:`create_descriptor`
    :param accession_full: accession (used only in error messages)
    :raises ValueError: if frictionless reports any metadata errors
    """
    errors = list(Package.metadata_validate(descriptor))
    if errors:
        error_details = "; ".join(str(e) for e in errors)
        msg = f"Frictionless validation failed for {accession_full}: {error_details}"
        raise ValueError(msg)
    logger.debug("Frictionless descriptor valid for %s", accession_full)


def upload_descriptor(
    descriptor: dict[str, Any],
    assembly_dir: PurePosixPath,
    bucket: PurePosixPath,
    key_prefix: PurePosixPath,
    *,
    dry_run: bool = False,
) -> PurePosixPath:
    """Serialise and upload a descriptor to the live ``metadata/`` path.

    :param descriptor: descriptor dict from :func:`create_descriptor`
    :param assembly_dir: full assembly directory name
    :param bucket: S3 bucket name
    :param key_prefix: Lakehouse key prefix
    :param dry_run: if True, log without uploading
    :return: S3 key the descriptor was (or would be) written to
    """
    key = build_descriptor_key(assembly_dir, key_prefix)

    if dry_run:
        logger.debug("[dry-run] would upload descriptor: s3://%s/%s", bucket, key)
        return key

    s3_client = client.get_s3_client()
    body = json.dumps(descriptor, indent=2).encode()

    with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as tmp:
        tmp_path = tmp.name
        tmp.write(body)

    try:
        s3_client.upload_file(Filename=tmp_path, Bucket=str(bucket), Key=str(key))
        logger.debug("Uploaded descriptor: s3://%s/%s", bucket, key)
    finally:
        Path(tmp_path).unlink()

    return key


def archive_descriptor(  # noqa: PLR0913
    assembly_dir: PurePosixPath,
    bucket: PurePosixPath,
    key_prefix: PurePosixPath,
    release_tag: str,
    *,
    archive_reason: str = "unknown",
    dry_run: bool = False,
) -> bool:
    """Copy the live descriptor to the archive path.

    If the live descriptor does not yet exist (e.g. archival is triggered
    before the first promote), logs a warning and returns ``False``.

    :param assembly_dir: full assembly directory name
    :param bucket: S3 bucket name
    :param key_prefix: Lakehouse key prefix
    :param release_tag: NCBI release tag for the archive path
    :param archive_reason: metadata value describing why archived (matches raw data metadata)
    :param dry_run: if True, log without copying
    :return: ``True`` if the descriptor was (or would be) archived; ``False`` if not found
    """
    source_key = build_descriptor_key(assembly_dir, key_prefix)
    archive_key = build_archive_descriptor_key(assembly_dir, release_tag, key_prefix, archive_reason)

    if dry_run:
        logger.debug("[dry-run] would archive descriptor: s3://%s/%s -> %s", bucket, source_key, archive_key)
        return True

    s3_client = client.get_s3_client()
    try:
        s3_client.head_object(Bucket=str(bucket), Key=str(source_key))
    except s3_client.exceptions.NoSuchKey:
        logger.warning("Descriptor not found, skipping archive: s3://%s/%s", bucket, source_key)
        return False
    except Exception as e:
        # head_object raises ClientError with 404 when key is absent
        if hasattr(e, "response") and e.response.get("Error", {}).get("Code") in ("404", "NoSuchKey"):  # type: ignore[union-attr]
            logger.warning("Descriptor not found, skipping archive: s3://%s/%s", bucket, source_key)
            return False
        raise

    _ = copy_object(
        str(bucket / source_key),
        str(bucket / archive_key),
    )
    logger.debug("Archived descriptor: s3://%s/%s -> %s", bucket, source_key, archive_key)
    return True
