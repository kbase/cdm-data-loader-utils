"""Tests for ncbi_ftp.promote module — S3 promote, archive, manifest trimming."""

import hashlib
import logging
from http import HTTPStatus
from pathlib import Path, PurePosixPath
from typing import Protocol, cast
from unittest.mock import patch

import botocore.client
import pytest

from cdm_data_loaders.ncbi_ftp.promote import (
    _archive_assemblies,
    _archive_objects,
    _dry_run_output,
    _get_accession_path_prefix,
    _get_source_dest_pairs_for_accession,
    _trim_manifest,
    promote_from_s3,
)
from tests.cdm_data_loaders.ncbi_ftp.conftest import ACC_PATH_215, ACC_PATH_845, TEST_BUCKET

DEFAULT_LAKEHOUSE_KEY_PREFIX: PurePosixPath = PurePosixPath("tenant-general-warehouse/kbase/datasets/ncbi")

# Promotion test constants

_STAGE_PREFIX: PurePosixPath = PurePosixPath("staging") / "run1"

# Assembly 1
_ACC1: str = "GCF_000001215.4"
_STG1: PurePosixPath = _STAGE_PREFIX / "raw_data" / ACC_PATH_215
_LKH1: PurePosixPath = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215

# Assembly 2
_ACC2: str = "GCF_000005845.2"
_STG2: PurePosixPath = _STAGE_PREFIX / "raw_data" / ACC_PATH_845
_LKH2: PurePosixPath = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_845


class DownloadFileClient(Protocol):
    """Protocol for dynamic BaseClient with download_file function."""

    def download_file(self, Bucket: str, Key: str, Filename: str, **kw: object) -> None:  #  noqa: N803
        """Download a file from an S3 store."""


def _stage(
    s3: botocore.client.BaseClient,
    staging_base: PurePosixPath,
    files: dict[PurePosixPath, bytes],
    *,
    with_md5: bool = True,
    with_crc64: bool = False,
) -> list[PurePosixPath]:
    """Stage files at *staging_base*, optionally adding .md5 / .crc64nvme sidecars.

    Returns list of all staged keys (data files only, not sidecars).
    """
    keys: list[PurePosixPath] = []
    for fname, content in files.items():
        key: PurePosixPath = staging_base / fname
        s3.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=content)
        keys.append(key)
        if with_md5:
            s3.put_object(
                Bucket=str(TEST_BUCKET),
                Key=str(key.with_name(f"{key.name}.md5")),
                Body=hashlib.md5(content).hexdigest().encode(),  # noqa: S324
            )
        if with_crc64:
            s3.put_object(Bucket=str(TEST_BUCKET), Key=str(key.with_name(f"{key.name}.crc64nvme")), Body=b"fake-crc")
    return keys


def _stage_files(s3_client: botocore.client.BaseClient, prefix: PurePosixPath) -> None:
    """Upload sample staged files to mock S3."""
    for key in [
        prefix / "raw_data" / ACC_PATH_215 / "GCF_000001215.4_genomic.fna.gz",
        prefix / "raw_data" / ACC_PATH_215 / "GCF_000001215.4_genomic.fna.gz.md5",
        prefix / "download_report.json",
    ]:
        body = b"md5hash123" if key.match("**.md5") else b"data"
        s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=body)


@pytest.mark.s3
def test_promote_dry_run_no_writes(mock_s3_client: botocore.client.BaseClient) -> None:
    """Verify dry_run does not write any objects."""
    prefix = PurePosixPath("staging/run1")
    _stage_files(mock_s3_client, prefix)

    report = promote_from_s3(
        staging_key_prefix=prefix,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        dry_run=True,
    )
    assert report["promoted"] == 1
    assert report["dry_run"] is True

    final_key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215 / "GCF_000001215.4_genomic.fna.gz"
    assert mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(final_key)).get("KeyCount", 0) == 0


@pytest.mark.s3
def test_promote_with_metadata(mock_s3_client: botocore.client.BaseClient) -> None:
    """Objects are promoted with MD5 metadata; download_report.json is skipped."""
    prefix = PurePosixPath("staging/run1")
    _stage_files(mock_s3_client, prefix)

    report = promote_from_s3(
        staging_key_prefix=prefix,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )
    assert report["promoted"] == 1  # only .fna.gz, not download_report.json
    assert report["failed"] == 0

    final_key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215 / "GCF_000001215.4_genomic.fna.gz"
    resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(final_key))
    assert resp["Metadata"].get("md5") == "md5hash123"


@pytest.mark.s3
@pytest.mark.parametrize(
    ("manifest_body", "promoted_set", "expected_present", "expected_absent"),
    [
        pytest.param(
            "/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6/\n"
            "/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/\n",
            {"GCF_000001215.4"},
            ["GCF_000001405.40"],
            ["GCF_000001215.4"],
            id="partial",
        ),
        pytest.param(
            "/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6/\n",
            {"GCF_000001215.4"},
            [],
            ["GCF_000001215.4"],
            id="all",
        ),
    ],
)
def test_trim_manifest(
    mock_s3_client: botocore.client.BaseClient,
    manifest_body: str,
    promoted_set: set[str],
    expected_present: list[str],
    expected_absent: list[str],
) -> None:
    """Promoted accessions are removed; others remain (partial) or the manifest empties (all)."""
    manifest_key = PurePosixPath("manifests/transfer_manifest.txt")
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(manifest_key), Body=manifest_body.encode())
    _trim_manifest(manifest_key, TEST_BUCKET, promoted_set)
    remaining = mock_s3_client.get_object(Bucket=str(TEST_BUCKET), Key=str(manifest_key))["Body"].read().decode()
    for acc in expected_present:
        assert acc in remaining
    for acc in expected_absent:
        assert acc not in remaining


# Helpers for _archive_assemblies tests


def _mock_list_objects(path: str) -> list[dict[str, str]]:
    bucket = PurePosixPath("some-bucket")
    prefix = PurePosixPath("some/prefix")
    p = PurePosixPath(path)
    if p == bucket / prefix / "raw_data" / ACC_PATH_215:
        return [
            {"Key": f"{prefix / 'raw_data' / ACC_PATH_215 / 'GCF_000001215.4_genomic.fna.gz'}"},
            {"Key": f"{prefix / 'raw_data' / ACC_PATH_215 / 'GCF_000001215.4_protein.faa.gz'}"},
        ]
    if p == bucket / prefix / "raw_data" / ACC_PATH_845:
        return [
            {"Key": f"{prefix / 'raw_data' / ACC_PATH_845 / 'GCF_000005845.2_genomic.fna.gz'}"},
        ]
    return []


@pytest.mark.parametrize(
    ("accession", "prefix", "expected"),
    [
        pytest.param(
            "GCF_012345678.90_Some_description",
            PurePosixPath("some") / "prefix",
            PurePosixPath("some")
            / "prefix"
            / "raw_data"
            / "GCF"
            / "012"
            / "345"
            / "678"
            / "GCF_012345678.90_Some_description",
            id="standard",
        ),
        pytest.param(
            "GCF_000001215.4_Release_6_plus_ISO1_MT",
            PurePosixPath("another") / "prefix",
            PurePosixPath("another") / "prefix" / "raw_data" / ACC_PATH_215,
            id="standard-2",
        ),
        pytest.param("INVALID_ACCESSION", PurePosixPath("prefix"), None, id="invalid-format"),
    ],
)
def test_get_accession_path_prefix(
    accession: str,
    prefix: PurePosixPath,
    expected: PurePosixPath | None,
) -> None:
    """get_accession_path_prefix returns correct path for valid accessions, None for invalid."""
    result = _get_accession_path_prefix(accession, prefix)
    assert result == expected


@pytest.mark.parametrize(
    ("accession", "bucket", "prefix", "release_tag", "archive_reason", "expected"),
    [
        pytest.param(
            "GCF_000001215.4_Release_6_plus_ISO1_MT",
            "some-bucket",
            PurePosixPath("some/prefix"),
            "2024-01",
            "test-reason",
            [
                (
                    PurePosixPath("some") / "prefix" / "raw_data" / ACC_PATH_215 / "GCF_000001215.4_genomic.fna.gz",
                    PurePosixPath("some")
                    / "prefix"
                    / "archive"
                    / "2024-01"
                    / "test-reason"
                    / "raw_data"
                    / ACC_PATH_215
                    / "GCF_000001215.4_genomic.fna.gz",
                ),
                (
                    PurePosixPath("some") / "prefix" / "raw_data" / ACC_PATH_215 / "GCF_000001215.4_protein.faa.gz",
                    PurePosixPath("some")
                    / "prefix"
                    / "archive"
                    / "2024-01"
                    / "test-reason"
                    / "raw_data"
                    / ACC_PATH_215
                    / "GCF_000001215.4_protein.faa.gz",
                ),
            ],
            id="standard",
        ),
        pytest.param(
            "GCF_000005845.2_ASM584v2",
            "some-bucket",
            PurePosixPath("some/prefix"),
            "2024-01",
            "test-reason",
            [
                (
                    PurePosixPath("some") / "prefix" / "raw_data" / ACC_PATH_845 / "GCF_000005845.2_genomic.fna.gz",
                    PurePosixPath("some")
                    / "prefix"
                    / "archive"
                    / "2024-01"
                    / "test-reason"
                    / "raw_data"
                    / ACC_PATH_845
                    / "GCF_000005845.2_genomic.fna.gz",
                ),
            ],
            id="standard-2",
        ),
        pytest.param(
            "GCF_000001405.39_GRCh38.p14",
            "some-bucket",
            PurePosixPath("some") / "prefix",
            "2024-01",
            "test-reason",
            [],
            id="accession-not-found",
        ),
        pytest.param(
            "INVALID_ACCESSION",
            "some-bucket",
            PurePosixPath("some") / "prefix",
            "2024-01",
            "test-reason",
            [],
            id="invalid-format",
        ),
    ],
)
def test_get_source_dest_pairs_for_accession(
    accession: str,
    bucket: PurePosixPath,
    prefix: PurePosixPath,
    release_tag: str,
    archive_reason: str,
    expected: list[tuple[str, str]],
) -> None:
    """get_source_dest_pairs_for_accession returns correct source-dest pairs for valid accessions, empty list for invalid."""
    with patch("cdm_data_loaders.ncbi_ftp.promote.list_objects", side_effect=_mock_list_objects):
        result = _get_source_dest_pairs_for_accession(
            accession,
            bucket,
            prefix,
            release_tag,
            archive_reason,
        )
        assert result == expected


@pytest.mark.parametrize(
    ("key_pairs", "log_count", "expected", "info_log_strings"),
    [
        pytest.param(
            [(PurePosixPath("source1"), PurePosixPath("dest1")), (PurePosixPath("source2"), PurePosixPath("dest2"))],
            0,
            2,
            ["[dry-run] would archive: source1 -> dest1", "[dry-run] would archive: source2 -> dest2"],
            id="standard",
        ),
        pytest.param(
            [
                (PurePosixPath("source1"), PurePosixPath("dest1")),
                (PurePosixPath("source2"), PurePosixPath("dest2")),
                (PurePosixPath("source3"), PurePosixPath("dest3")),
            ],
            9,
            12,
            ["[dry-run] would archive: source1 -> dest1"],
            id="exceeds-log-cutoff",
        ),
        pytest.param(
            [],
            0,
            0,
            [],
            id="empty",
        ),
    ],
)
def test_dry_run_output(
    key_pairs: list[tuple[PurePosixPath, PurePosixPath]],
    log_count: int,
    expected: int,
    info_log_strings: list[str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    """In dry_run mode, _archive_assemblies logs source-dest pairs but does not copy."""
    with caplog.at_level(logging.INFO):
        result = _dry_run_output(key_pairs, log_count)
        assert result == expected
        for log_string in info_log_strings:
            assert log_string in caplog.text


@pytest.mark.parametrize(
    ("key_pairs", "bucket", "delete_source", "expected", "existing_objects"),
    [
        pytest.param(
            [(PurePosixPath("source1"), PurePosixPath("dest1")), (PurePosixPath("source2"), PurePosixPath("dest2"))],
            TEST_BUCKET,
            False,
            2,
            ["source1", "source2", "dest1", "dest2"],
            id="keep-source",
        ),
        pytest.param(
            [(PurePosixPath("source1"), PurePosixPath("dest1")), (PurePosixPath("source2"), PurePosixPath("dest2"))],
            TEST_BUCKET,
            True,
            2,
            ["dest1", "dest2"],
            id="delete-source",
        ),
        pytest.param(
            [],
            TEST_BUCKET,
            True,
            0,
            [],
            id="empty",
        ),
    ],
)
@pytest.mark.s3
def test_archive_objects(
    mock_s3_client: botocore.client.BaseClient,
    key_pairs: list[tuple[PurePosixPath, PurePosixPath]],
    bucket: PurePosixPath,
    delete_source: bool,
    expected: int,
    existing_objects: list[str],
) -> None:
    """_archive_assemblies copies source to dest for each pair, and deletes source if delete_source=True."""
    for source, _ in key_pairs:
        mock_s3_client.put_object(Bucket=str(bucket), Key=str(source), Body=b"data")
    assert _archive_objects(key_pairs, bucket, delete_source=delete_source) == expected
    assert mock_s3_client.list_objects_v2(Bucket=str(bucket)).get("KeyCount", 0) == len(existing_objects)
    for obj in mock_s3_client.list_objects_v2(Bucket=str(bucket)).get("Contents", []):
        assert obj["Key"] in existing_objects, f"Unexpected object in bucket: {obj['Key']}"


@pytest.mark.s3
def test_archive_assemblies_removed(mock_s3_client: botocore.client.BaseClient, tmp_path: Path) -> None:
    """Removed accessions are archived and originals deleted."""
    accession = "GCF_000005845.2"
    key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_845 / f"{accession}_genomic.fna.gz"
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"data")

    manifest = tmp_path / "removed.txt"
    manifest.write_text(f"{accession}\n")

    assert (
        _archive_assemblies(
            manifest,
            lakehouse_bucket=TEST_BUCKET,
            lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
            ncbi_release="2024-01",
            archive_reason="replaced_or_suppressed",
            delete_source=True,
        )
        == 1
    )
    assert mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(key)).get("KeyCount", 0) == 0

    archive_key = (
        DEFAULT_LAKEHOUSE_KEY_PREFIX
        / "archive"
        / "2024-01"
        / "replaced_or_suppressed"
        / "raw_data"
        / ACC_PATH_845
        / f"{accession}_genomic.fna.gz"
    )
    assert mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(archive_key)).get("KeyCount", 0) == 1


@pytest.mark.s3
def test_archive_assemblies_updated_no_delete(mock_s3_client: botocore.client.BaseClient, tmp_path: Path) -> None:
    """Updated accessions are archived but originals remain."""
    accession = "GCF_000001215.4"
    key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215 / f"{accession}_genomic.fna.gz"
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"original-data")

    manifest = tmp_path / "updated.txt"
    manifest.write_text(f"{accession}\n")

    assert (
        _archive_assemblies(
            manifest,
            lakehouse_bucket=TEST_BUCKET,
            lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
            ncbi_release="2024-06",
            archive_reason="updated",
            delete_source=False,
        )
        == 1
    )
    assert mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(key)).get("KeyCount", 0) == 1

    archive_key = (
        DEFAULT_LAKEHOUSE_KEY_PREFIX
        / "archive"
        / "2024-06"
        / "updated"
        / "raw_data"
        / ACC_PATH_215
        / f"{accession}_genomic.fna.gz"
    )
    resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(archive_key))
    assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_archive_assemblies_multiple_releases_no_collision(
    mock_s3_client: botocore.client.BaseClient, tmp_path: Path
) -> None:
    """Archiving the same accession in different releases creates distinct folders."""
    accession = "GCF_000001215.4"
    key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215 / f"{accession}_genomic.fna.gz"
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"v1-data")

    manifest = tmp_path / "updated.txt"
    manifest.write_text(f"{accession}\n")

    _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
    )
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"v2-data")
    _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-06",
        archive_reason="updated",
    )

    archive_key_1 = (
        DEFAULT_LAKEHOUSE_KEY_PREFIX
        / "archive"
        / "2024-01"
        / "updated"
        / "raw_data"
        / ACC_PATH_215
        / f"{accession}_genomic.fna.gz"
    )
    archive_key_2 = (
        DEFAULT_LAKEHOUSE_KEY_PREFIX
        / "archive"
        / "2024-06"
        / "updated"
        / "raw_data"
        / ACC_PATH_215
        / f"{accession}_genomic.fna.gz"
    )
    assert mock_s3_client.get_object(Bucket=str(TEST_BUCKET), Key=str(archive_key_1))["Body"].read() == b"v1-data"
    assert mock_s3_client.get_object(Bucket=str(TEST_BUCKET), Key=str(archive_key_2))["Body"].read() == b"v2-data"


@pytest.mark.s3
def test_archive_assemblies_dry_run(mock_s3_client: botocore.client.BaseClient, tmp_path: Path) -> None:
    """dry_run does not copy or delete anything."""
    accession = "GCF_000005845.2"
    key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_845 / f"{accession}_genomic.fna.gz"
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"data")

    manifest = tmp_path / "removed.txt"
    manifest.write_text(f"{accession}\n")

    assert (
        _archive_assemblies(
            manifest,
            lakehouse_bucket=TEST_BUCKET,
            lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
            ncbi_release="2024-01",
            archive_reason="replaced_or_suppressed",
            delete_source=True,
            dry_run=True,
        )
        == 1
    )
    assert mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(key)).get("KeyCount", 0) == 1

    archive_prefix = DEFAULT_LAKEHOUSE_KEY_PREFIX / "archive" / "2024-01"
    assert mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(archive_prefix)).get("KeyCount", 0) == 0


@pytest.mark.s3
def test_archive_assemblies_no_objects_skips(
    mock_s3_client: botocore.client.BaseClient,  # noqa: ARG001
    tmp_path: Path,
) -> None:
    """Accessions with no existing S3 objects are silently skipped."""
    manifest = tmp_path / "updated.txt"
    manifest.write_text("GCF_000001215.4\n")
    assert (
        _archive_assemblies(
            manifest,
            lakehouse_bucket=TEST_BUCKET,
            lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
            ncbi_release="2024-01",
        )
        == 0
    )


@pytest.mark.s3
def test_archive_assemblies_unknown_release_fallback(
    mock_s3_client: botocore.client.BaseClient, tmp_path: Path
) -> None:
    """ncbi_release=None falls back to 'unknown' in the archive path."""
    accession = "GCF_000001215.4"
    key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215 / f"{accession}_genomic.fna.gz"
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"data")

    manifest = tmp_path / "updated.txt"
    manifest.write_text(f"{accession}\n")

    assert (
        _archive_assemblies(
            manifest, lakehouse_bucket=TEST_BUCKET, lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX, ncbi_release=None
        )
        == 1
    )

    archive_key = (
        DEFAULT_LAKEHOUSE_KEY_PREFIX
        / "archive"
        / "unknown"
        / "unknown"
        / "raw_data"
        / ACC_PATH_215
        / f"{accession}_genomic.fna.gz"
    )
    assert mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(archive_key)).get("KeyCount", 0) == 1


# Concurrent / multi-file archive (new behaviour)


@pytest.mark.s3
def test_archive_assemblies_multi_file_all_copied(mock_s3_client: botocore.client.BaseClient, tmp_path: Path) -> None:
    """All files for an accession are copied concurrently — none missed."""
    accession = "GCF_000001215.4"
    base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215
    file_names = [
        f"{accession}_genomic.fna.gz",
        f"{accession}_protein.faa.gz",
        f"{accession}_rna.fna.gz",
        f"{accession}_assembly_report.txt",
        f"{accession}_assembly_stats.txt",
    ]
    for fname in file_names:
        mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / fname), Body=fname.encode())

    manifest = tmp_path / "updated.txt"
    manifest.write_text(f"{accession}\n")

    archived = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
        delete_source=False,
    )

    assert archived == len(file_names)
    archive_base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "archive" / "2024-01" / "updated" / "raw_data" / ACC_PATH_215
    for fname in file_names:
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(archive_base / fname))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_archive_assemblies_multi_file_content_preserved(
    mock_s3_client: botocore.client.BaseClient, tmp_path: Path
) -> None:
    """Archive copies preserve byte-for-byte content of each file."""
    accession = "GCF_000001215.4"
    base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215
    files = {
        f"{accession}_genomic.fna.gz": b"\x1f\x8bGENOMIC",
        f"{accession}_protein.faa.gz": b"\x1f\x8bPROTEIN",
    }
    for fname, body in files.items():
        mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / fname), Body=body)

    manifest = tmp_path / "updated.txt"
    manifest.write_text(f"{accession}\n")

    _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
        delete_source=False,
    )

    archive_base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "archive" / "2024-01" / "updated" / "raw_data" / ACC_PATH_215
    for fname, original_body in files.items():
        obj = mock_s3_client.get_object(Bucket=str(TEST_BUCKET), Key=str(archive_base / fname))
        assert obj["Body"].read() == original_body, f"Content mismatch for {fname}"


@pytest.mark.s3
def test_archive_assemblies_multi_file_delete_all(mock_s3_client: botocore.client.BaseClient, tmp_path: Path) -> None:
    """Batch delete removes ALL source files when delete_source=True."""
    accession = "GCF_000005845.2"
    base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_845
    file_names = [
        f"{accession}_genomic.fna.gz",
        f"{accession}_protein.faa.gz",
        f"{accession}_assembly_report.txt",
    ]
    for fname in file_names:
        mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / fname), Body=b"data")

    manifest = tmp_path / "removed.txt"
    manifest.write_text(f"{accession}\n")

    archived = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-03",
        archive_reason="replaced_or_suppressed",
        delete_source=True,
    )

    assert archived == len(file_names)
    # All sources deleted
    for fname in file_names:
        result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(base / fname))
        assert result.get("KeyCount", 0) == 0, f"Source not deleted: {fname}"
    # All archives present
    archive_base = (
        DEFAULT_LAKEHOUSE_KEY_PREFIX / "archive" / "2024-03" / "replaced_or_suppressed" / "raw_data" / ACC_PATH_845
    )
    for fname in file_names:
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(archive_base / fname))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


# Partial-archive idempotency


@pytest.mark.s3
def test_archive_assemblies_partial_already_archived_overwritten(
    mock_s3_client: botocore.client.BaseClient, tmp_path: Path
) -> None:
    """Re-running archive after a partial run overwrites the already-archived files.

    Simulates a partial failure: file_a was archived, file_b was not.
    The second run should archive both file_a (overwrite) and file_b.
    """
    accession = "GCF_000001215.4"
    base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215
    archive_base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "archive" / "2024-01" / "updated" / "raw_data" / ACC_PATH_215

    file_a = f"{accession}_genomic.fna.gz"
    file_b = f"{accession}_protein.faa.gz"

    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / file_a), Body=b"new-genomic")
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / file_b), Body=b"new-protein")
    # Simulate partial prior run: file_a already archived with stale content
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=f"{archive_base / file_a}", Body=b"stale-genomic")

    manifest = tmp_path / "updated.txt"
    manifest.write_text(f"{accession}\n")

    archived = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
        delete_source=False,
    )

    assert archived == 2  # noqa: PLR2004
    # file_a should now have the current content (overwritten)
    obj_a = mock_s3_client.get_object(Bucket=str(TEST_BUCKET), Key=str(archive_base / file_a))
    assert obj_a["Body"].read() == b"new-genomic", "Re-run should overwrite stale archive"
    # file_b should now be archived
    obj_b = mock_s3_client.get_object(Bucket=str(TEST_BUCKET), Key=str(archive_base / file_b))
    assert obj_b["Body"].read() == b"new-protein"


@pytest.mark.s3
def test_archive_assemblies_partial_delete_resumes(mock_s3_client: botocore.client.BaseClient, tmp_path: Path) -> None:
    """Re-running replaced_or_suppressed archive after partial copy+delete is safe.

    Simulates: file_a was copied+deleted, file_b was copied but NOT deleted,
    file_c was not touched. The re-run finds only file_b and file_c present
    (file_a is gone), archives both, and deletes both.
    """
    accession = "GCF_000005845.2"
    base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_845
    archive_base = (
        DEFAULT_LAKEHOUSE_KEY_PREFIX / "archive" / "2024-03" / "replaced_or_suppressed" / "raw_data" / ACC_PATH_845
    )

    file_b = f"{accession}_protein.faa.gz"
    file_c = f"{accession}_assembly_report.txt"

    # file_a already gone (deleted in first partial run)
    # file_b present at source (not yet deleted from first partial run)
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / file_b), Body=b"protein")
    # file_c present at source (not touched at all)
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / file_c), Body=b"report")
    # file_a already at archive destination
    mock_s3_client.put_object(
        Bucket=str(TEST_BUCKET), Key=f"{archive_base / accession}_genomic.fna.gz", Body=b"genomic"
    )

    manifest = tmp_path / "removed.txt"
    manifest.write_text(f"{accession}\n")

    archived = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-03",
        archive_reason="replaced_or_suppressed",
        delete_source=True,
    )

    # Only the 2 remaining source files were archived
    assert archived == 2  # noqa: PLR2004
    # Both now gone from source
    for fname in (file_b, file_c):
        result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(base / fname))
        assert result.get("KeyCount", 0) == 0, f"Expected {fname} deleted"
    # file_a archive still intact (not touched by re-run)
    resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=f"{archive_base / accession}_genomic.fna.gz")
    assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_archive_assemblies_idempotent_updated_reruns_cleanly(
    mock_s3_client: botocore.client.BaseClient, tmp_path: Path
) -> None:
    """Running updated archive twice on the same data produces the same result."""
    accession = "GCF_000001215.4"
    base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215
    file_names = [f"{accession}_genomic.fna.gz", f"{accession}_protein.faa.gz"]
    for fname in file_names:
        mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / fname), Body=b"content")

    manifest = tmp_path / "updated.txt"
    manifest.write_text(f"{accession}\n")

    archived_1 = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
    )
    archived_2 = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
    )

    assert archived_1 == len(file_names)
    assert archived_2 == len(file_names)
    # Sources still present after both runs (delete_source=False)
    for fname in file_names:
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(base / fname))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_archive_assemblies_multi_accession_manifest(
    mock_s3_client: botocore.client.BaseClient, tmp_path: Path
) -> None:
    """Multiple accessions in a single manifest are all archived."""
    accessions = [
        ("GCF_000001215.4", "GCF_000001215.4_Release_6", "GCF/000/001/215"),
        ("GCF_000005845.2", "GCF_000005845.2_ASM584v2", "GCF/000/005/845"),
    ]
    for accession, asm_dir, path in accessions:
        key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / path / asm_dir / f"{accession}_genomic.fna.gz"
        mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"data")

    manifest = tmp_path / "updated.txt"
    manifest.write_text("\n".join(acc for acc, _, _ in accessions) + "\n")

    archived = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
    )

    assert archived == len(accessions)
    for accession, asm_dir, path in accessions:
        archive_key = (
            DEFAULT_LAKEHOUSE_KEY_PREFIX
            / "archive"
            / "2024-01"
            / "updated"
            / "raw_data"
            / path
            / asm_dir
            / f"{accession}_genomic.fna.gz"
        )
        result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(archive_key))
        assert result.get("KeyCount", 0) == 1, f"Archive missing for {accession}"


@pytest.mark.s3
def test_archive_assemblies_dry_run_multi_file(mock_s3_client: botocore.client.BaseClient, tmp_path: Path) -> None:
    """dry_run with multiple files per accession makes no copies and no deletes."""
    accession = "GCF_000005845.2"
    base = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_845
    file_names = [f"{accession}_genomic.fna.gz", f"{accession}_protein.faa.gz", f"{accession}_rna.fna.gz"]
    for fname in file_names:
        mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(base / fname), Body=b"data")

    manifest = tmp_path / "removed.txt"
    manifest.write_text(f"{accession}\n")

    archived = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="replaced_or_suppressed",
        delete_source=True,
        dry_run=True,
    )

    # Reported count matches
    assert archived == len(file_names)
    # No actual archive keys created
    archive_prefix = DEFAULT_LAKEHOUSE_KEY_PREFIX / "archive"
    result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(archive_prefix))
    assert result.get("KeyCount", 0) == 0
    # Sources untouched
    for fname in file_names:
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(base / fname))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_archive_assemblies_invalid_accession_skipped(
    mock_s3_client: botocore.client.BaseClient, tmp_path: Path
) -> None:
    """Malformed accession lines are skipped; valid ones still archived."""
    accession = "GCF_000001215.4"
    key = DEFAULT_LAKEHOUSE_KEY_PREFIX / "raw_data" / ACC_PATH_215 / f"{accession}_genomic.fna.gz"
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(key), Body=b"data")

    manifest = tmp_path / "mixed.txt"
    manifest.write_text("NOT_AN_ACCESSION\n\n   \n" + f"{accession}\n")

    archived = _archive_assemblies(
        manifest,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        ncbi_release="2024-01",
        archive_reason="updated",
    )
    assert archived == 1


# Concurrent / multi-file promotion (new behaviour)


@pytest.mark.s3
def test_promote_multi_file_all_land_at_final_path(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """All files for an assembly are promoted concurrently — none missed."""
    file_names = [
        f"{_ACC1}_genomic.fna.gz",
        f"{_ACC1}_protein.faa.gz",
        f"{_ACC1}_rna.fna.gz",
        f"{_ACC1}_assembly_report.txt",
        f"{_ACC1}_assembly_stats.txt",
    ]
    _stage(mock_s3_client, _STG1, {PurePosixPath(f): f.encode() for f in file_names})

    report = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    assert report["promoted"] == len(file_names)
    assert report["failed"] == 0
    for fname in file_names:
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=f"{_LKH1 / fname}")
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_promote_multi_file_content_preserved(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Content at the final key is byte-identical to the staged content."""
    files = {
        PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"\x1f\x8bGENOMIC",
        PurePosixPath(f"{_ACC1}_protein.faa.gz"): b"\x1f\x8bPROTEIN",
        PurePosixPath(f"{_ACC1}_rna.fna.gz"): b"\x1f\x8bRNA",
    }
    _stage(mock_s3_client, _STG1, files)

    promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    for fname, expected in files.items():
        obj = mock_s3_client.get_object(Bucket=str(TEST_BUCKET), Key=f"{_LKH1 / fname}")
        assert obj["Body"].read() == expected, f"Content mismatch for {fname}"


@pytest.mark.s3
def test_promote_md5_metadata_set_from_sidecar(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """MD5 metadata on the promoted object matches the .md5 sidecar value."""
    content = b"\x1f\x8bGENOMIC"
    fname = PurePosixPath(f"{_ACC1}_genomic.fna.gz")
    _stage(mock_s3_client, _STG1, {fname: content}, with_md5=True)
    expected_md5 = hashlib.md5(content).hexdigest()  # noqa: S324

    promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(_LKH1 / fname))
    assert resp["Metadata"].get("md5") == expected_md5


@pytest.mark.s3
def test_promote_no_sidecar_no_md5_metadata(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """A file staged without a .md5 sidecar is promoted but carries no md5 metadata."""
    fname = PurePosixPath(f"{_ACC1}_genomic.fna.gz")
    _stage(mock_s3_client, _STG1, {fname: b"data"}, with_md5=False)

    promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(_LKH1 / fname))
    assert resp["Metadata"].get("md5") is None


@pytest.mark.s3
def test_promote_staging_data_files_deleted_after_promote(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Staged data files are deleted from staging after a fully successful assembly promote."""
    files = {
        PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"genomic",
        PurePosixPath(f"{_ACC1}_protein.faa.gz"): b"protein",
    }
    staged_keys = _stage(mock_s3_client, _STG1, files)

    promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    for key in staged_keys:
        result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(key))
        assert result.get("KeyCount", 0) == 0, f"Staged data file not deleted: {key}"


@pytest.mark.s3
def test_promote_md5_sidecars_deleted_after_promote(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Staged .md5 sidecar files are deleted from staging after a successful promote."""
    files = {
        PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"genomic",
        PurePosixPath(f"{_ACC1}_protein.faa.gz"): b"protein",
    }
    staged_keys = _stage(mock_s3_client, _STG1, files, with_md5=True)

    promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    for key in staged_keys:
        for sidecar_key in (f"{key}.md5", f"{key}.crc64nvme"):
            result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(sidecar_key))
            assert result.get("KeyCount", 0) == 0, f"Sidecar not deleted: {sidecar_key}"


@pytest.mark.s3
def test_promote_crc64nvme_sidecars_deleted_after_promote(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Staged .crc64nvme sidecar files are also batch-deleted after a successful promote."""
    fname = PurePosixPath(f"{_ACC1}_genomic.fna.gz")
    _stage(mock_s3_client, _STG1, {fname: b"data"}, with_md5=True, with_crc64=True)

    promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    staged_key = f"{_STG1}{fname}"
    for sidecar_key in (f"{staged_key}.md5", f"{staged_key}.crc64nvme"):
        result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(sidecar_key))
        assert result.get("KeyCount", 0) == 0, f"Sidecar not deleted: {sidecar_key}"


@pytest.mark.s3
def test_promote_partial_failure_staging_not_cleaned(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """When one file in an assembly fails, NO staged files for that assembly are deleted.

    Preserving staging on partial failure lets an operator re-run without
    re-staging and without losing the partially-promoted state.
    """
    file_ok = PurePosixPath(f"{_ACC1}_genomic.fna.gz")
    file_fail = PurePosixPath(f"{_ACC1}_protein.faa.gz")
    _stage(mock_s3_client, _STG1, {file_ok: b"ok", file_fail: b"fail"})
    staged_ok = _STG1 / file_ok
    staged_fail = _STG1 / file_fail

    # Make download_file raise for exactly the failing key
    dynamic_client = cast("DownloadFileClient", mock_s3_client)
    original_download = dynamic_client.download_file

    def _download_one_fail(Bucket: str, Key: str, Filename: str, **kw: object) -> None:  # noqa: N803
        if Key == f"{staged_fail}":
            msg = "simulated download failure"
            raise RuntimeError(msg)
        return original_download(Bucket=Bucket, Key=Key, Filename=Filename, **kw)

    dynamic_client.download_file = _download_one_fail

    report = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    assert report["failed"] == 1
    # Staging files must still be present (cleanup skipped due to failure)
    for key in (staged_ok, staged_fail):
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(key))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK, (
            f"Expected staged file to survive partial failure: {key}"
        )


@pytest.mark.s3
def test_promote_partial_failure_failed_count(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    r"""report[\"failed\"] reflects the number of files that could not be promoted."""
    file_names = [
        PurePosixPath(f"{_ACC1}_genomic.fna.gz"),
        PurePosixPath(f"{_ACC1}_protein.faa.gz"),
        PurePosixPath(f"{_ACC1}_rna.fna.gz"),
    ]
    _stage(mock_s3_client, _STG1, dict.fromkeys(file_names, b"data"))

    failing_key = f"{_STG1 / file_names[1]}"
    dynamic_client = cast("DownloadFileClient", mock_s3_client)
    original_download = dynamic_client.download_file

    def _download_middle_fail(Bucket: str, Key: str, Filename: str, **kw: object) -> None:  # noqa: N803
        if Key == str(failing_key):
            msg = "simulated failure"
            raise RuntimeError(msg)
        return original_download(Bucket=Bucket, Key=Key, Filename=Filename, **kw)

    dynamic_client.download_file = _download_middle_fail

    report = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    assert report["failed"] == 1
    assert report["promoted"] == 2  # noqa: PLR2004


@pytest.mark.s3
def test_promote_two_assemblies_independent_cleanup(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """A fully successful assembly cleans up its staging even when another assembly partially fails.

    Assembly 1 fully succeeds → staging cleared.
    Assembly 2 has one failing file → staging NOT cleared.
    """
    # Assembly 1: two files, both succeed
    _stage(
        mock_s3_client,
        _STG1,
        {PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"g1", PurePosixPath(f"{_ACC1}_protein.faa.gz"): b"p1"},
    )
    # Assembly 2: two files, one will fail
    _stage(
        mock_s3_client,
        _STG2,
        {PurePosixPath(f"{_ACC2}_genomic.fna.gz"): b"g2", PurePosixPath(f"{_ACC2}_protein.faa.gz"): b"p2"},
    )
    failing_key = f"{_STG2 / _ACC2}_protein.faa.gz"
    dynamic_client = cast("DownloadFileClient", mock_s3_client)
    original_download = dynamic_client.download_file

    def _patched(Bucket: str, Key: str, Filename: str, **kw: object) -> None:  # noqa: N803
        if Key == f"{failing_key}":
            msg = "simulated failure"
            raise RuntimeError(msg)
        return original_download(Bucket=Bucket, Key=Key, Filename=Filename, **kw)

    dynamic_client.download_file = _patched

    report = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    assert report["failed"] == 1

    # Assembly 1 staging must be gone
    for fname in (f"{_ACC1}_genomic.fna.gz", f"{_ACC1}_protein.faa.gz"):
        result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(_STG1 / fname))
        assert result.get("KeyCount", 0) == 0, f"Assembly 1 staging should be cleaned: {fname}"

    # Assembly 2 staging must remain (partial failure)
    for fname in (f"{_ACC2}_genomic.fna.gz", f"{_ACC2}_protein.faa.gz"):
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(_STG2 / fname))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK, (
            f"Assembly 2 staging must survive partial failure: {fname}"
        )


@pytest.mark.s3
def test_promote_multi_assembly_all_succeed_all_cleaned(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Two assemblies both fully succeed → all staged files removed for both."""
    _stage(mock_s3_client, _STG1, {PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"g1"})
    _stage(mock_s3_client, _STG2, {PurePosixPath(f"{_ACC2}_genomic.fna.gz"): b"g2"})

    report = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    assert report["promoted"] == 2  # noqa: PLR2004
    assert report["failed"] == 0

    for stg, fname, lkh in (
        (_STG1, f"{_ACC1}_genomic.fna.gz", _LKH1),
        (_STG2, f"{_ACC2}_genomic.fna.gz", _LKH2),
    ):
        result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(stg / fname))
        assert result.get("KeyCount", 0) == 0, f"Staging not cleaned: {fname}"
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(lkh / fname))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_promote_dry_run_multi_file_no_writes_no_cleanup(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """dry_run with multiple files writes nothing to final path and does not delete staging."""
    file_names = [f"{_ACC1}_genomic.fna.gz", f"{_ACC1}_protein.faa.gz", f"{_ACC1}_rna.fna.gz"]
    staged_keys = _stage(mock_s3_client, _STG1, {PurePosixPath(f): f.encode() for f in file_names})

    report = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
        dry_run=True,
    )

    assert report["promoted"] == len(file_names)
    assert report["dry_run"] is True

    # Final path must be empty
    result = mock_s3_client.list_objects_v2(Bucket=str(TEST_BUCKET), Prefix=str(_LKH1))
    assert result.get("KeyCount", 0) == 0

    # Staging keys must survive
    for key in staged_keys:
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(key))
        assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK, f"Staging deleted during dry-run: {key}"


@pytest.mark.s3
def test_promote_skips_non_raw_data_paths(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Files outside raw_data/ (e.g. download_report.json) are silently skipped."""
    # Stage a real data file alongside non-promotable files
    _stage(mock_s3_client, _STG1, {PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"data"})
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(_STAGE_PREFIX / "download_report.json"), Body=b"{}")
    mock_s3_client.put_object(Bucket=str(TEST_BUCKET), Key=str(_STAGE_PREFIX / "logs/run.log"), Body=b"logs")

    report = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    assert report["promoted"] == 1  # only the .fna.gz
    assert report["failed"] == 0


@pytest.mark.s3
def test_promote_idempotent_second_run_on_empty_staging(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Second promote run after staging has been cleaned promotes 0 files without error."""
    _stage(mock_s3_client, _STG1, {PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"data"})

    report1 = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )
    report2 = promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    assert report1["promoted"] == 1
    assert report2["promoted"] == 0
    assert report2["failed"] == 0

    # Final key still present after second run
    resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=f"{_LKH1 / _ACC1}_genomic.fna.gz")
    assert resp["ResponseMetadata"]["HTTPStatusCode"] == HTTPStatus.OK


@pytest.mark.s3
def test_promote_multi_file_md5_per_file(
    mock_s3_client: botocore.client.BaseClient,
) -> None:
    """Each promoted file carries the MD5 matching its own content, not another file's."""
    files = {
        PurePosixPath(f"{_ACC1}_genomic.fna.gz"): b"GENOMIC_UNIQUE",
        PurePosixPath(f"{_ACC1}_protein.faa.gz"): b"PROTEIN_UNIQUE",
        PurePosixPath(f"{_ACC1}_rna.fna.gz"): b"RNA_UNIQUE",
    }
    _stage(mock_s3_client, _STG1, files, with_md5=True)

    promote_from_s3(
        staging_key_prefix=_STAGE_PREFIX,
        staging_bucket=TEST_BUCKET,
        lakehouse_bucket=TEST_BUCKET,
        lakehouse_key_prefix=DEFAULT_LAKEHOUSE_KEY_PREFIX,
    )

    for fname, content in files.items():
        expected_md5 = hashlib.md5(content).hexdigest()  # noqa: S324
        resp = mock_s3_client.head_object(Bucket=str(TEST_BUCKET), Key=str(_LKH1 / fname))
        assert resp["Metadata"].get("md5") == expected_md5, f"Wrong MD5 on {fname}"
