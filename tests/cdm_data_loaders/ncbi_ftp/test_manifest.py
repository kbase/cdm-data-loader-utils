"""Tests for ncbi_ftp.manifest module — assembly summary parsing, diff, filtering, writing."""

import json
from datetime import UTC, datetime
from pathlib import Path, PurePosixPath
from types import SimpleNamespace
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
from botocore.exceptions import ClientError

from cdm_data_loaders.ncbi_ftp.manifest import (
    AssemblyRecord,
    DiffResult,
    _extract_accession_dir_and_id_from_s3_key,
    _ftp_dir_from_url,
    accession_prefix,
    compute_diff,
    filter_by_prefix_range,
    get_latest_assembly_paths,
    parse_assembly_summary,
    scan_store_to_synthetic_summary,
    verify_transfer_candidates,
    write_diff_summary,
    write_removed_manifest,
    write_transfer_manifest,
    write_updated_manifest,
)

from .conftest import ACC_PATH_215, ACC_PATH_405, SAMPLE_SUMMARY

_EXPECTED_TWO = 2


# parse_assembly_summary


_EXPECTED_ASSEMBLIES = {
    "GCF_000001215.4": AssemblyRecord(
        accession="GCF_000001215.4",
        status="latest",
        seq_rel_date="2014/10/21",
        ftp_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT",
        assembly_dir=PurePosixPath("GCF_000001215.4_Release_6_plus_ISO1_MT"),
    ),
    "GCF_000001405.40": AssemblyRecord(
        accession="GCF_000001405.40",
        status="latest",
        seq_rel_date="2022/02/03",
        ftp_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14",
        assembly_dir=PurePosixPath("GCF_000001405.40_GRCh38.p14"),
    ),
    "GCF_000005845.2": AssemblyRecord(
        accession="GCF_000005845.2",
        status="replaced",
        seq_rel_date="2013/09/26",
        ftp_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2",
        assembly_dir=PurePosixPath("GCF_000005845.2_ASM584v2"),
    ),
    "GCF_000009999.1": AssemblyRecord(
        accession="GCF_000009999.1",
        status="suppressed",
        seq_rel_date="2010/01/01",
        ftp_url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/999/GCF_000009999.1_ASM999v1",
        assembly_dir=PurePosixPath("GCF_000009999.1_ASM999v1"),
    ),
    # GCF_000099999.1 is excluded because ftp_path == "na"
}


def test_parse_assembly_summary() -> None:
    """SAMPLE_SUMMARY is parsed to the expected assemblies."""
    assert parse_assembly_summary(SAMPLE_SUMMARY) == _EXPECTED_ASSEMBLIES


def test_parse_assembly_summary_empty() -> None:
    """Comment-only input returns empty dict."""
    assert parse_assembly_summary("# comment only\n") == {}


@pytest.mark.parametrize("source", ["file", "file_str", "list_of_lines"])
def test_parse_assembly_summary_input_types(source: str, tmp_path: Path) -> None:
    """Parsing works from a file path, string path, and list of lines."""
    if source == "list_of_lines":
        arg: str | Path | list[str] = list(SAMPLE_SUMMARY.splitlines(keepends=True))
    else:
        f = tmp_path / "summary.tsv"
        f.write_text(SAMPLE_SUMMARY)
        arg = f if source == "file" else str(f)
    assert parse_assembly_summary(arg) == _EXPECTED_ASSEMBLIES


# get_latest_assembly_paths


def test_get_latest_assembly_paths() -> None:
    """Only 'latest' assemblies appear; paths are FTP directories with trailing slash."""
    assert dict(get_latest_assembly_paths(parse_assembly_summary(SAMPLE_SUMMARY))) == {
        "GCF_000001215.4": PurePosixPath("/") / "genomes" / "all" / ACC_PATH_215,
        "GCF_000001405.40": PurePosixPath("/") / "genomes" / "all" / ACC_PATH_405,
    }


def test_get_latest_assembly_paths_empty() -> None:
    """Empty input returns empty list."""
    assert get_latest_assembly_paths(parse_assembly_summary("# empty\n")) == []


# compute_diff


def test_compute_diff_new() -> None:
    """All latest assemblies are new with no previous state; result is sorted."""
    current = parse_assembly_summary(SAMPLE_SUMMARY)
    diff = compute_diff(current, previous_accessions=set())
    assert "GCF_000001215.4" in diff.new
    assert "GCF_000001405.40" in diff.new
    assert "GCF_000005845.2" not in diff.new  # replaced
    assert diff.new == sorted(diff.new)


def test_compute_diff_updated() -> None:
    """seq_rel_date moving forward marks updated; moving backward does not."""
    current = parse_assembly_summary(SAMPLE_SUMMARY)
    previous = parse_assembly_summary(SAMPLE_SUMMARY)
    previous["GCF_000001215.4"].seq_rel_date = "2010/01/01"
    assert "GCF_000001215.4" in compute_diff(current, previous_assemblies=previous).updated

    previous2 = parse_assembly_summary(SAMPLE_SUMMARY)
    previous2["GCF_000001215.4"].seq_rel_date = "2099/12/31"
    assert "GCF_000001215.4" not in compute_diff(current, previous_assemblies=previous2).updated


def test_compute_diff_removed() -> None:
    """Replaced, suppressed, and entirely-absent accessions are classified correctly."""
    current = parse_assembly_summary(SAMPLE_SUMMARY)
    assert "GCF_000005845.2" in compute_diff(current, previous_accessions={"GCF_000005845.2"}).replaced
    assert "GCF_000009999.1" in compute_diff(current, previous_accessions={"GCF_000009999.1"}).suppressed
    # Accession absent from current entirely → suppressed
    assert (
        "GCF_000001215.4"
        in compute_diff(parse_assembly_summary("# empty\n"), previous_accessions={"GCF_000001215.4"}).suppressed
    )


def test_compute_diff_scan_store_fallback() -> None:
    """Known accessions are not marked new; unknown ones are."""
    current = parse_assembly_summary(SAMPLE_SUMMARY)
    diff = compute_diff(current, previous_accessions={"GCF_000001215.4"})
    assert "GCF_000001215.4" not in diff.new
    assert "GCF_000001405.40" in diff.new


# accession_prefix & filter_by_prefix_range


@pytest.mark.parametrize(
    ("accession", "expected"),
    [
        pytest.param("GCF_000001215.4", "000", id="three_zeros"),
        pytest.param("GCF_123456789.1", "123", id="non_zero"),
        pytest.param("invalid", None, id="invalid"),
    ],
)
def test_accession_prefix(accession: str, expected: str | None) -> None:
    """3-digit prefix is extracted from the accession; invalid input raises ValueError."""
    if expected is None:
        with pytest.raises(ValueError, match="Could not parse accession"):
            accession_prefix(accession)
    else:
        assert accession_prefix(accession) == expected


def test_filter_by_prefix_range() -> None:
    """Range filter is inclusive; out-of-range excluded; no range returns all."""
    assemblies = parse_assembly_summary(SAMPLE_SUMMARY)
    assert len(filter_by_prefix_range(assemblies, "000", "000")) == len(assemblies)
    assert len(filter_by_prefix_range(assemblies, "001", "999")) == 0
    assert len(filter_by_prefix_range(assemblies)) == len(assemblies)


# Manifest writing


def test_write_transfer_manifest(tmp_path: Path) -> None:
    """Transfer manifest is written with FTP paths that start with /genomes/ and end with /."""
    current = parse_assembly_summary(SAMPLE_SUMMARY)
    diff = compute_diff(current, previous_accessions=set())
    manifest_file = tmp_path / "transfer.txt"
    paths = write_transfer_manifest(diff, current, manifest_file)
    assert len(paths) > 0
    for path in paths:
        assert path.is_relative_to("/genomes")
        parts = path.parts
        assert len(parts) == 8  #  noqa: PLR2004  /genomes/all/GCF/012/345/678/GCF_012345678_blah_blah (initial '/' counts as a part)


def test_write_removed_manifest(tmp_path: Path) -> None:
    """Removed manifest lists replaced and suppressed accessions."""
    current = parse_assembly_summary(SAMPLE_SUMMARY)
    diff = compute_diff(current, previous_accessions={"GCF_000005845.2", "GCF_000009999.1"})
    removed_file = tmp_path / "removed.txt"
    removed = write_removed_manifest(diff, removed_file)
    assert len(removed) == _EXPECTED_TWO
    lines = [line.strip() for line in removed_file.read_text().splitlines() if line.strip()]
    assert len(lines) == _EXPECTED_TWO


def test_write_updated_manifest(tmp_path: Path) -> None:
    """Updated manifest lists only updated accessions, sorted."""
    diff = DiffResult(new=["GCF_000001215.4"], updated=["GCF_000005845.2", "GCF_000001405.40"])
    updated_file = tmp_path / "updated.txt"
    updated = write_updated_manifest(diff, updated_file)
    assert len(updated) == _EXPECTED_TWO
    lines = [line.strip() for line in updated_file.read_text().splitlines() if line.strip()]
    assert lines == ["GCF_000001405.40", "GCF_000005845.2"]


def test_write_diff_summary(tmp_path: Path) -> None:
    """Diff summary JSON is written with correct counts, prefix range, and database."""
    diff = DiffResult(new=["a"], updated=["b"], replaced=["c"], suppressed=[])
    summary_file = tmp_path / "summary.json"
    summary = write_diff_summary(diff, summary_file, "refseq", "000", "003")
    assert json.loads(summary_file.read_text()) == summary
    assert {k: summary[k] for k in ("database", "counts", "prefix_range", "accessions")} == {
        "database": "refseq",
        "counts": {
            "new": 1,
            "updated": 1,
            "replaced": 1,
            "suppressed": 0,
            "total_to_transfer": 2,
            "total_to_remove": 1,
        },
        "prefix_range": {"from": "000", "to": "003"},
        "accessions": {"new": ["a"], "updated": ["b"], "replaced": ["c"], "suppressed": []},
    }


# _ftp_dir_from_url


@pytest.mark.parametrize(
    ("url", "expected", "kwargs"),
    [
        pytest.param(
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT",
            PurePosixPath("/") / "genomes" / "all" / ACC_PATH_215,
            {},
            id="https_url",
        ),
        pytest.param(
            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT",
            PurePosixPath("/") / "genomes" / "all" / ACC_PATH_215,
            {},
            id="ftp_url",
        ),
        pytest.param(
            "/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT",
            None,
            {},
            id="bare_path",
        ),
        pytest.param(
            "ftp://custom.host.example.com/genomes/all/GCF/000/001/215",
            PurePosixPath("/") / "genomes" / "all" / "GCF" / "000" / "001" / "215",
            {"ftp_host": "custom.host.example.com"},
            id="custom_ftp_host",
        ),
    ],
)
def test_ftp_dir_from_url(url: str, expected: str | None, kwargs: dict) -> None:
    """Wrapper for _ftp_dir_from_url."""
    if expected is None:
        with pytest.raises(ValueError, match="Could not parse FTP URL"):
            _ftp_dir_from_url(url, **kwargs)
    else:
        assert _ftp_dir_from_url(url, **kwargs) == expected


# verify_transfer_candidates


_BUCKET: PurePosixPath = PurePosixPath("cdm-lake")
_KEY_PREFIX: PurePosixPath = PurePosixPath("tenant-general-warehouse") / "kbase" / "datasets " / "ncbi"


def _assemblies() -> dict:
    return parse_assembly_summary(SAMPLE_SUMMARY)


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.head_object")
@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_prunes_when_all_match(
    mock_connect: MagicMock,
    mock_head: MagicMock,
) -> None:
    """Assemblies where every file matches S3 are pruned from the list."""
    mock_connect.return_value = MagicMock()

    def head_side_effect(s3_path: str) -> dict | None:
        if "_genomic.fna.gz" in s3_path:
            return {
                "ContentLength": 100,
                "Metadata": {"md5": "d41d8cd98f00b204e9800998ecf8427e"},
                "ChecksumCRC64NVME": None,
            }
        if "_protein.faa.gz" in s3_path:
            return {
                "ContentLength": 100,
                "Metadata": {"md5": "a1b2c3d4e5f6a1b2c3d4e5f6a1b2c3d4"},
                "ChecksumCRC64NVME": None,
            }
        if "_assembly_report.txt" in s3_path:
            return {
                "ContentLength": 100,
                "Metadata": {"md5": "ffffffffffffffffffffffffffffffff"},
                "ChecksumCRC64NVME": None,
            }
        return None

    mock_head.side_effect = head_side_effect
    assert verify_transfer_candidates(["GCF_000001215.4"], _assemblies(), _BUCKET, _KEY_PREFIX) == []


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.head_object")
@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_keeps_when_md5_differs(
    mock_connect: MagicMock,
    mock_head: MagicMock,
) -> None:
    """Assembly is kept when at least one file has a different MD5."""
    mock_connect.return_value = MagicMock()
    mock_head.return_value = {"ContentLength": 100, "Metadata": {"md5": "WRONG"}, "ChecksumCRC64NVME": None}
    assert verify_transfer_candidates(["GCF_000001215.4"], _assemblies(), _BUCKET, _KEY_PREFIX) == ["GCF_000001215.4"]


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.head_object")
@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_keeps_when_s3_object_missing(
    mock_connect: MagicMock,
    mock_head: MagicMock,
) -> None:
    """Assembly is kept when at least one file doesn't exist in S3."""
    mock_connect.return_value = MagicMock()

    def head_side_effect(s3_path: str) -> dict[str, Any]:
        if "GCF_000001215.4" in s3_path:
            raise ClientError({"Error": {"Code": "404", "Message": "Not found"}}, "HeadObject")
        return {}

    mock_head.side_effect = head_side_effect
    assert verify_transfer_candidates(["GCF_000001215.4"], _assemblies(), _BUCKET, _KEY_PREFIX) == ["GCF_000001215.4"]


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.head_object")
@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_keeps_when_no_md5_metadata(
    mock_connect: MagicMock,
    mock_head: MagicMock,
) -> None:
    """Assembly is kept when S3 object exists but has no md5 metadata."""
    mock_connect.return_value = MagicMock()
    mock_head.return_value = {"ContentLength": 100, "Metadata": {}, "ChecksumCRC64NVME": None}
    assert verify_transfer_candidates(["GCF_000001215.4"], _assemblies(), _BUCKET, _KEY_PREFIX) == ["GCF_000001215.4"]


@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_keeps_when_ftp_fails(
    mock_connect: MagicMock, manifest_transfer_mocks: SimpleNamespace
) -> None:
    """Assembly is kept (conservative) when md5checksums.txt cannot be fetched."""
    mock_connect.return_value = MagicMock()
    manifest_transfer_mocks.retrieve.side_effect = Exception("FTP error")
    assert verify_transfer_candidates(["GCF_000001215.4"], _assemblies(), _BUCKET, _KEY_PREFIX) == ["GCF_000001215.4"]


def test_verify_transfer_candidates_empty_input() -> None:
    """Empty accession list returns empty result without connecting."""
    with patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp", return_value=MagicMock()):
        assert verify_transfer_candidates([], {}, _BUCKET, PurePosixPath("prefix")) == []


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects", return_value=["one key"])
def test_verify_transfer_candidates_unknown_accession_kept(
    mock_connect: MagicMock,
) -> None:
    """Accessions not in assemblies dict are kept (conservative)."""
    mock_connect.return_value = MagicMock()
    assert verify_transfer_candidates(["GCF_999999999.1"], {}, _BUCKET, PurePosixPath("prefix")) == ["GCF_999999999.1"]


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.head_object")
@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_short_circuits_on_first_mismatch(
    mock_connect: MagicMock,
    mock_head: MagicMock,
) -> None:
    """Verification stops checking after the first missing/mismatched file."""
    mock_connect.return_value = MagicMock()

    def head_side_effects(s3_path: str) -> dict[str, Any]:
        if "GCF_000001215.4" in s3_path:
            raise ClientError({"Error": {"Code": "404", "Message": "Not found"}}, "HeadObject")
        return {}

    mock_head.side_effect = head_side_effects
    verify_transfer_candidates(["GCF_000001215.4"], _assemblies(), _BUCKET, _KEY_PREFIX)
    assert mock_head.call_count == 1


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.head_object")
@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_mixed(
    mock_connect: MagicMock,
    mock_head: MagicMock,
) -> None:
    """A mix of matching and non-matching assemblies: matched pruned, unmatched kept."""
    mock_connect.return_value = MagicMock()

    def head_side_effect(s3_path: str) -> dict | None:
        if "GCF_000001215.4_Release_6_plus_ISO1_MT/" in s3_path:
            if "_genomic.fna.gz" in s3_path:
                return {
                    "ContentLength": 1,
                    "Metadata": {"md5": "d41d8cd98f00b204e9800998ecf8427e"},
                    "ChecksumCRC64NVME": None,
                }
            if "_protein.faa.gz" in s3_path:
                return {
                    "ContentLength": 1,
                    "Metadata": {"md5": "a1b2c3d4e5f6a1b2c3d4e5f6a1b2c3d4"},
                    "ChecksumCRC64NVME": None,
                }
            if "_assembly_report.txt" in s3_path:
                return {
                    "ContentLength": 1,
                    "Metadata": {"md5": "ffffffffffffffffffffffffffffffff"},
                    "ChecksumCRC64NVME": None,
                }
        raise ClientError({"Error": {"Code": "404", "Message": "Not found"}}, "HeadObject")

    mock_head.side_effect = head_side_effect
    result = verify_transfer_candidates(["GCF_000001215.4", "GCF_000001405.40"], _assemblies(), _BUCKET, _KEY_PREFIX)
    assert result == ["GCF_000001405.40"]


@pytest.mark.usefixtures("manifest_transfer_mocks")
@patch("cdm_data_loaders.ncbi_ftp.manifest.connect_ftp")
def test_verify_transfer_candidates_skips_ftp_when_folder_missing(
    mock_connect: MagicMock,
    manifest_transfer_mocks: SimpleNamespace,
) -> None:
    """Accessions with no objects in S3 are confirmed without FTP round-trip."""
    manifest_transfer_mocks.list.return_value = []
    result = verify_transfer_candidates(["GCF_000001215.4"], _assemblies(), _BUCKET, _KEY_PREFIX)
    assert result == ["GCF_000001215.4"]
    mock_connect.assert_not_called()


# Synthetic summary from S3 store scan


@pytest.mark.parametrize(
    ("key", "expected"),
    [
        pytest.param(
            PurePosixPath("tenant-general-warehouse")
            / "kbase"
            / "datasets"
            / "ncbi"
            / "refseq"
            / "GCF_000001215.4_Release_6_plus_ISO1_MT"
            / "file.gz",
            ("GCF_000001215.4_Release_6_plus_ISO1_MT", "GCF_000001215.4"),
            id="long_path",
        ),
        pytest.param(
            PurePosixPath("some") / "path" / "GCA_999999999.1_whatever" / "data.txt",
            ("GCA_999999999.1_whatever", "GCA_999999999.1"),
            id="short_path",
        ),
        pytest.param(
            PurePosixPath("prefix") / "GCA_999999999.1_assembly_name" / "subdir" / "data.txt",
            ("GCA_999999999.1_assembly_name", "GCA_999999999.1"),
            id="subdir",
        ),
        pytest.param(PurePosixPath("some") / "random" / "path", None, id="no_accession"),
        pytest.param(PurePosixPath(""), None, id="empty"),
    ],
)
def test_extract_accession_dir_and_id_from_s3_key(key: PurePosixPath, expected: tuple[str, str] | None) -> None:
    """Accession is extracted from S3 key paths; invalid/empty paths return None."""
    if expected is None:
        with pytest.raises(ValueError, match="Could not parse S3 key for accession info"):
            _extract_accession_dir_and_id_from_s3_key(key)
    else:
        assert _extract_accession_dir_and_id_from_s3_key(key) == expected


def _make_mock_list_object_return_value() -> list[dict[str, Any]]:
    """Return a response for list objects with two assemblies (GCF_000001215.4, GCF_000005845.2)."""
    return [
        {
            "Key": PurePosixPath("tenant-general-warehouse")
            / "kbase"
            / "datasets"
            / "ncbi"
            / "refseq"
            / "GCF_000001215.4_Release_6/file1.gz",
            "LastModified": datetime(2024, 1, 15, tzinfo=UTC),
        },
        {
            "Key": PurePosixPath("tenant-general-warehouse")
            / "kbase"
            / "datasets"
            / "ncbi"
            / "refseq"
            / "GCF_000001215.4_Release_6/file2.gz",
            "LastModified": datetime(2024, 1, 16, tzinfo=UTC),
        },
        {
            "Key": PurePosixPath("tenant-general-warehouse")
            / "kbase"
            / "datasets"
            / "ncbi"
            / "refseq"
            / "GCF_000005845.2_Assembly"
            / "file.gz",
            "LastModified": datetime(2024, 2, 20, tzinfo=UTC),
        },
    ]


@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects")
def test_scan_store_builds_summary(mock_list_objects: MagicMock) -> None:
    """Synthetic summary is built correctly with provided release_date for all assemblies."""
    mock_list_objects.return_value = _make_mock_list_object_return_value()
    assert scan_store_to_synthetic_summary(PurePosixPath("test-bucket"), PurePosixPath("prefix"), "2024/01/31") == {
        "GCF_000001215.4": AssemblyRecord(
            accession="GCF_000001215.4",
            status="latest",
            seq_rel_date="2024/01/31",
            ftp_url="https://ftp.ncbi.nlm.nih.gov/synthetic/GCF_000001215.4_Release_6",
            assembly_dir=PurePosixPath("GCF_000001215.4_Release_6"),
        ),
        "GCF_000005845.2": AssemblyRecord(
            accession="GCF_000005845.2",
            status="latest",
            seq_rel_date="2024/01/31",
            ftp_url="https://ftp.ncbi.nlm.nih.gov/synthetic/GCF_000005845.2_Assembly",
            assembly_dir=PurePosixPath("GCF_000005845.2_Assembly"),
        ),
    }


@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects")
def test_scan_store_applies_release_date_to_all(mock_list_objects: MagicMock) -> None:
    """Provided release_date is used even when files have different LastModified dates."""
    mock_list_objects.return_value = [
        {
            "Key": "prefix/GCF_000001215.4_v1/file_newer.gz",
            "LastModified": datetime(2024, 3, 20, tzinfo=UTC),
        },
        {
            "Key": "prefix/GCF_000001215.4_v1/file_older.gz",
            "LastModified": datetime(2024, 1, 10, tzinfo=UTC),
        },
    ]
    assert (
        scan_store_to_synthetic_summary(PurePosixPath("test-bucket"), PurePosixPath("prefix"), "2024/03/31")[
            "GCF_000001215.4"
        ].seq_rel_date
        == "2024/03/31"
    )


@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects")
def test_scan_store_raises_for_invalid_release_date(mock_list_objects: MagicMock) -> None:
    """Invalid release_date format is rejected."""
    mock_list_objects.return_value = _make_mock_list_object_return_value()
    with pytest.raises(ValueError, match="Invalid release_date"):
        scan_store_to_synthetic_summary(PurePosixPath("test-bucket"), PurePosixPath("prefix"), "2024-03-31")


@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects")
def test_scan_store_invokes_progress_callback(mock_list_objects: MagicMock) -> None:
    """Progress callback is called once per unique assembly discovered."""
    mock_list_objects.return_value = _make_mock_list_object_return_value()
    calls: list[tuple[int, str]] = []
    scan_store_to_synthetic_summary(
        PurePosixPath("test-bucket"),
        PurePosixPath("prefix"),
        "2024/01/31",
        progress_callback=lambda n, a: calls.append((n, a)),
    )
    assert len(calls) == 2  # noqa: PLR2004
    assert calls[0][0] == 1
    assert calls[1][0] == 2  # noqa: PLR2004


@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects")
def test_scan_store_handles_empty_store(mock_list_objects: MagicMock) -> None:
    """Empty store returns empty dict."""
    mock_list_objects.return_value = []
    assert scan_store_to_synthetic_summary(PurePosixPath("test-bucket"), PurePosixPath("prefix"), "2024/01/31") == {}


@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects")
def test_scan_store_skips_objects_without_accession(mock_list_objects: MagicMock) -> None:
    """Objects without valid accessions in the key are skipped."""
    mock_list_objects.return_value = [
        {
            "Key": PurePosixPath("prefix") / "some" / "random" / "file.txt",
            "LastModified": datetime(2024, 1, 1, tzinfo=UTC),
        },
        {
            "Key": PurePosixPath("prefix") / "GCF_000001215.4_Assembly" / "valid_file.gz",
            "LastModified": datetime(2024, 2, 1, tzinfo=UTC),
        },
    ]
    result = scan_store_to_synthetic_summary(PurePosixPath("test-bucket"), PurePosixPath("prefix"), "2024/01/31")
    assert len(result) == 1
    assert "GCF_000001215.4" in result


@patch("cdm_data_loaders.ncbi_ftp.manifest.list_objects")
def test_scan_store_assembly_dir_survives_round_trip(mock_list_objects: MagicMock, tmp_path: Path) -> None:
    """assembly_dir is preserved after save-to-file / parse-back round-trip.

    Regression: previously ftp_path was written as "" causing assembly_dir=""
    and compute_diff flagging every assembly as updated.
    """
    mock_list_objects.return_value = [
        {
            "Key": "prefix/GCF_000001215.4_Release_6_plus_ISO1_MT/file.gz",
            "LastModified": datetime(2024, 3, 10, tzinfo=UTC),
        }
    ]

    synthetic = scan_store_to_synthetic_summary(PurePosixPath("test-bucket"), PurePosixPath("prefix"), "2024/03/10")

    # save scanned records to synthetic summary file
    out_file = tmp_path / "synthetic_summary.txt"
    with out_file.open("w") as f:
        for acc in sorted(synthetic.keys()):
            rec = synthetic[acc]
            f.write(
                f"{rec.accession}\t.\t.\t.\t.\t.\t.\t.\t.\t.\t{rec.status}\t.\t.\t.\t{rec.seq_rel_date}\t.\t.\t.\t.\t{rec.ftp_url}\t.\n"
            )

    # ensure generated summary file can be parsed to regenerate original sythentic summary records
    reparsed = parse_assembly_summary(out_file)
    assert "GCF_000001215.4" in reparsed
    reparsed_rec = reparsed["GCF_000001215.4"]
    original_rec = synthetic["GCF_000001215.4"]
    assert reparsed_rec.assembly_dir == original_rec.assembly_dir
    assert reparsed_rec.seq_rel_date == original_rec.seq_rel_date
    assert reparsed_rec.status == original_rec.status
