"""End-to-end tests for the PDB manifest-generating pipeline.

Tests in ``TestRunManifestGenerationE2E`` and
``TestSnapshotManifestRoundTripCeph`` run the full ``run_manifest_generation``
pipeline against a real CEPH test store (PDB HTTP downloads mocked) and are
marked ``requires_ceph`` + ``slow_test``.
"""

from datetime import date
from pathlib import Path, PurePosixPath
from types import SimpleNamespace
from typing import Final, cast

import pytest
from types_boto3_s3 import S3Client

import cdm_data_loaders.pipelines.pdb_manifest as pdb_manifest_mod
from cdm_data_loaders.pdb.constants import (
    DEFAULT_DESTINATION_PREFIX,
    HoldingsFileTypes,
    PDBRecord,
)
from cdm_data_loaders.pipelines.pdb_manifest import (
    PdbManfestSettings,
    _generate_snapshot_from_s3_state,
    _save_holdings_snapshot,
    run_manifest_generation,
)
from tests.conftest import TEST_BUCKET

# A prefix under which we seed fake PDB objects in CEPH
_PDB_KEY_PREFIX = DEFAULT_DESTINATION_PREFIX / "raw_data"


def _seed_fake_pdb_objects(
    s3: S3Client,
    bucket: PurePosixPath,
    pdb_ids: list[str],
    prefix: PurePosixPath = _PDB_KEY_PREFIX,
) -> None:
    """Upload a minimal placeholder object for each PDB ID into CEPH."""
    for pdb_id in pdb_ids:
        key = str(prefix / pdb_id / f"{pdb_id}_model.cif.gz")
        s3.put_object(Bucket=str(bucket), Key=key, Body=b"placeholder")


# Full pipeline end-to-end tests

# Holdings data shared across E2E tests.
# pdb_aaaaaaaa is intentionally absent from LAST_MODIFIED; should show up in missing_dates.
_E2E_CURRENT: dict[str, PDBRecord] = {
    "pdb_00001abc": PDBRecord(id="pdb_00001abc"),
    "pdb_00001def": PDBRecord(id="pdb_00001def"),
    "pdb_aaaaaaaa": PDBRecord(id="pdb_aaaaaaaa"),
}
_E2E_LAST_MODIFIED: dict[str, PDBRecord] = {
    "pdb_00001abc": PDBRecord(id="pdb_00001abc", last_modified="2024-01-15"),
    "pdb_00001def": PDBRecord(id="pdb_00001def", last_modified="2024-02-20"),
}
_E2E_REMOVED: dict[str, PDBRecord] = {}

_SNAPSHOT_PREFIX = PurePosixPath("e2e-pdb-test")
_SNAPSHOT_FILENAME = PurePosixPath("current_snapshot.json.gz")


@pytest.fixture
def pdb_test_bucket() -> PurePosixPath:
    """Test bucket as a PPP."""
    return PurePosixPath(TEST_BUCKET)


@pytest.fixture(autouse=True)
def _mock_holdings_download(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace live PDB downloads with controlled in-memory data."""
    mock_raw = {
        HoldingsFileTypes.CURRENT: dict(_E2E_CURRENT),
        HoldingsFileTypes.LAST_MODIFIED: dict(_E2E_LAST_MODIFIED),
        HoldingsFileTypes.REMOVED: dict(_E2E_REMOVED),
    }
    monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_files", lambda: mock_raw)


def _make_config(
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
    **kwargs: object,
) -> PdbManfestSettings:
    defaults: dict = {
        "bootstrap_date": None,
        "skip_diff": False,
        "regex_filter": None,
        "destination_bucket": pdb_test_bucket,
        "destination_prefix": _PDB_KEY_PREFIX,
        "holdings_snapshot_path": _SNAPSHOT_FILENAME,
        "output_path": tmp_path,
    }
    defaults.update(kwargs)
    return cast("PdbManfestSettings", SimpleNamespace(**defaults))


def _read_manifest(tmp_path: Path, filename: str) -> list[str]:
    return (tmp_path / filename).read_text().splitlines()


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_skip_diff_all_current_records_are_new(
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
) -> None:
    """With skip_diff=True, every current ID appears as new in the transfer manifest."""
    config = _make_config(pdb_test_bucket, tmp_path, skip_diff=True)
    run_manifest_generation(config)

    transfer = sorted(_read_manifest(tmp_path, "transfer_manifest.txt"))
    assert transfer == sorted(_E2E_CURRENT.keys())
    assert _read_manifest(tmp_path, "updated_manifest.txt") == []
    assert _read_manifest(tmp_path, "removed_manifest.txt") == []
    assert _read_manifest(tmp_path, "missing_dates.txt") == ["pdb_aaaaaaaa"]


@pytest.mark.s3
def test_bootstrap_date_uses_s3_state_as_previous(
    mock_s3_client: S3Client,
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
) -> None:
    """With bootstrap_date, the S3 store determines the previous snapshot."""
    # pdb_00001def: in S3, bootstrap date older than current: updated.
    # pdb_00009999: in S3 but not in current holdings: removed.
    _seed_fake_pdb_objects(mock_s3_client, pdb_test_bucket, ["pdb_00001def", "pdb_00009999"])
    bootstrap_date = date(2020, 1, 1)

    config = _make_config(pdb_test_bucket, tmp_path, bootstrap_date=bootstrap_date)
    run_manifest_generation(config)

    # pdb_00001abc and pdb_aaaaaaaa are not in the store: new
    transfer = _read_manifest(tmp_path, "transfer_manifest.txt")
    assert "pdb_00001abc" in transfer
    assert "pdb_00001def" in transfer  # also in transfer (updated IDs included)
    assert "pdb_aaaaaaaa" in transfer

    assert _read_manifest(tmp_path, "updated_manifest.txt") == ["pdb_00001def"]
    assert _read_manifest(tmp_path, "removed_manifest.txt") == ["pdb_00009999"]
    assert _read_manifest(tmp_path, "missing_dates.txt") == ["pdb_aaaaaaaa"]


@pytest.mark.s3
def test_with_existing_snapshot_incremental_diff(
    mock_s3_client: S3Client,
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
) -> None:
    """An S3 snapshot drives incremental diffing; only changed records appear."""
    # pdb_00001def: in snapshot with older date: updated.
    # pdb_00009999: in snapshot but not in current holdings: removed.
    previous = {
        "pdb_00001def": PDBRecord(id="pdb_00001def", last_modified="2023-06-01"),
        "pdb_00009999": PDBRecord(id="pdb_00009999", last_modified="2022-01-01"),
    }
    local_snapshot = tmp_path / "input_snapshot.json.gz"
    _save_holdings_snapshot(previous, local_snapshot)
    snapshot_key = _SNAPSHOT_PREFIX / _SNAPSHOT_FILENAME
    mock_s3_client.upload_file(
        Filename=str(local_snapshot),
        Bucket=str(pdb_test_bucket),
        Key=str(snapshot_key),
    )

    config = _make_config(
        pdb_test_bucket,
        tmp_path,
        destination_prefix=_SNAPSHOT_PREFIX,
    )
    run_manifest_generation(config)

    # pdb_00001abc and pdb_aaaaaaaa were not in snapshot: new
    transfer = _read_manifest(tmp_path, "transfer_manifest.txt")
    assert "pdb_00001abc" in transfer
    assert "pdb_00001def" in transfer  # also in transfer (updated IDs included)
    assert "pdb_aaaaaaaa" in transfer

    assert _read_manifest(tmp_path, "updated_manifest.txt") == ["pdb_00001def"]
    assert _read_manifest(tmp_path, "removed_manifest.txt") == ["pdb_00009999"]
    assert _read_manifest(tmp_path, "missing_dates.txt") == ["pdb_aaaaaaaa"]


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize(
    ("regex_filter", "expected_transfer", "expected_missing"),
    [
        pytest.param(
            "pdb_00001",
            ["pdb_00001abc", "pdb_00001def"],
            [],
            id="prefix-filter",
        ),
        pytest.param(
            "pdb_aaa",
            ["pdb_aaaaaaaa"],
            ["pdb_aaaaaaaa"],
            id="single-match-still-missing-date",
        ),
    ],
)
def test_regex_filter(
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
    regex_filter: str,
    expected_transfer: list[str],
    expected_missing: list[str],
) -> None:
    """Regex filter restricts which IDs appear across all output manifest files."""
    config = _make_config(pdb_test_bucket, tmp_path, skip_diff=True, regex_filter=regex_filter)
    run_manifest_generation(config)

    transfer = sorted(_read_manifest(tmp_path, "transfer_manifest.txt"))
    assert transfer == sorted(expected_transfer)
    assert _read_manifest(tmp_path, "updated_manifest.txt") == []
    assert _read_manifest(tmp_path, "removed_manifest.txt") == []
    assert sorted(_read_manifest(tmp_path, "missing_dates.txt")) == sorted(expected_missing)


# ---------------------------------------------------------------------------
# Round-trip: generate snapshot from S3 state then run manifest generation
# ---------------------------------------------------------------------------

# Snapshot filename used by the round-trip tests (distinct from the other E2E tests).
_ROUNDTRIP_SNAPSHOT_FILENAME = PurePosixPath("roundtrip_snapshot.json.gz")

# IDs seeded into CEPH for the round-trip tests.
_ROUNDTRIP_IDS: list[str] = ["pdb_00001abc", "pdb_00001def", "pdb_aaaaaaaa"]

# A date well in the past so we can easily manufacture "newer" dates in holdings.
_BOOTSTRAP_DATE: Final[date] = date(2020, 1, 1)
_NEWER_DATE: Final[str] = "2024-06-15"


def _make_roundtrip_config(
    pdb_test_bucket: PurePosixPath,
    output_path: Path,
    **kwargs: object,
) -> PdbManfestSettings:
    defaults: dict = {
        "bootstrap_date": None,
        "skip_diff": False,
        "regex_filter": None,
        "destination_bucket": pdb_test_bucket,
        "destination_prefix": _PDB_KEY_PREFIX,
        "holdings_snapshot_path": _ROUNDTRIP_SNAPSHOT_FILENAME,
        "output_path": output_path,
    }
    defaults.update(kwargs)
    return cast("PdbManfestSettings", SimpleNamespace(**defaults))


def _build_and_upload_snapshot(
    s3: S3Client,
    bucket: PurePosixPath,
    snap_date: date,
    local_path: Path,
) -> dict[str, PDBRecord]:
    """Scan the CEPH store, save a snapshot locally, and upload it.

    Returns the generated snapshot dict so callers can inspect its contents.
    The snapshot is uploaded to the key read by ``_download_holdings_snapshot``
    when using the config produced by ``_make_roundtrip_config``.
    """
    snapshot = _generate_snapshot_from_s3_state(
        bucket=bucket,
        key_prefix=_PDB_KEY_PREFIX,
        date=snap_date,
    )
    _save_holdings_snapshot(snapshot, local_path)
    s3.upload_file(
        Filename=str(local_path),
        Bucket=str(bucket),
        Key=str(_PDB_KEY_PREFIX / _ROUNDTRIP_SNAPSHOT_FILENAME),
    )
    return snapshot


def _mock_holdings(
    monkeypatch: pytest.MonkeyPatch,
    current: dict[str, PDBRecord],
    last_modified: dict[str, PDBRecord] | None = None,
    removed: dict[str, PDBRecord] | None = None,
) -> None:
    monkeypatch.setattr(
        pdb_manifest_mod,
        "_download_holdings_files",
        lambda: {
            HoldingsFileTypes.CURRENT: dict(current),
            HoldingsFileTypes.LAST_MODIFIED: dict(last_modified or {}),
            HoldingsFileTypes.REMOVED: dict(removed or {}),
        },
    )


def test_no_changes_when_store_matches_holdings(
    mock_s3_client: S3Client,
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When the snapshot exactly matches the current holdings, all manifests are empty."""
    _seed_fake_pdb_objects(mock_s3_client, pdb_test_bucket, _ROUNDTRIP_IDS)
    _build_and_upload_snapshot(mock_s3_client, pdb_test_bucket, _BOOTSTRAP_DATE, tmp_path / "snap.json.gz")

    # Holdings current = same IDs; last-modified dates match the bootstrap date.
    current = {id_: PDBRecord(id=id_) for id_ in _ROUNDTRIP_IDS}
    last_modified = {id_: PDBRecord(id=id_, last_modified=_BOOTSTRAP_DATE.isoformat()) for id_ in _ROUNDTRIP_IDS}
    _mock_holdings(monkeypatch, current, last_modified)

    run_manifest_generation(_make_roundtrip_config(pdb_test_bucket, tmp_path / "out"))

    assert _read_manifest(tmp_path / "out", "transfer_manifest.txt") == []
    assert _read_manifest(tmp_path / "out", "updated_manifest.txt") == []
    assert _read_manifest(tmp_path / "out", "removed_manifest.txt") == []
    assert _read_manifest(tmp_path / "out", "missing_dates.txt") == []


def test_records_removed_from_store_reappear_in_transfer_manifest(
    mock_s3_client: S3Client,
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """IDs dropped from the S3 store between runs re-emerge as new in transfer_manifest.

    Workflow:
    1. Seed all IDs; generate snapshot. expect no changes on first run.
    2. Delete some IDs from the store; regenerate snapshot, expect deleted IDs
       are absent from the new snapshot but still present in the mocked
       holdings: they appear as *new* in ``transfer_manifest.txt``.
    """
    removed_ids = ["pdb_00001def", "pdb_aaaaaaaa"]
    current = {id_: PDBRecord(id=id_) for id_ in _ROUNDTRIP_IDS}
    last_modified = {id_: PDBRecord(id=id_, last_modified=_BOOTSTRAP_DATE.isoformat()) for id_ in _ROUNDTRIP_IDS}
    _mock_holdings(monkeypatch, current, last_modified)

    # First run: baseline:nothing should appear in any manifest
    _seed_fake_pdb_objects(mock_s3_client, pdb_test_bucket, _ROUNDTRIP_IDS)
    _build_and_upload_snapshot(mock_s3_client, pdb_test_bucket, _BOOTSTRAP_DATE, tmp_path / "snap1.json.gz")
    run_manifest_generation(_make_roundtrip_config(pdb_test_bucket, tmp_path / "run1"))
    assert _read_manifest(tmp_path / "run1", "transfer_manifest.txt") == []

    # Remove some records from the store
    for pdb_id in removed_ids:
        key = str(_PDB_KEY_PREFIX / pdb_id / f"{pdb_id}_model.cif.gz")
        mock_s3_client.delete_object(Bucket=str(pdb_test_bucket), Key=key)

    # Second run: removed IDs are absent from new snapshot: should appear in transfer_manifest.txt
    _build_and_upload_snapshot(mock_s3_client, pdb_test_bucket, _BOOTSTRAP_DATE, tmp_path / "snap2.json.gz")
    run_manifest_generation(_make_roundtrip_config(pdb_test_bucket, tmp_path / "run2"))

    transfer = sorted(_read_manifest(tmp_path / "run2", "transfer_manifest.txt"))
    assert transfer == sorted(removed_ids)
    assert _read_manifest(tmp_path / "run2", "updated_manifest.txt") == []
    assert _read_manifest(tmp_path / "run2", "removed_manifest.txt") == []


def test_holdings_newer_date_triggers_updated_manifest(
    mock_s3_client: S3Client,
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Records whose holdings last-modified date is newer than the snapshot appear in updated_manifest."""
    ids = ["pdb_00001abc", "pdb_00001def"]
    _seed_fake_pdb_objects(mock_s3_client, pdb_test_bucket, ids)
    _build_and_upload_snapshot(mock_s3_client, pdb_test_bucket, _BOOTSTRAP_DATE, tmp_path / "snap.json.gz")

    # pdb_00001abc has a newer date in holdings: should appear in updated
    # pdb_00001def has the same date as the snapshot: should not appear in updated
    current = {id_: PDBRecord(id=id_) for id_ in ids}
    last_modified = {
        "pdb_00001abc": PDBRecord(id="pdb_00001abc", last_modified=_NEWER_DATE),
        "pdb_00001def": PDBRecord(id="pdb_00001def", last_modified=_BOOTSTRAP_DATE.isoformat()),
    }
    _mock_holdings(monkeypatch, current, last_modified)

    run_manifest_generation(_make_roundtrip_config(pdb_test_bucket, tmp_path / "out"))

    assert _read_manifest(tmp_path / "out", "updated_manifest.txt") == ["pdb_00001abc"]
    transfer = _read_manifest(tmp_path / "out", "transfer_manifest.txt")
    assert "pdb_00001abc" in transfer
    assert "pdb_00001def" not in transfer


def test_id_absent_from_holdings_appears_in_removed_manifest(
    mock_s3_client: S3Client,
    pdb_test_bucket: PurePosixPath,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An ID present in the snapshot but absent from current holdings appears in removed_manifest."""
    kept_in_holdings = ["pdb_00001abc", "pdb_00001def"]
    dropped_from_holdings = ["pdb_aaaaaaaa"]

    # Seed all IDs so the snapshot contains all three.
    _seed_fake_pdb_objects(mock_s3_client, pdb_test_bucket, _ROUNDTRIP_IDS)
    _build_and_upload_snapshot(mock_s3_client, pdb_test_bucket, _BOOTSTRAP_DATE, tmp_path / "snap.json.gz")

    # Holdings no longer lists pdb_aaaaaaaa as current.
    current = {id_: PDBRecord(id=id_) for id_ in kept_in_holdings}
    last_modified = {id_: PDBRecord(id=id_, last_modified=_BOOTSTRAP_DATE.isoformat()) for id_ in kept_in_holdings}
    _mock_holdings(monkeypatch, current, last_modified)

    run_manifest_generation(_make_roundtrip_config(pdb_test_bucket, tmp_path / "out"))

    assert _read_manifest(tmp_path / "out", "removed_manifest.txt") == dropped_from_holdings
    assert _read_manifest(tmp_path / "out", "transfer_manifest.txt") == []
    assert _read_manifest(tmp_path / "out", "updated_manifest.txt") == []
