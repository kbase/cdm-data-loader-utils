"""Unit tests for the PDB manifest-generating pipeline."""

import gzip
import json
import tempfile
from collections.abc import Callable
from datetime import date
from pathlib import Path, PurePosixPath
from types import SimpleNamespace
from typing import IO, ClassVar, cast
from unittest.mock import ANY, MagicMock, patch

import pytest
from botocore.exceptions import ClientError

import cdm_data_loaders.pipelines.pdb_manifest as pdb_manifest_mod
from cdm_data_loaders.pdb.constants import (
    HoldingsFileTypes,
    ManifestData,
    PDBRecord,
)
from cdm_data_loaders.pipelines.pdb_manifest import (
    PdbManfestSettings,
    _download_holdings_files,
    _download_holdings_snapshot,
    _extract_id_from_s3_key,
    _generate_manifest_data,
    _generate_snapshot_from_s3_state,
    _load_holdings_snapshot,
    _save_holdings_snapshot,
    _save_manifest_files,
    _save_summary_file,
    run_manifest_generation,
)
from cdm_data_loaders.utils.file_transfer.s3 import client

# ---------------------------------------------------------------------------
# Shared test data
# ---------------------------------------------------------------------------

_ID_A = "pdb_00001abc"
_ID_B = "pdb_00001def"
_ID_C = "pdb_aaaaaaaa"
_ID_D = "pdb_03942abd"

_DATE_OLD = "2023-01-01"
_DATE_NEW = "2024-06-15"

_MANIFEST_ONE_OF_EACH = ManifestData(new=[_ID_A], updated=[_ID_B], removed=[_ID_C], missing_dates=[_ID_D])
_MANIFEST_THREE_NEW = ManifestData(new=[_ID_A, _ID_B, _ID_C], updated=[], removed=[], missing_dates=[])


def _make_gzipped_response(data: dict | list[str]) -> MagicMock:
    """Return a mock context-manager response yielding gzipped JSON bytes."""
    compressed = gzip.compress(json.dumps(data).encode())
    resp = MagicMock()
    resp.read.return_value = compressed
    resp.__enter__ = lambda s: s
    resp.__exit__ = MagicMock(return_value=False)
    return resp


# ---------------------------------------------------------------------------
# _generate_manifest_data
# ---------------------------------------------------------------------------


class TestGenerateManifestData:
    """Test for _generate_manifest_data()."""

    @pytest.mark.parametrize(
        ("current", "previous", "expected_new"),
        [
            pytest.param(
                {
                    _ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW),
                    _ID_B: PDBRecord(id=_ID_B, last_modified=_DATE_NEW),
                },
                {},
                [_ID_A, _ID_B],
                id="all-new-no-previous",
            ),
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW)},
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_OLD)},
                [],
                id="already-in-previous",
            ),
        ],
    )
    def test_new_field(
        self,
        current: dict[str, PDBRecord],
        previous: dict[str, PDBRecord],
        expected_new: list[str],
    ) -> None:
        """Tests new records are properly determined."""
        result = _generate_manifest_data(current, removed=set(), previous=previous)
        assert sorted(result.new) == sorted(expected_new)

    @pytest.mark.parametrize(
        ("current", "previous", "expect_in_updated"),
        [
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW)},
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_OLD)},
                True,
                id="newer-date",
            ),
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW)},
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW)},
                False,
                id="same-date",
            ),
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_OLD)},
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW)},
                False,
                id="older-date",
            ),
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW)},
                {_ID_A: PDBRecord(id=_ID_A, last_modified="")},
                True,
                id="prev-no-date",
            ),
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A, last_modified="")},
                {_ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_OLD)},
                False,
                id="current-no-date",
            ),
        ],
    )
    def test_updated_field(
        self,
        current: dict[str, PDBRecord],
        previous: dict[str, PDBRecord],
        expect_in_updated: bool,
    ) -> None:
        """Ensures various combinations of specified and unspecified dates are handled."""
        result = _generate_manifest_data(current, removed=set(), previous=previous)
        assert (_ID_A in result.updated) is expect_in_updated

    @pytest.mark.parametrize(
        ("current", "removed_set", "previous", "expected_removed"),
        [
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A)},
                set(),
                {_ID_A: PDBRecord(id=_ID_A), _ID_B: PDBRecord(id=_ID_B)},
                [_ID_B],
                id="vanished-from-current",
            ),
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A)},
                {_ID_B},
                {_ID_A: PDBRecord(id=_ID_A), _ID_B: PDBRecord(id=_ID_B)},
                [_ID_B],
                id="explicit-removed-in-previous",
            ),
            pytest.param(
                {_ID_A: PDBRecord(id=_ID_A)},
                {_ID_C},
                {_ID_A: PDBRecord(id=_ID_A)},
                [],
                id="explicit-removed-not-in-previous",
            ),
            pytest.param(
                {},
                set(),
                {_ID_B: PDBRecord(id=_ID_B), _ID_A: PDBRecord(id=_ID_A)},
                sorted([_ID_A, _ID_B]),
                id="all-removed-sorted",
            ),
        ],
    )
    def test_removed_field(
        self,
        current: dict[str, PDBRecord],
        removed_set: set[str],
        previous: dict[str, PDBRecord],
        expected_removed: list[str],
    ) -> None:
        """Ensures removed records are correctly determined based on holdings files and previous state."""
        result = _generate_manifest_data(current, removed=removed_set, previous=previous)
        assert result.removed == expected_removed

    def test_missing_dates_propagated(self) -> None:
        """Ensures missing dates are properly set."""
        result = _generate_manifest_data({}, removed=set(), missing_dates=[_ID_A, _ID_B])
        assert result.missing_dates == [_ID_A, _ID_B]

    def test_new_and_updated_are_disjoint(self) -> None:
        """An entry cannot be both new and updated."""
        current = {
            _ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW),  # new (not in previous)
            _ID_B: PDBRecord(id=_ID_B, last_modified=_DATE_NEW),  # updated (in previous, newer date)
        }
        previous = {_ID_B: PDBRecord(id=_ID_B, last_modified=_DATE_OLD)}
        result = _generate_manifest_data(current, removed=set(), previous=previous)
        assert set(result.new) & set(result.updated) == set()


# ---------------------------------------------------------------------------
# _extract_id_from_s3_key
# ---------------------------------------------------------------------------


class TestExtractIdFromS3Key:
    """Tests for _extract_id_from_s3_key()."""

    @pytest.mark.parametrize(
        ("key", "expected"),
        [
            pytest.param(
                "tenant-general-warehouse/pdb/raw_data/pdb_00001abc/file.cif.gz",
                "pdb_00001abc",
                id="nested-path",
            ),
            pytest.param("some/unrelated/path.txt", None, id="no-match"),
            pytest.param("data/PDB_00001ABC/file.cif", "pdb_00001abc", id="uppercase-normalized"),
            pytest.param(
                "pdb_00001abc/subdir/pdb_00001def/file",
                "pdb_00001abc",
                id="first-match-returned",
            ),
        ],
    )
    def test_extract_id(self, key: str, expected: str | None) -> None:
        """Ensures ids are properly extracted when present."""
        assert _extract_id_from_s3_key(key) == expected


# ---------------------------------------------------------------------------
# _save_holdings_snapshot / _load_holdings_snapshot
# ---------------------------------------------------------------------------


class TestHoldingsSnapshotRoundTrip:
    """Tests for IO of snapshop data."""

    @pytest.mark.parametrize(
        "records",
        [
            pytest.param(
                {
                    _ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW),
                    _ID_B: PDBRecord(id=_ID_B, last_modified=_DATE_OLD),
                },
                id="two-records",
            ),
            pytest.param({}, id="empty"),
        ],
    )
    def test_round_trip(self, tmp_path: Path, records: dict[str, PDBRecord]) -> None:
        """Tests save/load of snapshot data."""
        path = tmp_path / "snapshot.json.gz"
        _save_holdings_snapshot(records, path)
        assert _load_holdings_snapshot(path) == records

    def test_missing_last_modified_defaults_to_empty_string(self, tmp_path: Path) -> None:
        """Snapshot entries with no last_modified field are loaded with empty string."""
        payload = {_ID_A: PDBRecord(id=_ID_A)}  # no "last_modified" key
        path = tmp_path / "partial.json.gz"
        _save_holdings_snapshot(payload, path)
        assert _load_holdings_snapshot(path)[_ID_A].last_modified == ""


# ---------------------------------------------------------------------------
# _generate_snapshot_from_s3_state
# ---------------------------------------------------------------------------


class TestGenerateSnapshotFromS3StateUnit:
    """Unit tests for _generate_snapshot_from_s3_state() with a mocked S3 client."""

    def _build_mock_s3(self, pages: list[list[str]]) -> MagicMock:
        """Return a mock S3 client whose paginator yields one page per list of key strings."""
        mock_s3 = MagicMock()
        paginator = MagicMock()
        mock_s3.get_paginator.return_value = paginator
        page_dicts = [{"Contents": [{"Key": k} for k in keys]} if keys else {} for keys in pages]
        paginator.paginate.return_value = page_dicts
        return mock_s3

    def test_empty_store_returns_empty_dict(self) -> None:
        """Empty S3 store produces an empty snapshot."""
        mock_s3 = self._build_mock_s3([[]])
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _generate_snapshot_from_s3_state(
                PurePosixPath("bucket"), PurePosixPath("prefix"), date(2024, 1, 1)
            )
        assert result == {}

    def test_returns_record_for_each_valid_key(self) -> None:
        """Objects with extractable PDB IDs produce one record per ID."""
        keys = [
            "prefix/pdb_00001abc/pdb_00001abc_model.cif.gz",
            "prefix/pdb_00001def/pdb_00001def_data.cif.gz",
        ]
        mock_s3 = self._build_mock_s3([keys])
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _generate_snapshot_from_s3_state(
                PurePosixPath("bucket"), PurePosixPath("prefix"), date(2024, 1, 1)
            )
        assert set(result.keys()) == {"pdb_00001abc", "pdb_00001def"}

    def test_sets_provided_date_on_all_records(self) -> None:
        """The bootstrap date is written to every record's last_modified field."""
        keys = [
            "prefix/pdb_00001abc/model.cif.gz",
            "prefix/pdb_00001def/model.cif.gz",
        ]
        mock_s3 = self._build_mock_s3([keys])
        bootstrap = date(2023, 6, 15)
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _generate_snapshot_from_s3_state(PurePosixPath("bucket"), PurePosixPath("prefix"), bootstrap)
        for rec in result.values():
            assert rec.last_modified == bootstrap.isoformat()

    def test_deduplicates_multiple_objects_per_pdb_id(self) -> None:
        """Multiple S3 objects sharing a PDB ID produce exactly one snapshot entry."""
        pdb_id = "pdb_00001abc"
        keys = [
            f"prefix/{pdb_id}/{pdb_id}_model.cif.gz",
            f"prefix/{pdb_id}/{pdb_id}_data.cif.gz",
            f"prefix/{pdb_id}/{pdb_id}_info.json",
        ]
        mock_s3 = self._build_mock_s3([keys])
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _generate_snapshot_from_s3_state(
                PurePosixPath("bucket"), PurePosixPath("prefix"), date(2024, 1, 1)
            )
        assert len(result) == 1
        assert pdb_id in result

    def test_skips_keys_without_valid_pdb_id(self) -> None:
        """Keys that contain no recognisable PDB ID are silently ignored."""
        keys = [
            "some/unrelated/path.txt",
            "prefix/logs/error.log",
            "prefix/pdb_00001abc/model.cif.gz",  # the only valid one
        ]
        mock_s3 = self._build_mock_s3([keys])
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _generate_snapshot_from_s3_state(
                PurePosixPath("bucket"), PurePosixPath("prefix"), date(2024, 1, 1)
            )
        assert set(result.keys()) == {"pdb_00001abc"}

    def test_aggregates_records_across_multiple_pages(self) -> None:
        """Records distributed across paginator pages are all collected."""
        mock_s3 = self._build_mock_s3(
            [
                ["prefix/pdb_00001abc/model.cif.gz"],
                ["prefix/pdb_00001def/model.cif.gz"],
                ["prefix/pdb_aaaaaaaa/model.cif.gz"],
            ]
        )
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _generate_snapshot_from_s3_state(
                PurePosixPath("bucket"), PurePosixPath("prefix"), date(2024, 1, 1)
            )
        assert set(result.keys()) == {"pdb_00001abc", "pdb_00001def", "pdb_aaaaaaaa"}

    def test_paginator_called_with_correct_bucket_and_prefix(self) -> None:
        """list_objects_v2 paginator receives the exact bucket name and prefix."""
        mock_s3 = self._build_mock_s3([[]])
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            _generate_snapshot_from_s3_state(PurePosixPath("my-bucket"), PurePosixPath("a/b/c"), date(2024, 1, 1))
        mock_s3.get_paginator.assert_called_once_with("list_objects_v2")
        mock_s3.get_paginator.return_value.paginate.assert_called_once_with(Bucket="my-bucket", Prefix="a/b/c")

    def test_page_with_no_contents_key_is_handled_gracefully(self) -> None:
        """A paginator page with no 'Contents' key doesn't raise and contributes no records."""
        mock_s3 = self._build_mock_s3(
            [
                [],  # renders as {} (no Contents key)
                ["prefix/pdb_00001abc/model.cif.gz"],
            ]
        )
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _generate_snapshot_from_s3_state(
                PurePosixPath("bucket"), PurePosixPath("prefix"), date(2024, 1, 1)
            )
        assert set(result.keys()) == {"pdb_00001abc"}


# ---------------------------------------------------------------------------
# _download_holdings_snapshot
# ---------------------------------------------------------------------------


class TestDownloadHoldingsSnapshotUnit:
    """Unit tests for _download_holdings_snapshot() with a mocked S3 client."""

    _RECORDS: ClassVar[dict[str, PDBRecord]] = {
        _ID_A: PDBRecord(id=_ID_A, last_modified=_DATE_NEW),
        _ID_B: PDBRecord(id=_ID_B, last_modified=_DATE_OLD),
    }

    def _writing_side_effect(self, records: dict[str, PDBRecord]) -> Callable[..., None]:
        """Return a download_file side_effect that writes a valid snapshot to Filename."""

        def _effect(**kwargs: str) -> None:
            _save_holdings_snapshot(records, Path(kwargs["Filename"]))

        return _effect

    def test_returns_parsed_snapshot_on_success(self) -> None:
        """A successful download returns the correctly deserialised snapshot dict."""
        mock_s3 = MagicMock()
        mock_s3.download_file.side_effect = self._writing_side_effect(self._RECORDS)
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            result = _download_holdings_snapshot(PurePosixPath("bucket"), PurePosixPath("path/snapshot.json.gz"))
        assert result == self._RECORDS

    def test_calls_download_file_with_correct_bucket_and_key(self) -> None:
        """download_file is invoked with the exact bucket and key passed to the function."""
        mock_s3 = MagicMock()
        mock_s3.download_file.side_effect = self._writing_side_effect({})
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            _download_holdings_snapshot(PurePosixPath("my-bucket"), PurePosixPath("some/key.json.gz"))
        mock_s3.download_file.assert_called_once_with(Bucket="my-bucket", Key="some/key.json.gz", Filename=ANY)

    def test_client_error_raises_value_error(self) -> None:
        """A ClientError from S3 is converted into a descriptive ValueError."""
        mock_s3 = MagicMock()
        mock_s3.download_file.side_effect = ClientError(
            {"Error": {"Code": "NoSuchKey", "Message": "Not found"}}, "GetObject"
        )
        with (
            patch.object(client, "get_s3_client", return_value=mock_s3),
            pytest.raises(ValueError, match="not found"),
        ):
            _download_holdings_snapshot(PurePosixPath("bucket"), PurePosixPath("missing/key.json.gz"))

    def test_error_message_includes_bucket_and_key(self) -> None:
        """The ValueError message contains both the bucket name and the key path."""
        mock_s3 = MagicMock()
        mock_s3.download_file.side_effect = ClientError({"Error": {"Code": "NoSuchKey", "Message": ""}}, "GetObject")
        with (
            patch.object(client, "get_s3_client", return_value=mock_s3),
            pytest.raises(ValueError, match="Snapshot file") as exc_info,
        ):
            _download_holdings_snapshot(PurePosixPath("my-bucket"), PurePosixPath("path/to/snap.json.gz"))
        msg = str(exc_info.value)
        assert "my-bucket" in msg
        assert "path/to/snap.json.gz" in msg

    def test_temp_file_is_cleaned_up_on_success(self) -> None:
        """The temporary download file is deleted after a successful call."""
        captured: list[str] = []
        mock_s3 = MagicMock()

        def write_and_capture(**kwargs: str) -> None:
            _save_holdings_snapshot({}, Path(kwargs["Filename"]))
            captured.append(kwargs["Filename"])

        mock_s3.download_file.side_effect = write_and_capture
        with patch.object(client, "get_s3_client", return_value=mock_s3):
            _download_holdings_snapshot(PurePosixPath("b"), PurePosixPath("k.json.gz"))

        assert captured, "download_file was never called"
        assert not Path(captured[0]).exists(), f"Temp file {captured[0]} was not cleaned up"

    def test_temp_file_is_cleaned_up_on_error(self) -> None:
        """The temporary download file is deleted even when a ClientError is raised."""
        captured: list[Path] = []
        real_ntf = tempfile.NamedTemporaryFile

        def recording_ntf(**kwargs: object) -> IO[bytes]:
            handle = real_ntf(**kwargs)  # type: ignore[arg-type]
            captured.append(Path(handle.name))
            return handle

        mock_s3 = MagicMock()
        mock_s3.download_file.side_effect = ClientError({"Error": {"Code": "NoSuchKey", "Message": ""}}, "GetObject")
        with (
            patch.object(client, "get_s3_client", return_value=mock_s3),
            patch("cdm_data_loaders.pipelines.pdb_manifest.tempfile.NamedTemporaryFile", side_effect=recording_ntf),
            pytest.raises(ValueError, match="not found"),
        ):
            _download_holdings_snapshot(PurePosixPath("b"), PurePosixPath("k.json.gz"))

        assert captured, "NamedTemporaryFile was never called"
        assert not captured[0].exists(), f"Temp file {captured[0]} was not cleaned up after error"


# ---------------------------------------------------------------------------
# _download_holdings_files
# ---------------------------------------------------------------------------


class TestDownloadHoldingsFiles:
    """Tests for _download_holdings_files()."""

    def test_downloads_and_parses_all_three_files(self) -> None:
        """Tests that each holdings file type is properly loaded."""
        raw_current = {"PDB_00001ABC": {}, "PDB_00001DEF": {}}
        raw_dates = {"PDB_00001ABC": "2024-01-15"}
        raw_removed = {"PDB_00009999": {}}

        responses = [
            _make_gzipped_response(raw_current),
            _make_gzipped_response(raw_dates),
            _make_gzipped_response(raw_removed),
        ]

        with patch("cdm_data_loaders.pipelines.pdb_manifest.urlopen", side_effect=responses):
            result = _download_holdings_files()

        assert HoldingsFileTypes.CURRENT in result
        assert HoldingsFileTypes.LAST_MODIFIED in result
        assert HoldingsFileTypes.REMOVED in result
        assert "pdb_00001abc" in result[HoldingsFileTypes.CURRENT]
        assert "pdb_00001def" in result[HoldingsFileTypes.CURRENT]
        assert result[HoldingsFileTypes.LAST_MODIFIED]["pdb_00001abc"].last_modified == "2024-01-15"
        assert "pdb_00009999" in result[HoldingsFileTypes.REMOVED]


# ---------------------------------------------------------------------------
# _save_manifest_files
# ---------------------------------------------------------------------------


class TestSaveManifestFiles:
    """Tests for _save_manifest_files()."""

    @pytest.mark.parametrize(
        ("data", "filename", "expected_lines"),
        [
            pytest.param(
                _MANIFEST_ONE_OF_EACH,
                "transfer_manifest.txt",
                [_ID_A, _ID_B],
                id="transfer-contains-new-and-updated",
            ),
            pytest.param(
                _MANIFEST_ONE_OF_EACH,
                "updated_manifest.txt",
                [_ID_B],
                id="updated-contains-only-updated",
            ),
            pytest.param(
                _MANIFEST_ONE_OF_EACH,
                "removed_manifest.txt",
                [_ID_C],
                id="removed-contains-only-removed",
            ),
            pytest.param(
                _MANIFEST_ONE_OF_EACH,
                "missing_dates.txt",
                [_ID_D],
                id="missing-dates-empty",
            ),
            pytest.param(
                _MANIFEST_THREE_NEW,
                "transfer_manifest.txt",
                [_ID_A, _ID_B, _ID_C],
                id="each-id-on-its-own-line",
            ),
        ],
    )
    def test_file_contents(self, tmp_path: Path, data: ManifestData, filename: str, expected_lines: list[str]) -> None:
        """Ensures manifest file contents are correct."""
        _save_manifest_files(data, tmp_path)
        lines = (tmp_path / filename).read_text().splitlines()
        assert sorted(lines) == sorted(expected_lines)


# ---------------------------------------------------------------------------
# _save_summary_file
# ---------------------------------------------------------------------------


class TestSaveSummaryFile:
    """Tests for _save_summary_file()."""

    @pytest.mark.parametrize(
        ("data", "regex_filter", "expected_values", "absent_keys"),
        [
            pytest.param(
                ManifestData(new=[_ID_A, _ID_B], updated=[_ID_C], removed=[], missing_dates=[_ID_A]),
                None,
                {"new": 2, "updated": 1, "removed": 0, "missing_dates": 1},
                ["regex_filter"],
                id="counts-correct-no-filter",
            ),
            pytest.param(
                ManifestData(new=[], updated=[], removed=[], missing_dates=[]),
                "pdb_0000",
                {"regex_filter": "pdb_0000"},
                [],
                id="regex-filter-included",
            ),
            pytest.param(
                ManifestData(new=[], updated=[], removed=[], missing_dates=[]),
                None,
                {},
                ["regex_filter"],
                id="regex-filter-omitted",
            ),
        ],
    )
    def test_summary_file(
        self,
        tmp_path: Path,
        data: ManifestData,
        regex_filter: str | None,
        expected_values: dict,
        absent_keys: list[str],
    ) -> None:
        """Ensure metrics are output in summary file."""
        _save_summary_file(data, regex_filter, tmp_path)
        summary = json.loads((tmp_path / "summary.json").read_text())
        for key, val in expected_values.items():
            assert summary[key] == val
        for key in absent_keys:
            assert key not in summary
        assert "generated_at" in summary
        assert isinstance(summary["saved_to"], str)


# ---------------------------------------------------------------------------
# run_manifest_generation
# ---------------------------------------------------------------------------


class TestRunManifestGeneration:
    """Unit tests for the run_manifest_generation orchestration.

    All HTTP and S3 calls are mocked so no external services are required.
    """

    # Controlled holdings data: pdb_aaaaaaaa is intentionally absent from
    # LAST_MODIFIED to exercise the missing-dates path.
    _CURRENT: ClassVar[dict[str, PDBRecord]] = {
        "pdb_00001abc": PDBRecord(id="pdb_00001abc"),
        "pdb_00001def": PDBRecord(id="pdb_00001def"),
        "pdb_aaaaaaaa": PDBRecord(id="pdb_aaaaaaaa"),
    }
    _LAST_MODIFIED: ClassVar[dict[str, PDBRecord]] = {
        "pdb_00001abc": PDBRecord(id="pdb_00001abc", last_modified="2024-01-15"),
        "pdb_00001def": PDBRecord(id="pdb_00001def", last_modified="2024-02-20"),
    }
    _REMOVED: ClassVar[dict[str, PDBRecord]] = {}

    @pytest.fixture(autouse=True)
    def _patch_holdings_download(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """Replace HTTP downloads with controlled in-memory data."""
        mock_raw = {
            HoldingsFileTypes.CURRENT: dict(self._CURRENT),
            HoldingsFileTypes.LAST_MODIFIED: dict(self._LAST_MODIFIED),
            HoldingsFileTypes.REMOVED: dict(self._REMOVED),
        }
        monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_files", lambda: mock_raw)

    def _make_config(self, tmp_path: Path, **kwargs: object) -> PdbManfestSettings:
        defaults: dict = {
            "bootstrap_date": None,
            "skip_diff": False,
            "regex_filter": None,
            "destination_bucket": PurePosixPath("test-bucket"),
            "destination_prefix": PurePosixPath("test/prefix"),
            "holdings_snapshot_path": PurePosixPath("snapshot.json.gz"),
            "output_path": tmp_path,
        }
        defaults.update(kwargs)
        return cast("PdbManfestSettings", SimpleNamespace(**defaults))

    def _read_manifest(self, tmp_path: Path, filename: str) -> list[str]:
        return (tmp_path / filename).read_text().splitlines()

    def test_skip_diff_all_records_are_new(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
        """With skip_diff=True no snapshot is loaded and every current ID is new."""
        mock_dl = MagicMock()
        monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_snapshot", mock_dl)

        run_manifest_generation(self._make_config(tmp_path, skip_diff=True))

        mock_dl.assert_not_called()
        transfer = self._read_manifest(tmp_path, "transfer_manifest.txt")
        assert sorted(transfer) == sorted(self._CURRENT.keys())
        assert self._read_manifest(tmp_path, "updated_manifest.txt") == []
        assert self._read_manifest(tmp_path, "removed_manifest.txt") == []

    def test_bootstrap_date_calls_generate_snapshot(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
        """With bootstrap_date set, _generate_snapshot_from_s3_state is invoked."""
        mock_gen = MagicMock(return_value={})
        monkeypatch.setattr(pdb_manifest_mod, "_generate_snapshot_from_s3_state", mock_gen)

        bootstrap_date = date(2020, 1, 1)
        run_manifest_generation(self._make_config(tmp_path, bootstrap_date=bootstrap_date))

        mock_gen.assert_called_once_with(
            PurePosixPath("test-bucket"),
            PurePosixPath("test") / "prefix",
            bootstrap_date,
        )

    def test_normal_mode_downloads_snapshot(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
        """Without skip_diff or bootstrap_date, _download_holdings_snapshot is invoked."""
        mock_dl = MagicMock(return_value={})
        monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_snapshot", mock_dl)

        run_manifest_generation(self._make_config(tmp_path))

        mock_dl.assert_called_once_with(
            bucket=PurePosixPath("test-bucket"),
            key=PurePosixPath("test") / "prefix" / "snapshot.json.gz",
        )

    @pytest.mark.parametrize(
        ("regex_filter", "expected_ids"),
        [
            pytest.param("pdb_00001", ["pdb_00001abc", "pdb_00001def"], id="prefix-match"),
            pytest.param("pdb_aaa", ["pdb_aaaaaaaa"], id="single-match"),
            pytest.param("pdb_xxxxx", [], id="no-match"),
        ],
    )
    def test_regex_filter_limits_transfer_manifest(
        self,
        monkeypatch: pytest.MonkeyPatch,
        tmp_path: Path,
        regex_filter: str,
        expected_ids: list[str],
    ) -> None:
        """Regex filter restricts which IDs appear in the output manifests."""
        monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_snapshot", MagicMock(return_value={}))

        run_manifest_generation(self._make_config(tmp_path, skip_diff=True, regex_filter=regex_filter))

        transfer = self._read_manifest(tmp_path, "transfer_manifest.txt")
        assert sorted(transfer) == sorted(expected_ids)

    def test_missing_dates_computed_as_current_minus_last_modified(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        """IDs present in CURRENT but absent from LAST_MODIFIED appear in missing_dates.txt."""
        monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_snapshot", MagicMock(return_value={}))

        run_manifest_generation(self._make_config(tmp_path, skip_diff=True))

        missing = self._read_manifest(tmp_path, "missing_dates.txt")
        assert missing == ["pdb_aaaaaaaa"]

    def test_summary_file_written(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
        """A summary.json is written with correct counts."""
        monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_snapshot", MagicMock(return_value={}))

        run_manifest_generation(self._make_config(tmp_path, skip_diff=True))

        summary = json.loads((tmp_path / "summary.json").read_text())
        total = summary["new"] + summary["updated"]
        assert total == len(self._CURRENT)
        assert summary["missing_dates"] == 1

    def test_output_directory_created_when_missing(self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
        """run_manifest_generation creates output_path and its parents when they don't exist."""
        monkeypatch.setattr(pdb_manifest_mod, "_download_holdings_snapshot", MagicMock(return_value={}))

        output_path = tmp_path / "nested" / "output"
        assert not output_path.exists()

        run_manifest_generation(self._make_config(tmp_path, skip_diff=True, output_path=output_path))

        assert output_path.is_dir()
        assert (output_path / "transfer_manifest.txt").exists()
