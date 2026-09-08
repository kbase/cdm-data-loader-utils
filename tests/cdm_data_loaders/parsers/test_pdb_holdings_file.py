"""
Tests of the PDB holdings file parser.
"""

import json
from pathlib import Path

import pytest

from cdm_data_loaders.parsers.pdb_holdings_file import _parse_pdb_record, parse_holdings_file
from cdm_data_loaders.pdb.constants import HoldingsFile, HoldingsFileSchemas, PDBRecord

_HOLDINGS_ID_ONLY = HoldingsFile(filename=Path("dummy.json.gz"), schema=HoldingsFileSchemas.ID_ONLY)
_HOLDINGS_ID_DATE = HoldingsFile(filename=Path("dummy.json.gz"), schema=HoldingsFileSchemas.ID_DATE)
_HOLDINGS_BAD = HoldingsFile(filename=Path("dummy.json.gz"), schema="bad_schema")  # type: ignore[arg-type]

# ---------------------------------------------------------------------------
# parse_holdings_file
# ---------------------------------------------------------------------------


class TestParseHoldingsFile:
    """Tests for parse_holdings_file()."""

    @pytest.mark.parametrize(
        ("file", "data", "expected_records"),
        [
            pytest.param(
                HoldingsFile(
                    filename=Path("current_file_holdings.json.gz"),
                    schema=HoldingsFileSchemas.ID_ONLY,
                ),
                json.dumps({"PDB_00001ABC": {}, "PDB_00001DEF": {}}).encode("utf-8"),
                {"pdb_00001abc": "", "pdb_00001def": ""},
            ),
            pytest.param(
                HoldingsFile(
                    filename=Path("released_structures_last_modified_dates.json.gz"),
                    schema=HoldingsFileSchemas.ID_DATE,
                ),
                json.dumps({"PDB_00001ABC": "2024-01-15"}).encode("utf-8"),
                {"pdb_00001abc": "2024-01-15"},
            ),
            pytest.param(
                HoldingsFile(
                    filename=Path("all_removed_entries.json.gz"),
                    schema=HoldingsFileSchemas.ID_ONLY,
                ),
                json.dumps({"PDB_00009999": {}}).encode("utf-8"),
                {"pdb_00009999": ""},
            ),
        ],
    )
    def test_downloads_and_parses_all_three_files(
        self, file: HoldingsFile, data: bytes, expected_records: dict[str, str | None]
    ) -> None:
        """Tests that each holdings file type is properly loaded."""
        result = parse_holdings_file(file, data)
        for key, date in expected_records.items():
            assert key in result
            assert result[key].id == key
            assert result[key].last_modified == date

    def test_non_dict_response_raises_type_error(self) -> None:
        """Ensures invalid JSON data raises TypeError."""
        file = HoldingsFile(
            filename=Path("invalid.json.gz"),
            schema=HoldingsFileSchemas.ID_ONLY,
        )
        data = json.dumps(["not", "a", "dict"]).encode("utf-8")
        with pytest.raises(TypeError, match="expected to contain a JSON dict"):
            parse_holdings_file(file, data)


# ---------------------------------------------------------------------------
# _parse_pdb_record
# ---------------------------------------------------------------------------


class TestParsePdbRecord:
    """Parameterized cases for testing _parse_pdb_record()."""

    @pytest.mark.parametrize(
        ("holdings", "raw", "expected"),
        [
            pytest.param(
                _HOLDINGS_ID_ONLY,
                {"PDB_00001ABC": {"status": "RELEASED"}, "PDB_00001DEF": {}},
                {"pdb_00001abc": PDBRecord(id="pdb_00001abc"), "pdb_00001def": PDBRecord(id="pdb_00001def")},
                id="id-only-extracts-keys",
            ),
            pytest.param(
                _HOLDINGS_ID_ONLY,
                {"PDB_AAAAAAAA": {}},
                {"pdb_aaaaaaaa": PDBRecord(id="pdb_aaaaaaaa")},
                id="id-only-normalizes-case",
            ),
            pytest.param(_HOLDINGS_ID_ONLY, {}, {}, id="id-only-empty"),
            pytest.param(
                _HOLDINGS_ID_DATE,
                {"PDB_00001ABC": "2024-01-15", "PDB_00001DEF": "2024-02-20"},
                {
                    "pdb_00001abc": PDBRecord(id="pdb_00001abc", last_modified="2024-01-15"),
                    "pdb_00001def": PDBRecord(id="pdb_00001def", last_modified="2024-02-20"),
                },
                id="id-date-extracts-key-and-date",
            ),
            pytest.param(
                _HOLDINGS_ID_DATE,
                {"PDB_AAAAAAAA": "2024-03-01"},
                {"pdb_aaaaaaaa": PDBRecord(id="pdb_aaaaaaaa", last_modified="2024-03-01")},
                id="id-date-normalizes-case",
            ),
            pytest.param(_HOLDINGS_ID_DATE, {}, {}, id="id-date-empty"),
        ],
    )
    def test_parse_pdb_record(self, holdings: HoldingsFile, raw: dict, expected: dict) -> None:
        """Tests parsed output is as expected."""
        assert _parse_pdb_record(holdings, raw) == expected

    def test_invalid_schema_raises_value_error(self) -> None:
        """Ensures ValueError is raised when an invalid schema is specified."""
        with pytest.raises(ValueError, match="Invalid holdings file schema"):
            _parse_pdb_record(_HOLDINGS_BAD, {"PDB_00001ABC": {}})
