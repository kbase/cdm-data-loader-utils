"""Tests of the PDB holdings file parser."""

import json
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.parsers.pdb_holdings_file import _parse_pdb_record, parse_holdings_file
from cdm_data_loaders.pdb.constants import HoldingsFile, HoldingsFileSchemas, PDBRecord

_HOLDINGS_ID_ONLY = HoldingsFile(filename=Path("dummy.json.gz"), schema=HoldingsFileSchemas.ID_ONLY)
_HOLDINGS_ID_DATE = HoldingsFile(filename=Path("dummy.json.gz"), schema=HoldingsFileSchemas.ID_DATE)
_HOLDINGS_BAD = HoldingsFile(filename=Path("dummy.json.gz"), schema="bad_schema")  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("holdings", "data", "expected"),
    [
        pytest.param(
            _HOLDINGS_ID_ONLY,
            {"PDB_00001ABC": {}, "PDB_00001DEF": {}},
            {
                "pdb_00001abc": PDBRecord(id="pdb_00001abc"),
                "pdb_00001def": PDBRecord(id="pdb_00001def"),
            },
            id="id-only-extracts-ids",
        ),
        pytest.param(
            _HOLDINGS_ID_DATE,
            {"PDB_00001ABC": "2024-01-15", "PDB_00001DEF": "2024-02-20"},
            {
                "pdb_00001abc": PDBRecord(id="pdb_00001abc", last_modified="2024-01-15"),
                "pdb_00001def": PDBRecord(id="pdb_00001def", last_modified="2024-02-20"),
            },
            id="id-date-extracts-ids-and-dates",
        ),
        pytest.param(_HOLDINGS_ID_ONLY, {}, {}, id="id-only-empty-file"),
        pytest.param(_HOLDINGS_ID_DATE, {}, {}, id="id-date-empty-file"),
    ],
)
def test_parse_holdings_file_pass_parses_records(
    holdings: HoldingsFile, data: dict[str, Any], expected: dict[str, PDBRecord]
) -> None:
    """Test that holdings file bytes parse into PDBRecord dicts keyed on lowercased PDB IDs."""
    assert parse_holdings_file(holdings, json.dumps(data).encode("utf-8")) == expected


def test_parse_holdings_file_fail_non_dict_json_raises_type_error() -> None:
    """Test that JSON data that is not a dict raises TypeError."""
    data = json.dumps(["not", "a", "dict"]).encode("utf-8")
    with pytest.raises(TypeError, match="expected to contain a JSON dict"):
        parse_holdings_file(_HOLDINGS_ID_ONLY, data)


@pytest.mark.parametrize(
    ("holdings", "records", "expected"),
    [
        pytest.param(
            _HOLDINGS_ID_ONLY,
            {"PDB_00001ABC": {"status": "RELEASED"}},
            {"pdb_00001abc": PDBRecord(id="pdb_00001abc")},
            id="id-only-extracts-keys-ignores-values",
        ),
        pytest.param(
            _HOLDINGS_ID_ONLY,
            {"PDB_AAAAAAAA": {}},
            {"pdb_aaaaaaaa": PDBRecord(id="pdb_aaaaaaaa")},
            id="id-only-lowercases-ids",
        ),
        pytest.param(_HOLDINGS_ID_ONLY, {}, {}, id="id-only-empty"),
        pytest.param(
            _HOLDINGS_ID_DATE,
            {"PDB_00001ABC": "2024-01-15", "PDB_00001DEF": "2024-02-20"},
            {
                "pdb_00001abc": PDBRecord(id="pdb_00001abc", last_modified="2024-01-15"),
                "pdb_00001def": PDBRecord(id="pdb_00001def", last_modified="2024-02-20"),
            },
            id="id-date-extracts-keys-and-dates",
        ),
        pytest.param(
            _HOLDINGS_ID_DATE,
            {"PDB_AAAAAAAA": "2024-03-01"},
            {"pdb_aaaaaaaa": PDBRecord(id="pdb_aaaaaaaa", last_modified="2024-03-01")},
            id="id-date-lowercases-ids",
        ),
        pytest.param(_HOLDINGS_ID_DATE, {}, {}, id="id-date-empty"),
    ],
)
def test_parse_pdb_record_pass_parses_records(
    holdings: HoldingsFile, records: dict[str, Any], expected: dict[str, PDBRecord]
) -> None:
    """Test that raw records parse into PDBRecords keyed on lowercased PDB IDs."""
    assert _parse_pdb_record(holdings, records) == expected


def test_parse_pdb_record_fail_invalid_schema_raises_value_error() -> None:
    """Test that an invalid holdings file schema raises ValueError."""
    with pytest.raises(ValueError, match="Invalid holdings file schema"):
        _parse_pdb_record(_HOLDINGS_BAD, {"PDB_00001ABC": {}})
