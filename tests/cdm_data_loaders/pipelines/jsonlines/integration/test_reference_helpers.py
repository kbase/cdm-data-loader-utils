"""Exact checks for the normalization used by JSONL reference comparisons."""

from collections.abc import Callable
from copy import deepcopy
from datetime import UTC, date, datetime
from typing import Any

import duckdb
import pytest

from tests.cdm_data_loaders.pipelines.jsonlines.integration.jsonlines_reference_helpers import (
    canonicalize_pipeline_entry,
    canonicalize_reference_entry,
    reconstructed_entries,
    reference_entries_from_duckdb,
)


@pytest.mark.parametrize(
    "canonicalize", [canonicalize_reference_entry, canonicalize_pipeline_entry], ids=["reference", "pipeline"]
)
@pytest.mark.parametrize(
    ("entry", "expected"),
    [
        pytest.param(
            {"null": None, "owner": {}, "nested": {"absent": None}, "zero": 0, "false": False, "empty": ""},
            {"zero": 0, "false": False, "empty": ""},
            id="nulls-empty-objects-and-falsy-scalars",
        ),
        pytest.param(
            {"values": ["second", "first", "second"], "empty_list": []},
            {"values": [{"value": "second"}, {"value": "first"}, {"value": "second"}], "empty_list": []},
            id="scalar-list-order-and-duplicates",
        ),
        pytest.param(
            {"values": [{"value": "first"}, {"accession": "a", "unused": None}]},
            {"values": [{"value": "first"}, {"accession": "a"}]},
            id="already-reconstructed-list",
        ),
        pytest.param(
            {"gc_percent": {"v_double": 42}, "unrelated": {"v_double": 7}},
            {"gc_percent": 42.0, "unrelated": {"v_double": 7}},
            id="only-known-numeric-variant",
        ),
        pytest.param({"gc_percent": True}, {"gc_percent": True}, id="bool-is-not-double"),
        pytest.param(
            {"submission_date": "2014-08-21T09:55:24.333", "release_date": "2024-01-01"},
            {"submission_date": datetime(2014, 8, 21, 9, 55, 24, 333000, tzinfo=UTC), "release_date": date(2024, 1, 1)},
            id="naive-timestamp-and-date",
        ),
        pytest.param(
            {"when": "2014-08-21T09:55:24.333+00:00"},
            {"when": datetime(2014, 8, 21, 9, 55, 24, 333000, tzinfo=UTC)},
            id="aware-timestamp",
        ),
        pytest.param(
            {"when": "not-a-timestamp", "release_date": "not-a-date", "name": "2024-01-01"},
            {"when": "not-a-timestamp", "release_date": "not-a-date", "name": "2024-01-01"},
            id="invalid-and-unrelated-dates-unchanged",
        ),
    ],
)
def test_canonicalize_entry_pass_exact_storage_normalization(
    canonicalize: Callable[[dict[str, Any]], dict[str, Any]], entry: dict[str, Any], expected: dict[str, Any]
) -> None:
    """Canonicalization handles known storage differences without mutation or repeated wrapping."""
    original = deepcopy(entry)
    result = canonicalize(entry)
    assert result == expected
    assert canonicalize(result) == expected
    assert entry == original
    if "gc_percent" in expected:
        assert type(result["gc_percent"]) is type(expected["gc_percent"])


def test_reconstructed_entries_pass_rebuilds_ordered_children_and_ignores_metadata() -> None:
    """Root selection excludes rejection and metadata tables while preserving child order."""
    tables = {
        "dataset": [{"_dlt_id": "root", "accession": "assembly", "organism__tax_id": 42}],
        "dataset__tags": [
            {"_dlt_id": "second", "_dlt_parent_id": "root", "_dlt_list_idx": 1, "value": "b"},
            {"_dlt_id": "first", "_dlt_parent_id": "root", "_dlt_list_idx": 0, "value": "a"},
        ],
        "dataset_rejected": [{"raw_record": "bad"}],
        "_dlt_loads": [{"load_id": "ignored"}],
    }
    assert reconstructed_entries(tables, "dataset") == [
        {"accession": "assembly", "organism": {"tax_id": 42}, "tags": [{"value": "a"}, {"value": "b"}]}
    ]


@pytest.mark.parametrize("tables", [{}, {"_dlt_loads": []}, {"dataset": []}], ids=["empty", "metadata-only", "no-rows"])
def test_reconstructed_entries_pass_empty_tables(tables: dict[str, list[dict[str, Any]]]) -> None:
    """Empty output has no reconstructed records."""
    assert reconstructed_entries(tables) == []


@pytest.mark.parametrize(
    "tables", [{"first": [], "second": []}, {"dataset__child": []}], ids=["ambiguous-roots", "missing-root"]
)
def test_reconstructed_entries_fail_invalid_root(tables: dict[str, list[dict[str, Any]]]) -> None:
    """Ambiguous or orphaned tables fail instead of silently producing a partial reference match."""
    with pytest.raises(ValueError, match="Expected exactly one top-level table"):
        reconstructed_entries(tables)


def test_reference_entries_from_duckdb_pass_restores_json_values() -> None:
    """DuckDB dates become source strings and padded null fields are removed recursively."""
    with duckdb.connect() as connection:
        connection.execute(
            "CREATE TABLE reference AS SELECT 'assembly' AS accession, "
            "{'release_date': DATE '2024-01-01', 'missing': NULL} AS assembly_info, "
            "[{'name': 'submitter', 'missing': NULL}] AS additional_submitters"
        )
        assert reference_entries_from_duckdb(connection) == [
            {
                "accession": "assembly",
                "assembly_info": {"release_date": "2024-01-01"},
                "additional_submitters": [{"name": "submitter"}],
            }
        ]
