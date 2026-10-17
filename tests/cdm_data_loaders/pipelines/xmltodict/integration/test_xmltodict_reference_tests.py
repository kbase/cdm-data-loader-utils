"""Integration tests for the xmltodict pipeline against the uniref reference dataset."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import dlt
import duckdb
import pytest

from tests.cdm_data_loaders.pipelines.helpers import (
    BUFFER_SIZES,
    LOADER_FILE_FORMATS,
    WORKER_CONFIGS,
    assert_dataset_matches_reference,
    assert_datasets_equal,
    read_pipeline_tables,
    reconstruct_entries,
    reference_entries_from_duckdb,
    run_and_read,
    sorted_json,
)

CHUNK_DIRS = ("chunk_5_el", "chunk_20_el", "chunk_100_el")
EXPECTED_ENTRY_COUNT = 100


# tests for the "reference" data structure
def test_xmltodict_reference_xml_as_jsonl_pass_matches_source(
    reference_xml_entries: list[dict[str, Any]],
    reference_xml_as_jsonl: Path,
) -> None:
    """The generated reference JSONL has 100 entries with the expected shape."""
    assert len(reference_xml_entries) == EXPECTED_ENTRY_COUNT
    with reference_xml_as_jsonl.open(encoding="utf-8") as fh:
        line_count = sum(1 for line in fh if line.strip())
    assert line_count == EXPECTED_ENTRY_COUNT
    for key in ("_id", "name", "property", "representativeMember"):
        assert key in reference_xml_entries[0].get("entry", {})


def test_xml_pass_reference_dataset_has_nested_structures(
    reference_xml_dataset: duckdb.DuckDBPyConnection,
) -> None:
    """The DuckDB reference stores nested JSON as STRUCT/LIST types, not JSON text columns."""
    reference_xml_table = reference_xml_dataset.table("reference")
    entry_type = reference_xml_table.types[0]
    entry_fields = dict(entry_type.children)

    property_type = entry_fields["property"]
    member_type = entry_fields["representativeMember"]
    assert property_type.id == "list"
    assert member_type.id == "struct"

    count = reference_xml_table.count("*").fetchone()[0]
    assert count == EXPECTED_ENTRY_COUNT


def test_xml_pass_diff_format_same_results(run_xmltodict_pipeline_fxt: Callable[..., Any]) -> None:
    """Ensure that the save format does not affect the output."""
    load_info_and_dir = {
        str(output_format): run_and_read(
            run_xmltodict_pipeline_fxt,
            "chunk_20_el",
            loader_file_format=output_format,
        )
        for output_format in LOADER_FILE_FORMATS
    }
    assert_datasets_equal(load_info_and_dir)


@pytest.mark.xfail(reason="not yet implemented")
@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xml_pass_nesting_levels(
    run_xmltodict_pipeline_fxt: Callable[..., Any],
    sorted_reference_xml_entries: list[str],
    loader_file_format: str,
) -> None:
    """Pipeline output matches the reference after reconstruction."""
    load_info_and_dir = {
        f"nesting_{option}": run_and_read(
            run_xmltodict_pipeline_fxt,
            "chunk_20_el",
            loader_file_format=loader_file_format,
            preserve_table_nesting=option,
        )
        for option in [True, False]
    }

    # if preserve_table_nesting is True, there should only be a single table
    datasets = {key: value[0].pipeline.dataset() for key, value in load_info_and_dir.items()}
    assert {t for t in datasets["nesting_True"].tables if not t.startswith("_dlt")} == {"entry"}
    nesting_false_tables = {t for t in datasets["nesting_False"].tables if not t.startswith("_dlt")}

    assert "entry" in nesting_false_tables
    assert len(nesting_false_tables) > 1

    assert_datasets_equal(load_info_and_dir)
    assert_dataset_matches_reference(load_info_and_dir, sorted_reference_xml_entries, loader_file_format)


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xml_pass_output_identical_across_chunk_dirs(
    run_xmltodict_pipeline_fxt: Callable[..., Any], sorted_reference_xml_entries: list[str], loader_file_format: str
) -> None:
    """Pipeline output is identical across all three chunk dirs and matches the reference."""
    load_info_and_dir = {
        f"{key}_{loader_file_format}": run_and_read(
            run_xmltodict_pipeline_fxt, key, loader_file_format=loader_file_format
        )
        for key in CHUNK_DIRS
    }

    assert_datasets_equal(load_info_and_dir)
    assert_dataset_matches_reference(load_info_and_dir, sorted_reference_xml_entries, loader_file_format)


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xml_pass_output_identical_across_buffer_sizes(
    run_xmltodict_pipeline_fxt: Callable[..., Any], sorted_reference_xml_entries: list[str], loader_file_format: str
) -> None:
    """List buffer size does not affect pipeline output for any chunk dir."""
    for chunk_dir in CHUNK_DIRS:
        load_info_and_dir = {
            f"{buffer_size}_{loader_file_format}": run_and_read(
                run_xmltodict_pipeline_fxt, chunk_dir, buffer_size=buffer_size, loader_file_format=loader_file_format
            )
            for buffer_size in BUFFER_SIZES
        }
        assert_datasets_equal(load_info_and_dir)
        assert_dataset_matches_reference(load_info_and_dir, sorted_reference_xml_entries, loader_file_format)


@pytest.mark.parametrize(
    ("extract_workers", "normalize_workers", "load_workers"),
    WORKER_CONFIGS,
    ids=[f"extract_{e}-normalize_{n}-load_{l}" for e, n, l in WORKER_CONFIGS],  # noqa: E741
)
@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_xml_pass_no_data_lost_across_worker_configs(
    *,
    run_xmltodict_pipeline_fxt: Callable[..., Any],
    sorted_reference_xml_entries: list[str],
    extract_workers: int,
    normalize_workers: int,
    load_workers: int,
    loader_file_format: str,
) -> None:
    """No data is lost regardless of how the job is split across extract/normalize/load workers."""
    load_info_and_dir = {}
    for worker_key in ("baseline", "current"):
        worker_config = (
            {}
            if worker_key == "baseline"
            else {
                "extract.workers": extract_workers,
                "normalize.workers": normalize_workers,
                "load.workers": load_workers,
            }
        )
        with dlt.config.values(worker_config):
            load_info_and_dir[worker_key] = run_and_read(
                run_xmltodict_pipeline_fxt, "chunk_5_el", loader_file_format=loader_file_format
            )

    assert_datasets_equal(load_info_and_dir)
    assert_dataset_matches_reference(load_info_and_dir, sorted_reference_xml_entries, loader_file_format)


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xml_pass_reconstructed_matches_duckdb_reference(
    run_xmltodict_pipeline_fxt: Callable[..., Any],
    reference_xml_dataset: duckdb.DuckDBPyConnection,
    loader_file_format: str,
) -> None:
    """Reconstructed pipeline entries equal the entries loaded from the DuckDB reference."""
    (load_info, output_dir) = run_and_read(
        run_xmltodict_pipeline_fxt, "chunk_100_el", loader_file_format=loader_file_format
    )
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    reconstructed = reconstruct_entries(tables)
    from_duckdb = reference_entries_from_duckdb(reference_xml_dataset)
    assert len(from_duckdb) == EXPECTED_ENTRY_COUNT
    assert sorted_json(reconstructed) == sorted_json(from_duckdb)
