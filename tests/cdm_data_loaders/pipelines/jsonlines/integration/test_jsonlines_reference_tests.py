"""Integration tests for all JSONL pipelines against the dataset_report reference dataset."""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import dlt
import duckdb
import pytest
from jsonschema import FormatChecker
from jsonschema.validators import validator_for

from tests.cdm_data_loaders.pipelines.helpers import (
    BUFFER_SIZES,
    LOADER_FILE_FORMATS,
    WORKER_CONFIGS,
    assert_dataset_matches_reference,
    assert_datasets_equal,
    read_pipeline_tables,
    run_and_read,
)
from tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_models import ENTITY_MODELS
from tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_schemas import ENTITY_SCHEMAS
from tests.cdm_data_loaders.pipelines.jsonlines.integration.jsonlines_reference_helpers import (
    EXPECTED_ENTRY_COUNT,
    canonicalize_pipeline_entry,
    canonicalize_reference_entry,
    reconstructed_entries,
    reference_entries_from_duckdb,
)

CHUNK_DIRS = ("chunk_4", "chunk_10", "chunk_52")


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
@pytest.mark.parametrize("pipeline_kind", ["extract_validate", "extract_jsonschema_validate"], ids=str)
def test_jsonlines_validation_pass_reference_records_load(
    run_reference_pipeline: Callable[..., Any], loader_file_format: str
) -> None:
    """All assembly reports load into the validated dataset without rejections."""
    load_info, output_dir = run_and_read(run_reference_pipeline, "chunk_52", loader_file_format=loader_file_format)
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    assert "dataset_rejected" not in tables
    assert len(tables["dataset"]) == EXPECTED_ENTRY_COUNT


def test_jsonlines_reference_pass_assembly_model_accepts_all_records(
    jsonl_reference_entities: list[dict[str, Any]],
) -> None:
    """The registered assembly model accepts every reference record without losing accessions."""
    model = ENTITY_MODELS["dataset"]
    validated = [model.model_validate(entry) for entry in jsonl_reference_entities]
    assert len(validated) == EXPECTED_ENTRY_COUNT
    assert [entry.model_dump()["accession"] for entry in validated] == [
        entry["accession"] for entry in jsonl_reference_entities
    ]
    assert [canonicalize_reference_entry(entry.model_dump(mode="json", exclude_unset=True)) for entry in validated] == [
        canonicalize_reference_entry(entry) for entry in jsonl_reference_entities
    ]


def test_jsonlines_reference_pass_assembly_schema_accepts_all_records(
    jsonl_reference_entities: list[dict[str, Any]],
) -> None:
    """The attached schema validates every individual report with all nested references resolved."""
    schema = ENTITY_SCHEMAS["dataset"]
    validator_class = validator_for(schema)
    validator_class.check_schema(schema)
    validator = validator_class(schema, format_checker=FormatChecker())
    assert len(jsonl_reference_entities) == EXPECTED_ENTRY_COUNT
    assert [list(validator.iter_errors(entry)) for entry in jsonl_reference_entities] == [
        [] for _ in jsonl_reference_entities
    ]


# tests for the "reference" data structure
def test_jsonlines_reference_jsonl_pass_matches_source(
    jsonl_reference_entities: list[dict[str, Any]],
    reference_jsonl: Path,
) -> None:
    """The generated reference JSONL has all expected entries in the expected structure."""
    assert len(jsonl_reference_entities) == EXPECTED_ENTRY_COUNT
    with reference_jsonl.open(encoding="utf-8") as fh:
        line_count = sum(1 for line in fh if line.strip())
    assert line_count == EXPECTED_ENTRY_COUNT
    for tag in ("accession", "current_accession", "paired_accession", "source_database", "organism"):
        assert tag in jsonl_reference_entities[0]


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
@pytest.mark.parametrize("source_dir", CHUNK_DIRS, ids=str)
def test_jsonlines_passed_dataset_matches_reference(
    run_reference_pipeline: Callable[..., Any],
    source_dir: str,
    loader_file_format: str,
    sorted_reference_entities: list[str],
) -> None:
    """Pipeline output matches the pipeline input."""
    load_info_and_dir = {
        source_dir: run_and_read(run_reference_pipeline, source_dir, loader_file_format=loader_file_format)
    }

    assert_dataset_matches_reference(
        load_info_and_dir,
        sorted_reference_entities,
        loader_file_format,
        reconstruct=reconstructed_entries,
        canonicalize=canonicalize_pipeline_entry,
    )


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_jsonlines_pass_output_identical_across_chunk_dirs(
    run_reference_pipeline: Callable[..., Any], sorted_reference_entities: list[str], loader_file_format: str
) -> None:
    """Pipeline output is identical across all three chunk dirs and matches the reference."""
    load_info_and_dir = {
        f"{key}_{loader_file_format}": run_and_read(run_reference_pipeline, key, loader_file_format=loader_file_format)
        for key in CHUNK_DIRS
    }

    assert_datasets_equal(load_info_and_dir)
    assert_dataset_matches_reference(
        load_info_and_dir,
        sorted_reference_entities,
        loader_file_format,
        reconstruct=reconstructed_entries,
        canonicalize=canonicalize_pipeline_entry,
    )


def test_jsonlines_pass_reference_dataset_has_nested_structures(
    reference_dataset: duckdb.DuckDBPyConnection,
) -> None:
    """The DuckDB reference stores nested JSON as STRUCT/LIST types, not JSON text columns."""
    column_types = dict(
        reference_dataset.execute(
            "SELECT column_name, data_type FROM information_schema.columns WHERE table_name = 'reference'"
        ).fetchall()
    )

    organism_type = column_types["organism"]
    assert organism_type.startswith("STRUCT")
    assembly_stats_type = column_types["assembly_stats"]
    assert assembly_stats_type.startswith("STRUCT")
    additional_submitters_type = column_types["additional_submitters"]
    assert additional_submitters_type.startswith("STRUCT(")
    assert additional_submitters_type.endswith("[]")

    count = reference_dataset.execute("SELECT COUNT(*) FROM reference").fetchone()
    assert count == (EXPECTED_ENTRY_COUNT,)


def test_jsonlines_pass_diff_format_same_results(
    run_reference_pipeline: Callable[..., Any], sorted_reference_entities: list[str]
) -> None:
    """Ensure that the save format does not affect the output."""
    load_info_and_dir = {
        str(output_format): run_and_read(
            run_reference_pipeline,
            "chunk_10",
            loader_file_format=output_format,
        )
        for output_format in LOADER_FILE_FORMATS
    }
    for output_format, result in load_info_and_dir.items():
        assert_dataset_matches_reference(
            {output_format: result},
            sorted_reference_entities,
            output_format,
            reconstruct=reconstructed_entries,
            canonicalize=canonicalize_pipeline_entry,
        )


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_jsonlines_pass_output_identical_across_buffer_sizes(
    run_reference_pipeline: Callable[..., Any], sorted_reference_entities: list[str], loader_file_format: str
) -> None:
    """List buffer size does not affect pipeline output for any chunk dir."""
    for chunk_dir in CHUNK_DIRS:
        load_info_and_dir = {
            f"{buffer_size}_{loader_file_format}": run_and_read(
                run_reference_pipeline, chunk_dir, buffer_size=buffer_size, loader_file_format=loader_file_format
            )
            for buffer_size in BUFFER_SIZES
        }

        assert_datasets_equal(load_info_and_dir)
        assert_dataset_matches_reference(
            load_info_and_dir,
            sorted_reference_entities,
            loader_file_format,
            reconstruct=reconstructed_entries,
            canonicalize=canonicalize_pipeline_entry,
        )


@pytest.mark.parametrize(
    ("extract_workers", "normalize_workers", "load_workers"),
    WORKER_CONFIGS,
    ids=[f"extract_{extract}-normalize_{normalize}-load_{load}" for extract, normalize, load in WORKER_CONFIGS],
)
@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_jsonlines_pass_no_data_lost_across_worker_configs(
    *,
    run_reference_pipeline: Callable[..., Any],
    sorted_reference_entities: list[str],
    extract_workers: int,
    normalize_workers: int,
    load_workers: int,
    loader_file_format: str,
) -> None:
    """No data is lost regardless of how the job is split across extract/normalize/load workers."""
    load_info_and_dir = {}
    with dlt.config.values(
        {
            "extract.workers": extract_workers,
            "normalize.workers": normalize_workers,
            "load.workers": load_workers,
        }
    ):
        load_info_and_dir["chunk_4"] = run_and_read(
            run_reference_pipeline, "chunk_4", loader_file_format=loader_file_format
        )

    assert_dataset_matches_reference(
        load_info_and_dir,
        sorted_reference_entities,
        output_format=loader_file_format,
        reconstruct=reconstructed_entries,
        canonicalize=canonicalize_pipeline_entry,
    )


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_jsonlines_pass_reconstructed_matches_duckdb_reference(
    run_reference_pipeline: Callable[..., Any],
    reference_dataset: duckdb.DuckDBPyConnection,
    loader_file_format: str,
    prepare_reference: Callable[[dict[str, Any]], dict[str, Any]],
) -> None:
    """Reconstructed pipeline entries equal the entries loaded from the DuckDB reference."""
    (load_info, output_dir) = run_and_read(run_reference_pipeline, "chunk_52", loader_file_format=loader_file_format)
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    reconstructed = reconstructed_entries(tables)
    from_duckdb = reference_entries_from_duckdb(reference_dataset)
    assert len(from_duckdb) == EXPECTED_ENTRY_COUNT
    from_duckdb_canonical = sorted(
        json.dumps(canonicalize_reference_entry(prepare_reference(entry)), sort_keys=True, default=str)
        for entry in from_duckdb
    )
    from_pipeline = sorted(
        json.dumps(canonicalize_pipeline_entry(entry), sort_keys=True, default=str) for entry in reconstructed
    )
    assert from_pipeline == from_duckdb_canonical


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_jsonlines_pass_gzip_input_matches_plain_text(
    *,
    run_reference_pipeline: Callable[..., Any],
    tmp_path: Path,
    write_gzip_jsonl_file: Callable[[Path, str, list[str]], Path],
    jsonl_reference_entities: list[dict[str, Any]],
    sorted_reference_entities: list[str],
    loader_file_format: str,
) -> None:
    """A gzip-compressed input file produces the same rows as the plain-text reference."""
    input_dir = tmp_path / "input_gzip"
    table_dir = input_dir / "dataset"
    write_gzip_jsonl_file(table_dir, "compressed.jsonl.gz", [json.dumps(entity) for entity in jsonl_reference_entities])

    result = run_and_read(
        run_reference_pipeline, "chunk_52", input_dir=str(input_dir), loader_file_format=loader_file_format
    )
    assert_dataset_matches_reference(
        {"gzip": result},
        sorted_reference_entities,
        loader_file_format,
        reconstruct=reconstructed_entries,
        canonicalize=canonicalize_pipeline_entry,
    )


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_jsonlines_pass_nesting_matches_pipeline_contract(
    run_reference_pipeline: Callable[..., Any],
    pipeline_kind: str,
    loader_file_format: str,
    sorted_reference_entities: list[str],
) -> None:
    """Extraction flattens nested fields; validation preserves them without changing their content."""
    result = run_and_read(run_reference_pipeline, "chunk_10", loader_file_format=loader_file_format)
    load_info, output_dir = result
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    data_tables = {name for name in tables if not name.startswith("_dlt")}
    columns = set().union(*(row.keys() for row in tables["dataset"]))
    if pipeline_kind != "extract":
        assert data_tables == {"dataset"}
        assert {"organism", "assembly_info", "additional_submitters"} <= columns
    else:
        assert {"dataset", "dataset__additional_submitters"} <= data_tables
        assert "organism__tax_id" in columns
    assert_dataset_matches_reference(
        {pipeline_kind: result},
        sorted_reference_entities,
        loader_file_format,
        reconstruct=reconstructed_entries,
        canonicalize=canonicalize_pipeline_entry,
    )
