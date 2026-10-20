"""Error and boundary tests for all JSONL pipelines using assembly-report contracts."""

import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.pipelines.jsonlines.extract_jsonschema_validate_pipeline import cli as jsonschema_cli
from cdm_data_loaders.pipelines.jsonlines.extract_pipeline import cli as extract_cli
from cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline import cli as validate_cli
from tests.cdm_data_loaders.pipelines.helpers import LOADER_FILE_FORMATS, read_pipeline_tables, run_and_read
from tests.cdm_data_loaders.pipelines.jsonlines.integration.jsonlines_reference_helpers import (
    canonicalize_pipeline_entry,
    canonicalize_reference_entry,
    reconstructed_entries,
)

pytestmark = pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)


@pytest.mark.parametrize("include_valid", [True, False], ids=["mixed", "all-invalid"])
def test_jsonlines_fail_malformed_rows_preserve_provenance(
    *,
    run_reference_pipeline: Callable[..., Any],
    pipeline_kind: str,
    tmp_path: Path,
    loader_file_format: str,
    include_valid: bool,
    jsonl_reference_entities: list[dict[str, Any]],
    prepare_reference: Callable[[dict[str, Any]], dict[str, Any]],
) -> None:
    """Malformed rows retain physical line numbers and file paths across files and buffer flushes."""
    input_dir = tmp_path / "mixed"
    table_dir = input_dir / "dataset"
    table_dir.mkdir(parents=True)
    entries = jsonl_reference_entities[:3] if include_valid else []
    valid_lines = [json.dumps(entry) for entry in entries]
    malformed = "{not valid json"
    (table_dir / "first.jsonl").write_text("\n \t\n" + malformed + "\n" + "\n".join(valid_lines), encoding="utf-8")
    (table_dir / "second.jsonl").write_text(malformed, encoding="utf-8")
    result, output_dir = run_and_read(
        run_reference_pipeline,
        "chunk_52",
        input_dir=str(input_dir),
        buffer_size=2,
        loader_file_format=loader_file_format,
    )
    tables = read_pipeline_tables(output_dir / result.dataset_name, loader_file_format)
    actual = reconstructed_entries(tables, "dataset")
    assert sorted(
        json.dumps(canonicalize_pipeline_entry(entry), sort_keys=True, default=str) for entry in actual
    ) == sorted(
        json.dumps(canonicalize_reference_entry(prepare_reference(entry)), sort_keys=True, default=str)
        for entry in entries
    )
    suffix, error_key = ("invalid", "parse_error") if pipeline_kind == "extract" else ("rejected", "error_detail")
    rejected = sorted(tables[f"dataset_{suffix}"], key=lambda row: row["source_file"])
    prefix = "dataset/" if pipeline_kind == "extract" else ""
    with pytest.raises(json.JSONDecodeError) as error_info:
        json.loads(malformed)
    expected_error = str(error_info.value)
    assert [
        (row["source_file"], row["line_no"], row["raw_record"], row.get("record"), row[error_key]) for row in rejected
    ] == [
        (f"{prefix}first.jsonl", 3, malformed, None, expected_error),
        (f"{prefix}second.jsonl", 1, malformed, None, expected_error),
    ]
    assert {name for name in tables if "__" not in name and not name.startswith("_dlt")} == (
        {"dataset", f"dataset_{suffix}"} if include_valid else {f"dataset_{suffix}"}
    )


@pytest.mark.parametrize("input_case", ["missing-entity", "empty-directory", "empty-file", "blank-lines"], ids=str)
def test_jsonlines_pass_empty_inputs_produce_no_data(
    run_reference_pipeline: Callable[..., Any], tmp_path: Path, loader_file_format: str, input_case: str
) -> None:
    """Missing entity directories, no files, empty files, and whitespace yield no data rows."""
    input_dir = tmp_path / "empty"
    input_dir.mkdir()
    if input_case != "missing-entity":
        table_dir = input_dir / "dataset"
        table_dir.mkdir()
        if input_case in {"empty-file", "blank-lines"}:
            content = "\n \t\n\r\n" if input_case == "blank-lines" else ""
            (table_dir / "empty.jsonl").write_text(content, encoding="utf-8")
    load_info, output_dir = run_and_read(
        run_reference_pipeline, "chunk_52", input_dir=str(input_dir), loader_file_format=loader_file_format
    )
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    assert {name: rows for name, rows in tables.items() if not name.startswith("_dlt")} == {}


def test_jsonlines_pass_duplicates_blanks_and_final_line_preserved(
    *,
    run_reference_pipeline: Callable[..., Any],
    tmp_path: Path,
    loader_file_format: str,
    jsonl_reference_entities: list[dict[str, Any]],
    prepare_reference: Callable[[dict[str, Any]], dict[str, Any]],
) -> None:
    """Blank lines are skipped without deduplicating records or losing an unterminated final line."""
    input_dir = tmp_path / "duplicates"
    table_dir = input_dir / "dataset"
    table_dir.mkdir(parents=True)
    entry = jsonl_reference_entities[0]
    line = json.dumps(entry)
    (table_dir / "records.jsonl").write_text(f"\n{line}\n\t \n{line}", encoding="utf-8")
    load_info, output_dir = run_and_read(
        run_reference_pipeline,
        "chunk_52",
        input_dir=str(input_dir),
        buffer_size=1,
        loader_file_format=loader_file_format,
    )
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    expected = canonicalize_reference_entry(prepare_reference(entry))
    assert [canonicalize_pipeline_entry(row) for row in reconstructed_entries(tables)] == [expected, expected]


def test_jsonlines_pass_multiple_entities_and_file_glob(
    *,
    run_reference_pipeline: Callable[..., Any],
    pipeline_kind: str,
    tmp_path: Path,
    loader_file_format: str,
    reference_registry_factory: Callable[[list[str]], dict[str, str]],
) -> None:
    """Matching files route to configured tables and a custom glob excludes malformed files."""
    input_dir = tmp_path / "entities"
    expected = {"dataset": "GCA_000003115.1", "other": "GCF_002271065.1"}
    for table_name, accession in expected.items():
        table_dir = input_dir / table_name
        table_dir.mkdir(parents=True)
        (table_dir / "keep.jsonl").write_text(json.dumps({"accession": accession}), encoding="utf-8")
        (table_dir / "ignore.jsonl").write_text("{bad json", encoding="utf-8")
    load_info, output_dir = run_and_read(
        run_reference_pipeline,
        "chunk_52",
        input_dir=str(input_dir),
        file_glob="**/keep.jsonl",
        loader_file_format=loader_file_format,
        **reference_registry_factory(list(expected)),
    )
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    expected_tables = (
        {"dataset": sorted(expected.values())}
        if pipeline_kind == "extract"
        else {table_name: [accession] for table_name, accession in expected.items()}
    )
    assert {
        name: sorted(row["accession"] for row in rows) for name, rows in tables.items() if not name.startswith("_dlt")
    } == expected_tables


@pytest.mark.parametrize(
    ("record", "error_type", "location", "error_message"),
    [
        pytest.param(
            {"source_database": "not-a-database"},
            "enum",
            ["source_database"],
            "'not-a-database' is not one of ['SOURCE_DATABASE_UNSPECIFIED', "
            "'SOURCE_DATABASE_GENBANK', 'SOURCE_DATABASE_REFSEQ']",
            id="invalid-enum",
        ),
        pytest.param(
            {"organism": {"tax_id": "not-an-integer"}},
            "int_parsing",
            ["organism", "tax_id"],
            "'not-an-integer' is not of type 'integer'",
            id="nested-type",
        ),
        pytest.param(None, "model_type", [], "None is not of type 'object'", id="json-null"),
        pytest.param([], "model_type", [], "[] is not of type 'object'", id="json-array"),
        pytest.param("text", "model_type", [], "'text' is not of type 'object'", id="json-string"),
    ],
)
@pytest.mark.parametrize("pipeline_kind", ["extract_validate", "extract_jsonschema_validate"], ids=str)
def test_jsonlines_validation_fail_schema_invalid_rows_are_rejected(
    *,
    run_reference_pipeline: Callable[..., Any],
    pipeline_kind: str,
    tmp_path: Path,
    loader_file_format: str,
    record: Any,
    error_type: str,
    location: list[str],
    error_message: str,
) -> None:
    """Both validators preserve rejected input and their error details without dropping valid neighbors."""
    input_dir = tmp_path / "invalid_model"
    table_dir = input_dir / "dataset"
    table_dir.mkdir(parents=True)
    valid = {"accession": "GCA_000003115.1"}
    raw_record = json.dumps(record)
    (table_dir / "records.jsonl").write_text(
        json.dumps(valid) + "\n" + raw_record + "\n" + json.dumps(valid), encoding="utf-8"
    )
    load_info, output_dir = run_and_read(
        run_reference_pipeline,
        "chunk_52",
        input_dir=str(input_dir),
        buffer_size=2,
        loader_file_format=loader_file_format,
    )
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    assert [row["accession"] for row in tables["dataset"]] == [valid["accession"], valid["accession"]]
    rejected = tables["dataset_rejected"]
    assert len(rejected) == 1
    row = rejected[0]
    assert (row["source_file"], row["line_no"]) == ("records.jsonl", 2)
    assert row["raw_record"] == raw_record
    errors = json.loads(row["error_detail"])
    if pipeline_kind == "extract_validate":
        assert reconstructed_entries(tables, "dataset_rejected")[0].get("record") == record
        assert len(errors) == 1
        assert errors[0]["type"] == error_type
        assert errors[0]["loc"] == location
    else:
        assert "record" not in row
        assert errors == [error_message]


@pytest.mark.parametrize("pipeline_kind", ["extract_validate", "extract_jsonschema_validate"], ids=str)
def test_jsonlines_validation_pass_table_selection_skips_other_entities(
    *,
    run_reference_pipeline: Callable[..., Any],
    reference_registry_factory: Callable[[list[str]], dict[str, str]],
    tmp_path: Path,
    loader_file_format: str,
) -> None:
    """An explicit table selection processes only that registered entity and uses the configured dataset name."""
    input_dir = tmp_path / "selected"
    for table_name in ("dataset", "other"):
        table_dir = input_dir / table_name
        table_dir.mkdir(parents=True)
        (table_dir / "records.jsonl").write_text('{"accession": "GCA_000003115.1"}', encoding="utf-8")
    load_info, output_dir = run_and_read(
        run_reference_pipeline,
        "chunk_52",
        input_dir=str(input_dir),
        table_names=["other"],
        dataset_name="selected",
        loader_file_format=loader_file_format,
        **reference_registry_factory(["dataset", "other"]),
    )
    assert load_info.dataset_name == "selected"
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    assert {name for name in tables if not name.startswith("_dlt")} == {"other"}
    assert [row["accession"] for row in tables["other"]] == ["GCA_000003115.1"]


def test_jsonlines_cli_pass_reference_dataset_loads(
    *,
    pipeline_kind: str,
    reference_data_dir: Path,
    tmp_path: Path,
    loader_file_format: str,
    dlt_destination_config: str,
    monkeypatch: pytest.MonkeyPatch,
    sorted_reference_entities: list[str],
) -> None:
    """All real CLI entry points load the reference dataset from command-line settings."""
    log_config = tmp_path / "logging.json"
    log_config.write_text('{"version": 1}', encoding="utf-8")
    output_dir = tmp_path / "bucket"
    registry_flag, registry_module = {
        "extract": ("table-name", "dataset"),
        "extract_validate": (
            "entity-models-module",
            "tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_models",
        ),
        "extract_jsonschema_validate": (
            "schema-files-module",
            "tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_schemas",
        ),
    }[pipeline_kind]
    argv = [
        "jsonlines",
        "--input-dir",
        str(reference_data_dir / "chunk_52"),
        "--output-dir",
        str(output_dir),
        "--log-config-file",
        str(log_config),
        "--use-destination",
        dlt_destination_config,
        "--dataset-name",
        "reference_cli",
        "--file-glob",
        "**/*.jsonl",
        "--loader-file-format",
        str(loader_file_format),
        "--use-output-dir-for-pipeline-metadata",
        "true",
    ]
    if registry_flag:
        argv.extend([f"--{registry_flag}", registry_module])
    monkeypatch.setattr(sys, "argv", argv)
    cli = {
        "extract": extract_cli,
        "extract_validate": validate_cli,
        "extract_jsonschema_validate": jsonschema_cli,
    }[pipeline_kind]
    load_info = cli.run()
    assert load_info is not None
    assert not load_info.has_failed_jobs
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    assert (
        sorted(
            json.dumps(canonicalize_pipeline_entry(entry), sort_keys=True, default=str)
            for entry in reconstructed_entries(tables)
        )
        == sorted_reference_entities
    )
