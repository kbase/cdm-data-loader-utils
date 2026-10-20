"""Shared fixtures for jsonlines pipeline tests."""

import json
from collections.abc import Callable, Generator
from pathlib import Path
from typing import Any

import duckdb
import pytest
from dlt.common.pipeline import LoadInfo
from pydantic import BaseModel

from cdm_data_loaders.pipelines.jsonlines.extract_jsonschema_validate_pipeline import (
    run_jsonlines_pipeline as run_jsonschema_pipeline,
)
from cdm_data_loaders.pipelines.jsonlines.extract_pipeline import run_jsonlines_pipeline
from cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline import (
    run_jsonlines_with_validation_pipeline,
)
from cdm_data_loaders.pipelines.jsonlines.settings import (
    JsonlExtractSettings,
    JsonlJsonschemaIngestSettings,
    JsonlPydanticIngestSettings,
)
from tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_models import ENTITY_MODELS
from tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_schemas import ENTITY_SCHEMAS
from tests.cdm_data_loaders.pipelines.jsonlines.integration.jsonlines_reference_helpers import (
    canonicalize_reference_entry,
)
from tests.conftest import TEST_DATA_DIR

type PipelineRunner = Callable[..., tuple[LoadInfo | None, Path]]
type PrepareReference = Callable[[dict[str, Any]], dict[str, Any]]


@pytest.fixture(params=["extract", "extract_validate", "extract_jsonschema_validate"], ids=str)
def pipeline_kind(request: pytest.FixtureRequest) -> str:
    """Run the same reference contract against all three JSONL pipelines."""
    return request.param


@pytest.fixture
def run_reference_pipeline(
    pipeline_kind: str,
    reference_data_dir: Path,
    dlt_destination_config: str,
    isolated_run_factory: Callable[..., Callable[..., tuple[LoadInfo | None, Path]]],
) -> PipelineRunner:
    """Run the selected JSONL pipeline on a chunk dir with per-run isolated output directories."""
    settings_cls, registry_kwarg, dataset_prefix, file_glob = {
        "extract": (JsonlExtractSettings, {"table_name": "dataset"}, "jsonl_extract_run", "**/*.jsonl*"),
        "extract_validate": (
            JsonlPydanticIngestSettings,
            {"entity_models_module": "tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_models"},
            "jsonl_extract_validate_run",
            "*.jsonl*",
        ),
        "extract_jsonschema_validate": (
            JsonlJsonschemaIngestSettings,
            {"schema_files_module": "tests.cdm_data_loaders.pipelines.jsonlines.integration.dataset_report_schemas"},
            "jsonschema_run",
            "*.jsonl*",
        ),
    }[pipeline_kind]
    run_fn = {
        "extract": run_jsonlines_pipeline,
        "extract_validate": run_jsonlines_with_validation_pipeline,
        "extract_jsonschema_validate": run_jsonschema_pipeline,
    }[pipeline_kind]
    base_kwargs: dict[str, Any] = {
        "buffer_size": 10,
        "dlt_dev_mode": False,
        "file_glob": file_glob,
        "settings_cls": settings_cls,
        "use_destination": dlt_destination_config,
        "use_output_dir_for_pipeline_metadata": True,
        **registry_kwarg,
    }
    runner = isolated_run_factory(dataset_prefix=dataset_prefix, run_fn=run_fn, base_kwargs=base_kwargs)

    def _factory(chunk_dir: str, **overrides: Any) -> tuple[LoadInfo | None, Path]:
        kwargs = {"input_dir": str(reference_data_dir / chunk_dir), **overrides}
        return runner(**kwargs)

    return _factory


@pytest.fixture
def prepare_reference(pipeline_kind: str) -> PrepareReference:
    """Apply Pydantic defaults only for model validation; other pipelines preserve input."""

    def _prepare(entry: dict[str, Any]) -> dict[str, Any]:
        if pipeline_kind == "extract_validate":
            return ENTITY_MODELS["dataset"].model_validate(entry).model_dump(mode="json")
        return entry

    return _prepare


@pytest.fixture
def reference_registry_factory(
    pipeline_kind: str,
    entity_models_module_factory: Callable[[dict[str, type[BaseModel]]], str],
    schema_module_factory: Callable[[object], str],
) -> Callable[[list[str]], dict[str, str]]:
    """Register the requested tables using the selected pipeline's assembly contract."""

    def _factory(table_names: list[str]) -> dict[str, str]:
        if pipeline_kind == "extract_validate":
            return {
                "entity_models_module": entity_models_module_factory(
                    dict.fromkeys(table_names, ENTITY_MODELS["dataset"])
                )
            }
        if pipeline_kind == "extract_jsonschema_validate":
            return {"schema_files_module": schema_module_factory(dict.fromkeys(table_names, ENTITY_SCHEMAS["dataset"]))}
        return {}

    return _factory


@pytest.fixture(scope="session")
def jsonl_reference_entities(reference_jsonl: Path) -> list[dict[str, Any]]:
    """Parsed JSONL file."""
    with reference_jsonl.open(encoding="utf-8") as fh:
        return [json.loads(line.strip()) for line in fh if line.strip()]


@pytest.fixture(scope="session")
def reference_jsonl() -> Path:
    """Path to the reference JSONL file."""
    return TEST_DATA_DIR / "ncbi_rest_api" / "dataset_report" / "chunk_52" / "dataset" / "dataset_reports_52.jsonl"


@pytest.fixture
def reference_data_dir() -> Path:
    """Return the directory containing all reference chunk configurations."""
    return TEST_DATA_DIR / "ncbi_rest_api" / "dataset_report"


@pytest.fixture(scope="session")
def reference_dataset(reference_jsonl: Path) -> Generator[duckdb.DuckDBPyConnection]:
    """Load the reference JSONL into an in-memory DuckDB connection with nested structures.

    read_json_auto infers native STRUCT/STRUCT[] types for the nested content.
    """
    connection = duckdb.connect()
    jsonl_path = str(reference_jsonl)
    connection.execute(
        "CREATE TABLE reference AS SELECT * FROM read_json_auto(?, sample_size=200)",
        [jsonl_path],
    )
    try:
        yield connection
    finally:
        connection.close()


@pytest.fixture
def sorted_reference_entities(
    jsonl_reference_entities: list[dict[str, Any]], prepare_reference: PrepareReference
) -> list[str]:
    """Canonicalized entries as sorted JSON strings for order-insensitive comparison."""
    lines_as_json = [
        json.dumps(canonicalize_reference_entry(prepare_reference(obj)), sort_keys=True, default=str)
        for obj in jsonl_reference_entities
    ]
    return sorted(lines_as_json)


@pytest.fixture
def settings_factory(
    tmp_path: Path, dlt_destination_config: str, entity_models_module: str
) -> Callable[..., JsonlPydanticIngestSettings]:
    """Return a factory that builds a valid JsonlPydanticIngestSettings, with fields open to override."""

    def _factory(**overrides: Any) -> JsonlPydanticIngestSettings:
        log_config_file = tmp_path / "logging.json"
        log_config_file.write_text('{"version": 1}')
        input_dir = tmp_path / "input"
        input_dir.mkdir(exist_ok=True)
        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        kwargs: dict[str, Any] = {
            "dataset_name": "test_dataset",
            "dlt_dev_mode": False,
            "entity_models_module": entity_models_module,
            "input_dir": str(input_dir),
            "log_config_file": str(log_config_file),
            "output_dir": str(output_dir),
            "use_destination": dlt_destination_config,
            "use_output_dir_for_pipeline_metadata": False,
        }
        kwargs.update(overrides)
        return JsonlPydanticIngestSettings(**kwargs)

    return _factory


@pytest.fixture
def write_gzip_jsonl_file(write_gzip_file: Callable[[Path, str, str], Path]) -> Callable[[Path, str, list[str]], Path]:
    """Return a function that writes a list of JSON lines to a gzip-compressed file."""

    def _write(directory: Path, filename: str, lines: list[str]) -> Path:
        return write_gzip_file(directory, filename, "\n".join(lines) + "\n")

    return _write
