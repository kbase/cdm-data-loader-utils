"""Integration tests for build_entity_resource: file reading, validation, and dlt wiring together."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import dlt
from pydantic import BaseModel

from cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline import build_entity_resource
from cdm_data_loaders.pipelines.jsonlines.settings import JsonlPydanticIngestSettings
from tests.conftest import LOCAL_FS


def run_pipeline(table_name: str, model: type[BaseModel], settings: JsonlPydanticIngestSettings) -> Any:
    """Run one entity resource through a local dlt pipeline. Return the resulting dataset."""
    resource = build_entity_resource(table_name, model, settings)

    pipeline = dlt.pipeline(
        pipeline_name=f"test_{table_name}",
        destination=dlt.destination(LOCAL_FS),
        dataset_name=f"test_{table_name}_dataset",
        dev_mode=True,
    )
    pipeline.run(resource, loader_file_format="parquet")
    return pipeline.dataset()


def test_build_entity_resource_pass_happy_path_all_valid(
    scenario_input_dir: Callable[[str], str],
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
    widget_model: type[BaseModel],
    tmp_path: Path,
) -> None:
    """Every valid record lands in the main table. No rejected table is created."""
    settings = settings_factory(input_dir=scenario_input_dir("happy_path"), output_dir=str(tmp_path))

    dataset = run_pipeline("widget", widget_model, settings)

    rows = dataset.widget.df()
    assert sorted(rows["widget_id"].tolist()) == ["a", "b", "c"]
    assert "widget_rejected" not in dataset.tables


def test_build_entity_resource_pass_mixed_batch_split_correctly(
    scenario_input_dir: Callable[[str], str],
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
    widget_model: type[BaseModel],
    tmp_path: Path,
) -> None:
    """A mixed batch puts valid rows in the main table and invalid rows in <table>_rejected."""
    settings = settings_factory(input_dir=scenario_input_dir("mixed"), output_dir=str(tmp_path))

    dataset = run_pipeline("widget", widget_model, settings)

    valid_rows = dataset.widget.df()
    rejected_rows = dataset.widget_rejected.df()
    assert sorted(valid_rows["widget_id"].tolist()) == ["a", "c"]
    assert len(rejected_rows) == 2


def test_build_entity_resource_fail_all_invalid_main_table_absent_or_empty(
    scenario_input_dir: Callable[[str], str],
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
    widget_model: type[BaseModel],
    tmp_path: Path,
) -> None:
    """When every record is invalid, the main table has zero rows or does not exist. Rejected has all rows."""
    settings = settings_factory(input_dir=scenario_input_dir("all_invalid"), output_dir=str(tmp_path))

    dataset = run_pipeline("widget", widget_model, settings)

    rejected_rows = dataset.widget_rejected.df()
    assert len(rejected_rows) == 3
    assert "widget" not in dataset.tables or len(dataset.widget.df()) == 0


def test_build_entity_resource_pass_empty_input_directory_yields_no_tables(
    tmp_path: Path, settings_factory: Callable[..., JsonlPydanticIngestSettings], widget_model: type[BaseModel]
) -> None:
    """An entity subdirectory with no matching files produces no tables."""
    input_dir = tmp_path / "empty_input"
    (input_dir / "widget").mkdir(parents=True)
    settings = settings_factory(input_dir=str(input_dir), output_dir=str(tmp_path))

    dataset = run_pipeline("widget", widget_model, settings)

    assert "widget" not in dataset.tables
    assert "widget_rejected" not in dataset.tables


def test_build_entity_resource_pass_missing_entity_directory_yields_no_tables(
    tmp_path: Path, settings_factory: Callable[..., JsonlPydanticIngestSettings], widget_model: type[BaseModel]
) -> None:
    """A missing entity subdirectory produces no tables and does not raise."""
    input_dir = tmp_path / "input_without_widget_dir"
    input_dir.mkdir()
    settings = settings_factory(input_dir=str(input_dir), output_dir=str(tmp_path))

    dataset = run_pipeline("widget", widget_model, settings)

    assert "widget" not in dataset.tables
    assert "widget_rejected" not in dataset.tables
