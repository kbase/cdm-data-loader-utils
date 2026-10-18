"""End-to-end tests for the jsonlines pipeline."""

import sys
from collections.abc import Callable
from pathlib import Path

import pytest

from cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline import (
    cli,
    load_entity_models,
    run_jsonlines_with_validation_pipeline,
)
from cdm_data_loaders.pipelines.jsonlines.settings import PYDANTIC_PIPELINE_NAME, JsonlPydanticIngestSettings


def test_run_jsonlines_pipeline_pass_writes_expected_parquet(
    scenario_input_dir: Callable[[str], str],
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
) -> None:
    """Running the pipeline on a mixed batch puts the correct rows in both output tables."""
    settings = settings_factory(input_dir=scenario_input_dir("mixed"))

    load_info = run_jsonlines_with_validation_pipeline(settings)
    dataset = load_info.pipeline.dataset()
    assert sorted(dataset.widget.df()["widget_id"].tolist()) == ["a", "c"]
    assert len(dataset.widget_rejected.df()) == 2

    # get stats on tables
    # dataset.row_counts(table_names=dataset.tables).fetchall()


def test_run_jsonlines_pipeline_fail_unknown_table_name_raises(
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
) -> None:
    """A table name in table_names that is not in the entity model registry raises ValueError."""
    settings = settings_factory(table_names=["not_a_real_table"])

    with pytest.raises(ValueError, match="Unknown table name"):
        run_jsonlines_with_validation_pipeline(settings)


def test_run_jsonlines_pipeline_pass_table_names_none_processes_every_registered_table(
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
) -> None:
    """table_names=None resolves to every entity in entity_models_module."""
    settings = settings_factory()

    entity_models = load_entity_models(settings.entity_models_module)
    resolved = settings.table_names or sorted(entity_models)

    assert resolved == ["widget"]


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    scenario_input_dir: Callable[[str], str],
    dlt_destination_config: str,
    entity_models_module: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """cli() reads command-line arguments and runs the pipeline, producing the expected output."""
    log_config_file = tmp_path / "logging.json"
    log_config_file.write_text('{"version": 1}')
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    argv = [
        PYDANTIC_PIPELINE_NAME,
        "--input-dir",
        scenario_input_dir("happy_path"),
        "--output-dir",
        str(output_dir),
        "--log-config-file",
        str(log_config_file),
        "--use-destination",
        dlt_destination_config,
        "--entity-models-module",
        entity_models_module,
        "--dataset-name",
        "some_cool_dataset",
        "--use-output-dir-for-pipeline-metadata",
        "true",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    load_info = cli.run()

    assert load_info is not None
    assert not load_info.has_failed_jobs
    assert load_info.pipeline.pipeline_name == PYDANTIC_PIPELINE_NAME
    dataset = load_info.pipeline.dataset()
    assert sorted(dataset.widget.df()["widget_id"].tolist()) == ["a", "b", "c"]
