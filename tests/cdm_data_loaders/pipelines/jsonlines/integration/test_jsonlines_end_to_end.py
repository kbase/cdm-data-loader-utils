"""End-to-end tests for the jsonlines pipeline."""

from collections.abc import Callable
from pathlib import Path

import pytest

from cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline import (
    load_entity_models,
)
from cdm_data_loaders.pipelines.jsonlines.settings import PYDANTIC_PIPELINE_NAME, JsonlPydanticIngestSettings
from tests.cdm_data_loaders.pipelines.framework import PipelineTestRunner
from tests.cdm_data_loaders.pipelines.helpers import (
    assert_pipeline_report_written,
    assert_run_ok,
)
from tests.helpers import settings_to_argv


def test_run_jsonlines_pipeline_pass_writes_expected_parquet(
    scenario_input_dir: Callable[[str], str],
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Running the pipeline on a mixed batch puts the correct rows in both output tables."""
    settings = settings_factory(input_dir=scenario_input_dir("mixed"))

    load_info = assert_run_ok(pipeline_runner.run(PYDANTIC_PIPELINE_NAME, settings=settings))
    dataset = load_info.pipeline.dataset()
    assert sorted(dataset.widget.df()["widget_id"].tolist()) == ["a", "c"]
    assert len(dataset.widget_rejected.df()) == 2
    assert_pipeline_report_written(Path(settings.output_dir), load_info)

    # FIXME: remove or improve
    # get stats on tables
    # dataset.row_counts(table_names=dataset.tables).fetchall()


def test_run_jsonlines_pipeline_fail_unknown_table_name_raises(
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """A table name in table_names that is not in the entity model registry raises ValueError."""
    settings = settings_factory(table_names=["not_a_real_table"])

    with pytest.raises(ValueError, match="Unknown table name"):
        pipeline_runner.run(PYDANTIC_PIPELINE_NAME, settings=settings)


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
    settings_factory: Callable[..., JsonlPydanticIngestSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """cli() reads command-line arguments and runs the pipeline, producing the expected output."""
    settings = settings_factory(
        input_dir=scenario_input_dir("happy_path"),
        dataset_name="some_cool_dataset",
        dlt_dev_mode=False,
        use_output_dir_for_pipeline_metadata=True,
    )

    load_info = assert_run_ok(pipeline_runner.run(PYDANTIC_PIPELINE_NAME, cli_args=settings_to_argv(settings)))
    assert load_info.pipeline.pipeline_name == PYDANTIC_PIPELINE_NAME
    dataset = load_info.pipeline.dataset()
    assert sorted(dataset.widget.df()["widget_id"].tolist()) == ["a", "b", "c"]
    assert_pipeline_report_written(Path(settings.output_dir), load_info, metadata_moved=True)
