"""Unit tests for run_jsonlines_pipeline and cli wiring."""

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any
from unittest.mock import patch

import dlt
import pytest
from dlt.common.pipeline import LoadInfo
from dlt.extract import DltResource

import cdm_data_loaders.pipelines.jsonlines.extract_pipeline as extract_pipeline_module
from cdm_data_loaders.core.fields import LoaderFileFormatEnum
from cdm_data_loaders.pipelines.jsonlines.extract_pipeline import cli, run_jsonlines_pipeline
from cdm_data_loaders.pipelines.jsonlines.settings import EXTRACT_PIPELINE_NAME, JsonlExtractSettings
from tests.cdm_data_loaders.pipelines.conftest import duckdb_pipeline

SIMPLE_JSONL = """{"widget_id": "a", "count": 1}
{"widget_id": "b", "count": 2}
"""


@pytest.mark.parametrize("loader_file_format", list(LoaderFileFormatEnum), ids=str)
def test_run_jsonlines_pipeline_pass_sets_core_run_pipeline_args_correctly(
    extract_settings_factory: Callable[..., JsonlExtractSettings],
    loader_file_format: str,
) -> None:
    """run_jsonlines_pipeline builds the reader and delegates to run_pipeline with the correct args."""
    settings = extract_settings_factory(loader_file_format=loader_file_format)

    with patch.object(extract_pipeline_module, "run_pipeline") as mock_run_pipeline:
        run_jsonlines_pipeline(settings)

    assert mock_run_pipeline.call_count == 1
    _, kwargs = mock_run_pipeline.call_args
    assert kwargs.keys() == {"settings", "resource", "pipeline_kwargs", "pipeline_run_kwargs"}
    assert kwargs["settings"] == settings
    assert kwargs["pipeline_kwargs"] == {
        "pipeline_name": EXTRACT_PIPELINE_NAME,
        "dataset_name": settings.dataset_name,
    }
    assert kwargs["pipeline_run_kwargs"] == {"loader_file_format": loader_file_format}
    assert isinstance(kwargs["resource"], DltResource)


def test_run_jsonlines_pipeline_pass_resource_is_reader_piped_into_source(
    extract_settings_factory: Callable[..., JsonlExtractSettings],
) -> None:
    """The resource passed to run_pipeline combines the filesystem source with the bound reader."""
    settings = extract_settings_factory()

    with patch.object(extract_pipeline_module, "run_pipeline") as mock_run_pipeline:
        run_jsonlines_pipeline(settings)

    resource = mock_run_pipeline.call_args.kwargs["resource"]
    assert isinstance(resource, dlt.sources.DltResource)


def test_run_jsonlines_pipeline_pass_routes_lines_by_directory(
    extract_settings_factory: Callable[..., JsonlExtractSettings],
) -> None:
    """Top-level and nested files all land in the table_name table regardless of directory."""
    settings = extract_settings_factory(file_glob="**/*.jsonl*")
    input_dir = Path(settings.input_dir)
    (input_dir / "widget").mkdir(parents=True)
    (input_dir / "widget" / "data.jsonl").write_text(SIMPLE_JSONL, encoding="utf-8")
    (input_dir / "gadget").mkdir()
    (input_dir / "gadget" / "data.jsonl").write_text('{"widget_id": "g", "count": 0}\n', encoding="utf-8")

    load_info = run_jsonlines_pipeline(settings)

    assert load_info is not None
    assert not load_info.has_failed_jobs
    dataset = load_info.pipeline.dataset()
    assert sorted(dataset.widget.df()["widget_id"].tolist()) == ["a", "b", "g"]


def test_run_jsonlines_pipeline_pass_invalid_rows_routed_to_invalid_table(
    extract_settings_factory: Callable[..., JsonlExtractSettings],
) -> None:
    """Lines that are not valid JSON land in <table_name>_invalid with parse error detail."""
    settings = extract_settings_factory(file_glob="**/*.jsonl*")
    input_dir = Path(settings.input_dir)
    (input_dir / "widget").mkdir(parents=True)
    (input_dir / "widget" / "data.jsonl").write_text('{"widget_id": "a"}\n{broken\n', encoding="utf-8")

    load_info = run_jsonlines_pipeline(settings)

    assert load_info is not None
    assert not load_info.has_failed_jobs
    dataset = load_info.pipeline.dataset()
    assert sorted(dataset.widget.df()["widget_id"].tolist()) == ["a"]
    invalid_df = dataset.widget_invalid.df()
    assert len(invalid_df) == 1
    assert invalid_df.iloc[0]["raw_record"] == "{broken"


def test_run_jsonlines_pipeline_pass_no_matching_files_yields_no_tables(
    extract_settings_factory: Callable[..., JsonlExtractSettings],
) -> None:
    """An input dir with no matching files does not fail the run and creates no data tables."""
    settings = extract_settings_factory()
    Path(settings.input_dir).mkdir(exist_ok=True)

    load_info = run_jsonlines_pipeline(settings)

    assert load_info is not None
    assert not load_info.has_failed_jobs
    dataset = load_info.pipeline.dataset()
    assert set(dataset.tables) == {"_dlt_version", "_dlt_loads", "_dlt_pipeline_state", "_load_info"}


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    extract_settings_factory: Callable[..., JsonlExtractSettings],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """cli() reads command-line arguments and runs the pipeline, producing the expected output."""
    settings = extract_settings_factory()
    input_dir = Path(settings.input_dir)
    (input_dir / "data.jsonl").write_text(SIMPLE_JSONL, encoding="utf-8")
    output_dir = tmp_path / "cli_output"
    output_dir.mkdir()
    log_config_file = tmp_path / "logging.json"
    log_config_file.write_text('{"version": 1}')

    argv = [
        "jsonlines",
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--log-config-file",
        str(log_config_file),
        "--use-destination",
        str(settings.use_destination),
        "--dataset-name",
        "cli_dataset",
        "--table-name",
        "widget",
        "--use-output-dir-for-pipeline-metadata",
        "false",
        "--dev-mode",
        "false",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    def fake_run_pipeline(
        *, resource: DltResource, pipeline_kwargs: dict[str, Any], **_: dict[str, Any]
    ) -> LoadInfo | None:
        assert pipeline_kwargs == {"pipeline_name": EXTRACT_PIPELINE_NAME, "dataset_name": "cli_dataset"}
        pipeline = duckdb_pipeline(tmp_path, pipeline_kwargs["pipeline_name"], pipeline_kwargs["dataset_name"])
        return pipeline.run(resource)

    with patch.object(extract_pipeline_module, "run_pipeline", fake_run_pipeline):
        load_info = cli.run()

    assert load_info is not None
    assert not load_info.has_failed_jobs
    # n.b. no _load_info due to this being a monkeypatched pipeline
    assert set(load_info.pipeline.dataset().tables) == {"widget", "_dlt_version", "_dlt_loads", "_dlt_pipeline_state"}
    assert sorted(load_info.pipeline.dataset().widget.df()["widget_id"].tolist()) == ["a", "b"]
