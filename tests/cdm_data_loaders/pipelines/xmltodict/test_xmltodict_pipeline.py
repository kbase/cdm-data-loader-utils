"""Unit tests for run_xmltodict_pipeline and cli wiring."""

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any
from unittest.mock import patch
from uuid import uuid4

import dlt
import pytest
from dlt.common.pipeline import LoadInfo
from dlt.extract import DltResource

import cdm_data_loaders.pipelines.xmltodict.pipeline as xmltodict_module
from cdm_data_loaders.core.fields import LoaderFileFormatEnum
from cdm_data_loaders.pipelines.xmltodict.pipeline import (
    cli,
    run_xmltodict_pipeline,
)
from cdm_data_loaders.pipelines.xmltodict.settings import PIPELINE_NAME, XmlToDictSettings
from tests.cdm_data_loaders.pipelines.conftest import duckdb_pipeline
from tests.cdm_data_loaders.pipelines.helpers import assert_run_ok, fetch_all_duckdb, settings_to_argv
from tests.xml_samples import TWO_BOOK_LIBRARY_XML


@pytest.mark.parametrize("loader_file_format", list(LoaderFileFormatEnum), ids=str)
def test_run_xmltodict_pipeline_pass_sets_core_run_pipeline_args_correctly(
    settings_factory: Callable[..., XmlToDictSettings],
    loader_file_format: str,
) -> None:
    """run_xmltodict_pipeline builds the reader and delegates to run_pipeline with the correct args."""
    settings = settings_factory(loader_file_format=loader_file_format)

    with patch.object(xmltodict_module, "run_pipeline") as mock_run_pipeline:
        run_xmltodict_pipeline(settings)

    assert mock_run_pipeline.call_count == 1
    _, kwargs = mock_run_pipeline.call_args
    assert kwargs.keys() == {"settings", "resource", "pipeline_kwargs", "pipeline_run_kwargs"}
    assert kwargs["settings"] == settings
    assert kwargs["pipeline_kwargs"] == {
        "pipeline_name": PIPELINE_NAME,
        "dataset_name": settings.dataset_name,
    }
    assert kwargs["pipeline_run_kwargs"] == {"loader_file_format": loader_file_format}
    assert isinstance(kwargs["resource"], DltResource)


def test_run_xmltodict_pipeline_pass_binds_reader_before_run_pipeline(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """The resource passed to run_pipeline is the bound reader piped into the filesystem source."""
    settings = settings_factory()

    with patch.object(xmltodict_module, "run_pipeline") as mock_run_pipeline:
        run_xmltodict_pipeline(settings)

    mock_run_pipeline.assert_called_once()
    resource = mock_run_pipeline.call_args.kwargs["resource"]
    assert resource.name == "xmltodict_reader"
    assert resource._pipe.parent.name == "filesystem"  # noqa: SLF001
    assert isinstance(resource, dlt.sources.DltResource)


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    tmp_path: Path,
    settings_factory: Callable[..., XmlToDictSettings],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """cli() reads command-line arguments and runs the pipeline via run_xmltodict_pipeline.

    The real XmlToDictSettings parses argv (seeded by the autouse dlt
    config isolation fixture), and core.run_pipeline is redirected to a real
    DuckDB pipeline so the full flow is validated.
    """
    settings = settings_factory(dataset_name="cli_dataset", dlt_dev_mode=False)
    (Path(settings.input_dir) / "library.xml").write_text(TWO_BOOK_LIBRARY_XML, encoding="utf-8")
    monkeypatch.setattr(sys, "argv", [PIPELINE_NAME, *settings_to_argv(settings)])

    # unique pipeline name: the duckdb database file is derived from it, so
    # repeated runs of the suite do not see rows from previous runs
    pipeline_name = f"test_xml_cli_pipeline_{uuid4().hex}"

    def fake_run_pipeline(
        *,
        resource: DltResource,
        pipeline_kwargs: dict[str, Any],
        **_: dict[str, Any],
    ) -> LoadInfo | None:
        assert pipeline_kwargs == {"pipeline_name": PIPELINE_NAME, "dataset_name": "cli_dataset"}
        pipeline = duckdb_pipeline(tmp_path, pipeline_name, "cli_dataset")
        return pipeline.run(resource)

    with patch.object(xmltodict_module, "run_pipeline", fake_run_pipeline):
        load_info = assert_run_ok(cli.run())

    assert fetch_all_duckdb(load_info.pipeline, "SELECT COUNT(*) FROM book") == [(2,)]


def test_cli_pass_writes_rows_for_matching_xml_files(
    tmp_path: Path,
    settings_factory: Callable[..., XmlToDictSettings],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """cli() with xml files in the input dir loads one row per matching element."""
    settings = settings_factory(dataset_name="cli_dataset", dlt_dev_mode=False)
    (Path(settings.input_dir) / "data.xml").write_text(TWO_BOOK_LIBRARY_XML, encoding="utf-8")
    monkeypatch.setattr(sys, "argv", [PIPELINE_NAME, *settings_to_argv(settings)])

    pipeline_name = f"test_xml_cli_rows_pipeline_{uuid4().hex}"

    def fake_run_pipeline(
        *,
        resource: Any,  # noqa: ANN401
        **_: Any,  # noqa: ANN401
    ) -> LoadInfo | None:
        pipeline = duckdb_pipeline(tmp_path, pipeline_name, "cli_dataset")
        return pipeline.run(resource)

    with patch.object(xmltodict_module, "run_pipeline", fake_run_pipeline):
        load_info = assert_run_ok(cli.run())

    assert fetch_all_duckdb(load_info.pipeline, "SELECT COUNT(*) FROM book") == [(2,)]


@pytest.mark.parametrize("missing_field", ["dataset_name", "table_name", "xml_tag"])
def test_cli_fail_missing_required_args_raise(
    settings_factory: Callable[..., XmlToDictSettings],
    monkeypatch: pytest.MonkeyPatch,
    missing_field: str,
) -> None:
    """cli() without required pipeline-specific args raises before running the pipeline."""
    settings = settings_factory()
    monkeypatch.setattr(sys, "argv", [PIPELINE_NAME, *settings_to_argv(settings, exclude=[missing_field])])

    with (
        patch.object(xmltodict_module, "run_pipeline") as mock_run_pipeline,
        pytest.raises(Exception, match=missing_field),
    ):
        cli()

    mock_run_pipeline.assert_not_called()
