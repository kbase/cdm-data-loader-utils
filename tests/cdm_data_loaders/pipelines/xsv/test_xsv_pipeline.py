"""Tests for _process_xsv_file, run_xsv_pipeline, and cli wiring.

These tests drive the real qsv binary (see the `qsv_cmd` fixture in conftest.py); tests are
skipped when qsv is not available on PATH or via the QSV_BIN env var.
"""

import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any
from unittest.mock import patch
from uuid import uuid4

import pytest

import cdm_data_loaders.pipelines.xsv.pipeline as xsv_pipeline_module
from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.pipelines.xsv.pipeline import (
    REJECTED_SUFFIX,
    XsvWorkPaths,
    _process_xsv_file,
    cli,
    run_xsv_pipeline,
)
from cdm_data_loaders.pipelines.xsv.settings import PIPELINE_NAME, XsvIngestSettings
from tests.cdm_data_loaders.pipelines.conftest import duckdb_pipeline
from tests.cdm_data_loaders.pipelines.helpers import (
    assert_run_ok,
    fetch_all_duckdb,
    make_file_item,
    settings_to_argv,
)
from tests.cdm_data_loaders.pipelines.xsv.conftest import VALID_SCHEMA_URI

SIMPLE_SCHEMA: dict[str, Any] = {
    "$schema": VALID_SCHEMA_URI,
    "required": ["number", "name"],
    "properties": {"number": {"type": "integer"}, "name": {"type": "string"}},
    "x-xsv-config": {"x-delimiter": ",", "x-has-header": True},
}


"""_process_xsv_file"""


def test_process_xsv_file_pass_routes_valid_rows_to_table(
    settings_factory: Callable[..., XsvIngestSettings],
    work_paths_factory: Callable[[XsvIngestSettings], XsvWorkPaths],
) -> None:
    """A clean, valid file is loaded, with every row cast to its declared schema types."""
    settings = settings_factory(SIMPLE_SCHEMA, {"data.csv": "number,name\n1,alice\n2,bob\n"})
    work_paths = work_paths_factory(settings)

    results = list(_process_xsv_file(make_file_item(Path(settings.input_dir) / "data.csv"), settings, work_paths))

    assert [(result.meta.table_name, result.data) for result in results] == [
        ("my_table", [{"number": 1, "name": "alice"}, {"number": 2, "name": "bob"}]),
    ]


def test_process_xsv_file_pass_replaces_configured_null_placeholder(
    settings_factory: Callable[..., XsvIngestSettings],
    work_paths_factory: Callable[[XsvIngestSettings], XsvWorkPaths],
) -> None:
    """A configured null_regex placeholder is replaced with None in the loaded rows."""
    schema = {**SIMPLE_SCHEMA, "properties": {"number": {"type": "integer"}, "name": {"type": ["string", "null"]}}}
    schema["x-xsv-config"] = {**schema["x-xsv-config"], "x-null-regex": "^NA$"}
    settings = settings_factory(schema, {"data.csv": "number,name\n1,NA\n2,bob\n"})
    work_paths = work_paths_factory(settings)

    results = list(_process_xsv_file(make_file_item(Path(settings.input_dir) / "data.csv"), settings, work_paths))

    assert results[0].data == [{"number": 1, "name": None}, {"number": 2, "name": "bob"}]


def test_process_xsv_file_pass_routes_total_failure_to_rejected_table(
    settings_factory: Callable[..., XsvIngestSettings],
    work_paths_factory: Callable[[XsvIngestSettings], XsvWorkPaths],
) -> None:
    """A file where every row is ragged (wrong column count) yields no valid rows and one rejected row."""
    settings = settings_factory(SIMPLE_SCHEMA, {"bad.csv": "number,name\n1\n2,bob,extra\n"})
    work_paths = work_paths_factory(settings)

    results = list(_process_xsv_file(make_file_item(Path(settings.input_dir) / "bad.csv"), settings, work_paths))

    table_names = [result.meta.table_name for result in results]
    assert table_names == [f"my_table{REJECTED_SUFFIX}"]
    (rejected,) = results[0].data
    assert rejected["file"] == "bad.csv"
    assert rejected["message"]


"""run_xsv_pipeline"""


@pytest.mark.usefixtures("qsv_cmd")
def test_run_xsv_pipeline_pass_end_to_end_via_duckdb(
    settings_factory: Callable[..., XsvIngestSettings],
    tmp_path: Path,
) -> None:
    """The pipeline discovers data files (excluding the schema file), skipping non-matching files.

    Runs against a real DuckDB destination via a patched run_pipeline, mirroring the xml2db
    pipeline's end-to-end test pattern.
    """
    settings = settings_factory(
        SIMPLE_SCHEMA,
        {
            "part_1.csv": "number,name\n1,alice\n",
            "part_2.csv": "number,name\n2,bob\n",
        },
    )

    pipeline_name = f"test_xsv_pipeline_{uuid4().hex}"
    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        *,
        resource: Any,  # noqa: ANN401
        pipeline_kwargs: dict[str, Any],
        **_: Any,  # noqa: ANN401
    ) -> None:
        assert pipeline_kwargs == {"pipeline_name": PIPELINE_NAME, "dataset_name": "xsv_test_dataset"}
        pipeline = duckdb_pipeline(tmp_path, pipeline_name, "xsv_test_dataset")
        captured["load_info"] = pipeline.run(resource)

    with patch.object(xsv_pipeline_module, "run_pipeline", fake_run_pipeline):
        run_xsv_pipeline(settings)

    assert_run_ok(captured["load_info"])

    pipeline = duckdb_pipeline(tmp_path, pipeline_name, "xsv_test_dataset")
    rows = fetch_all_duckdb(pipeline, "SELECT number, name FROM my_table ORDER BY number")
    assert [tuple(row) for row in rows] == [(1, "alice"), (2, "bob")]


@pytest.mark.usefixtures("qsv_cmd")
def test_run_xsv_pipeline_pass_no_matching_files_loads_nothing(
    settings_factory: Callable[..., XsvIngestSettings],
    tmp_path: Path,
) -> None:
    """When no data files are present besides the schema, no rows are extracted."""
    settings = settings_factory(SIMPLE_SCHEMA, {}, file_glob="*.csv")

    pipeline_name = f"test_xsv_pipeline_empty_{uuid4().hex}"
    captured: dict[str, Any] = {}

    def fake_run_pipeline(*, resource: Any, **_: Any) -> None:  # noqa: ANN401
        pipeline = duckdb_pipeline(tmp_path, pipeline_name, "xsv_test_dataset")
        captured["load_info"] = pipeline.run(resource)

    with patch.object(xsv_pipeline_module, "run_pipeline", fake_run_pipeline):
        run_xsv_pipeline(settings)

    load_info = assert_run_ok(captured["load_info"])
    assert "my_table" not in load_info.pipeline.default_schema.tables


"""cli"""


def test_cli_pass_calls_run_cli_with_settings_class_and_pipeline_fn() -> None:
    """Cli is a PipelineCli wired to XsvIngestSettings and run_xsv_pipeline."""
    assert isinstance(cli, PipelineCli)
    assert cli.settings_cls is XsvIngestSettings
    assert cli.pipeline_fn is run_xsv_pipeline


@pytest.mark.usefixtures("qsv_cmd")
def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    settings_factory: Callable[..., XsvIngestSettings],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """cli() reads command-line arguments and loads rows via a real DuckDB pipeline."""
    settings = settings_factory(
        SIMPLE_SCHEMA,
        {"data.csv": "number,name\n1,alice\n2,bob\n"},
        dataset_name="cli_dataset",
    )
    monkeypatch.setattr(sys, "argv", ["xsv", *settings_to_argv(settings)])

    pipeline_name = f"test_xsv_cli_pipeline_{uuid4().hex}"
    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        *,
        resource: Any,  # noqa: ANN401
        pipeline_kwargs: dict[str, Any],
        **_: Any,  # noqa: ANN401
    ) -> None:
        assert pipeline_kwargs == {"pipeline_name": PIPELINE_NAME, "dataset_name": "cli_dataset"}
        pipeline = duckdb_pipeline(tmp_path, pipeline_name, "cli_dataset")
        captured["load_info"] = pipeline.run(resource)

    with patch.object(xsv_pipeline_module, "run_pipeline", fake_run_pipeline):
        cli.run()

    assert_run_ok(captured["load_info"])

    pipeline = duckdb_pipeline(tmp_path, pipeline_name, "cli_dataset")
    ((row_count,),) = fetch_all_duckdb(pipeline, "SELECT COUNT(*) FROM my_table")
    assert row_count == 2  # noqa: PLR2004
