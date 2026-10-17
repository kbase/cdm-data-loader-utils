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

import dlt
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

VALID_SCHEMA_URI = "https://json-schema.org/draft/2020-12/schema"

SIMPLE_SCHEMA: dict[str, Any] = {
    "$schema": VALID_SCHEMA_URI,
    "required": ["number", "name"],
    "properties": {"number": {"type": "integer"}, "name": {"type": "string"}},
    "x-xsv-config": {"x-delimiter": ",", "x-has-header": True},
}


"""_process_xsv_file"""


import pendulum
from dlt.common.storages.fsspec_filesystem import FileItemDict


def make_file_item(file_path: Path) -> FileItemDict:
    """Build a dlt FileItemDict for a local file."""
    return FileItemDict(
        {
            "file_url": file_path.as_uri(),
            "file_name": file_path.name,
            "relative_path": file_path.name,
            "mime_type": "text/csv",
            "modification_date": pendulum.now(),
            "size_in_bytes": file_path.stat().st_size,
        }
    )


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
    output_dir = tmp_path / "duckdb_output"
    output_dir.mkdir()

    def fake_run_pipeline(
        *,
        resource: Any,  # noqa: ANN401
        pipeline_kwargs: dict[str, Any],
        **_: Any,  # noqa: ANN401
    ) -> None:
        assert pipeline_kwargs == {"pipeline_name": PIPELINE_NAME, "dataset_name": "xsv_test_dataset"}
        pipeline = dlt.pipeline(
            pipeline_name=pipeline_name,
            destination=dlt.destinations.duckdb(f"duckdb:///{output_dir!s}/{pipeline_name}.db"),
            dataset_name="xsv_test_dataset",
            pipelines_dir=str(tmp_path / "pipelines"),
        )
        captured["load_info"] = pipeline.run(resource)

    with patch.object(xsv_pipeline_module, "run_pipeline", fake_run_pipeline):
        run_xsv_pipeline(settings)

    load_info = captured["load_info"]
    assert not load_info.has_failed_jobs

    pipeline = dlt.pipeline(
        pipeline_name=pipeline_name,
        destination=dlt.destinations.duckdb(f"duckdb:///{output_dir!s}/{pipeline_name}.db"),
        dataset_name="xsv_test_dataset",
        pipelines_dir=str(tmp_path / "pipelines"),
    )
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT number, name FROM my_table ORDER BY number") as cur,
    ):
        rows = cur.fetchall()
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
        pipeline = dlt.pipeline(
            pipeline_name=pipeline_name,
            destination=dlt.destinations.duckdb(f"duckdb:///{settings.output_dir}/{pipeline_name}.db"),
            dataset_name="xsv_test_dataset",
            pipelines_dir=str(tmp_path / "pipelines"),
        )
        captured["load_info"] = pipeline.run(resource)

    with patch.object(xsv_pipeline_module, "run_pipeline", fake_run_pipeline):
        run_xsv_pipeline(settings)

    load_info = captured["load_info"]
    assert not load_info.has_failed_jobs
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
    settings = settings_factory(SIMPLE_SCHEMA, {"data.csv": "number,name\n1,alice\n2,bob\n"})

    argv = [
        "xsv",
        "--input-dir",
        str(settings.input_dir),
        "--output-dir",
        str(settings.output_dir),
        "--use-destination",
        "local_fs",
        "--dataset-name",
        "cli_dataset",
        "--schema-file",
        str(settings.schema_file),
        "--table-name",
        "my_table",
        "--use-output-dir-for-pipeline-metadata",
        "false",
        "--dlt-dev-mode",
        "false",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    pipeline_name = f"test_xsv_cli_pipeline_{uuid4().hex}"
    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        *,
        resource: Any,  # noqa: ANN401
        pipeline_kwargs: dict[str, Any],
        **_: Any,  # noqa: ANN401
    ) -> None:
        assert pipeline_kwargs == {"pipeline_name": PIPELINE_NAME, "dataset_name": "cli_dataset"}
        pipeline = dlt.pipeline(
            pipeline_name=pipeline_name,
            destination=dlt.destinations.duckdb(f"duckdb:///{settings.output_dir}/{pipeline_name}.db"),
            dataset_name="cli_dataset",
            pipelines_dir=str(tmp_path / "pipelines"),
        )
        captured["load_info"] = pipeline.run(resource)

    with patch.object(xsv_pipeline_module, "run_pipeline", fake_run_pipeline):
        cli.run()

    load_info = captured["load_info"]
    assert not load_info.has_failed_jobs

    pipeline = dlt.pipeline(
        pipeline_name=pipeline_name,
        destination=dlt.destinations.duckdb(f"duckdb:///{settings.output_dir}/{pipeline_name}.db"),
        dataset_name="cli_dataset",
        pipelines_dir=str(tmp_path / "pipelines"),
    )
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT COUNT(*) FROM my_table") as cur,
    ):
        (row_count,) = cur.fetchone()
    assert row_count == 2  # noqa: PLR2004
