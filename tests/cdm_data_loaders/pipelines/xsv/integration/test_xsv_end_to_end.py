"""End-to-end tests for the xsv pipeline, using a real qsv binary and local_fs destination."""

import sys
from collections.abc import Callable
from pathlib import Path
from uuid import uuid4

import pandas as pd
import pytest

from cdm_data_loaders.pipelines.xsv.pipeline import cli, run_xsv_pipeline
from cdm_data_loaders.pipelines.xsv.settings import PIPELINE_NAME, XsvIngestSettings
from tests.cdm_data_loaders.pipelines.conftest import DEFAULT_DLT_TABLES
from tests.cdm_data_loaders.pipelines.xsv.integration.conftest import (
    ALL_RAGGED_ROWS,
    COLUMNS,
    DEFAULT_XSV_SCHEMA,
    PARTIAL_RAGGED_ROWS,
    VALID_ROWS,
    rows_to_csv,
)

pytestmark = pytest.mark.usefixtures("qsv_cmd")

REJECTED_TABLE = "my_table_rejected"


def test_run_xsv_pipeline_pass_writes_expected_output(
    settings_factory: Callable[..., XsvIngestSettings],
    write_xsv_file: Callable[[XsvIngestSettings, str, str], Path],
) -> None:
    """A clean, valid file is loaded, with every row cast to its declared schema types."""
    settings = settings_factory()
    write_xsv_file(settings, "data.csv", rows_to_csv(VALID_ROWS))

    load_info = run_xsv_pipeline(settings)

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    dataset = load_info.pipeline.dataset()
    assert {"my_table", *DEFAULT_DLT_TABLES} == set(dataset.tables)

    records = dataset.table("my_table").df()[COLUMNS].sort_values("number").to_dict("records")
    assert records == [
        {"number": 2, "date": "2023-01-15", "float": "3.14", "boolean": True, "string": "key:value1"},
        {"number": 3, "date": "2023-02-20", "float": "1.11", "boolean": False, "string": "key:value2"},
    ]


def test_run_xsv_pipeline_pass_multiple_files_loaded_together(
    settings_factory: Callable[..., XsvIngestSettings],
    write_xsv_file: Callable[[XsvIngestSettings, str, str], Path],
) -> None:
    """Every file matching file_glob is cleaned, validated, and loaded into the same table."""
    settings = settings_factory()
    write_xsv_file(settings, "part_1.csv", rows_to_csv([VALID_ROWS[0]]))
    write_xsv_file(settings, "part_2.csv", rows_to_csv([VALID_ROWS[1]]))

    load_info = run_xsv_pipeline(settings)

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    df = load_info.pipeline.dataset().table("my_table").df()
    assert sorted(df["number"].tolist()) == [2, 3]


def test_run_xsv_pipeline_pass_partial_ragged_file_recovers_valid_rows(
    settings_factory: Callable[..., XsvIngestSettings],
    write_xsv_file: Callable[[XsvIngestSettings, str, str], Path],
) -> None:
    """A file with one ragged row loads its surviving valid rows AND records the row loss.

    qsv's first-pass validation splits ragged rows out with --split-ragged; the remaining valid
    subset is promoted through cleaning/validation and loaded, while the ragged-row count is also
    recorded as a rejected-table row for the file.
    """
    settings = settings_factory()
    write_xsv_file(settings, "data.csv", rows_to_csv(PARTIAL_RAGGED_ROWS))

    load_info = run_xsv_pipeline(settings)

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    dataset = load_info.pipeline.dataset()
    assert {"my_table", REJECTED_TABLE, *DEFAULT_DLT_TABLES} == set(dataset.tables)

    df = dataset.table("my_table").df()
    assert sorted(df["number"].tolist()) == [2, 4]

    rejected = dataset.table(REJECTED_TABLE).df()
    assert len(rejected) == 1
    assert rejected.iloc[0]["file"] == "data.csv"
    assert rejected.iloc[0]["message"] == "1 out of 3 records invalid."


def test_run_xsv_pipeline_pass_fully_ragged_file_is_entirely_rejected(
    settings_factory: Callable[..., XsvIngestSettings],
    write_xsv_file: Callable[[XsvIngestSettings, str, str], Path],
) -> None:
    """A file where every row is ragged yields no valid rows; the failure is recorded as rejected."""
    settings = settings_factory()
    write_xsv_file(settings, "bad.csv", rows_to_csv(ALL_RAGGED_ROWS))

    load_info = run_xsv_pipeline(settings)

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    dataset = load_info.pipeline.dataset()
    assert set(dataset.tables) == {REJECTED_TABLE, *DEFAULT_DLT_TABLES}

    rejected = dataset.table(REJECTED_TABLE).df()
    assert len(rejected) == 1
    assert rejected.iloc[0]["file"] == "bad.csv"
    assert rejected.iloc[0]["message"] == "2 out of 2 records invalid."


def test_run_xsv_pipeline_pass_no_matching_files_yields_no_data_table(
    settings_factory: Callable[..., XsvIngestSettings],
) -> None:
    """An input dir with no matching data files does not fail the run and creates no data table."""
    settings = settings_factory()

    load_info = run_xsv_pipeline(settings)

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    assert set(load_info.pipeline.dataset().tables) == DEFAULT_DLT_TABLES


def test_run_xsv_pipeline_pass_configured_null_regex_replaces_placeholder(
    settings_factory: Callable[..., XsvIngestSettings],
    write_xsv_file: Callable[[XsvIngestSettings, str, str], Path],
) -> None:
    """A schema-configured null_regex placeholder is replaced with a null value before loading."""
    schema = {
        **DEFAULT_XSV_SCHEMA,
        "x-xsv-config": {**DEFAULT_XSV_SCHEMA["x-xsv-config"], "x-null-regex": "^NA$"},
    }
    settings = settings_factory(schema=schema)
    write_xsv_file(
        settings,
        "data.csv",
        rows_to_csv(
            [["2", "2023-01-15", "NA", "true", "key:value1"], ["3", "2023-02-20", "1.11", "false", "key:value2"]]
        ),
    )

    load_info = run_xsv_pipeline(settings)

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    df = load_info.pipeline.dataset().table("my_table").df().sort_values("number")
    floats = df["float"].tolist()
    assert pd.isna(floats[0])
    assert floats[1] == "1.11"


def test_run_xsv_pipeline_pass_parquet_format_loads_expected_row_count(
    settings_factory: Callable[..., XsvIngestSettings],
    write_xsv_file: Callable[[XsvIngestSettings, str, str], Path],
) -> None:
    """The pipeline also loads successfully when loader_file_format is parquet.

    Only row count is checked here (not exact column values/types): dlt infers column types
    during parquet normalization (e.g. ISO-looking date strings may become typed date columns),
    which can differ from the jsonl path used by the other tests in this module.
    """
    settings = settings_factory(loader_file_format="parquet")
    write_xsv_file(settings, "data.csv", rows_to_csv(VALID_ROWS))

    load_info = run_xsv_pipeline(settings)

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    df = load_info.pipeline.dataset().table("my_table").df()
    assert len(df) == len(VALID_ROWS)


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    settings_factory: Callable[..., XsvIngestSettings],
    write_xsv_file: Callable[[XsvIngestSettings, str, str], Path],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """cli() reads command-line arguments and runs the pipeline, producing the expected output."""
    settings = settings_factory(dataset_name=f"cli_dataset_{uuid4().hex}")
    write_xsv_file(settings, "data.csv", rows_to_csv(VALID_ROWS))

    argv: list[str] = [
        PIPELINE_NAME,
        "--input-dir",
        settings.input_dir,
        "--output-dir",
        settings.output_dir,
        "--use-destination",
        settings.use_destination,
        "--dataset-name",
        settings.dataset_name,
        "--schema-file",
        str(settings.schema_file),
        "--table-name",
        "my_table",
        "--use-output-dir-for-pipeline-metadata",
        "false",
        "--dev-mode",
        "false",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    load_info = cli.run()

    assert load_info is not None
    assert load_info.has_failed_jobs is False
    df = load_info.pipeline.dataset().table("my_table").df()
    assert sorted(df["number"].tolist()) == [2, 3]
