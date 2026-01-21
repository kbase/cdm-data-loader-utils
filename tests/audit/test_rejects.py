"""Tests for the rejects audit table functions."""

import datetime
import json
import logging
from typing import Any

import pyspark.sql.functions as sf
import pytest
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import StringType, StructField, StructType

from cdm_data_loader_utils.audit.rejects import write_rejects
from cdm_data_loader_utils.audit.schema import (
    PARSED_ROW,
    PIPELINE,
    RAW_ROW,
    REJECTS,
    ROW_ERRORS,
    RUN_ID,
    SOURCE,
    TIMESTAMP,
)
from cdm_data_loader_utils.core.constants import INVALID_DATA_FIELD_NAME
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from tests.audit.conftest import create_table

n_rows_per_file = 4


@pytest.mark.requires_spark
def test_write_rejects(  # noqa: PLR0913
    spark: SparkSession,
    pipeline_run: PipelineRun,
    caplog: pytest.LogCaptureFixture,
    csv_schema: list[StructField],
    annotated_df_data: list[dict[str, Any]],
    annotated_df_schema: StructType,
) -> None:
    """Write out the rejects in the annotated_df dataframe."""
    create_table(spark, REJECTS)
    annotated_df = spark.createDataFrame(annotated_df_data, schema=annotated_df_schema)
    invalid_df = annotated_df.filter(sf.size(ROW_ERRORS) > 0)
    assert invalid_df.count() == n_rows_per_file * 4
    # filter out rows
    write_rejects(pipeline_run, annotated_df, csv_schema, INVALID_DATA_FIELD_NAME)

    rejects_df = spark.table(f"{pipeline_run.namespace}.{REJECTS}")
    rejects_rows = [r.asDict() for r in rejects_df.collect()]
    assert len(rejects_rows) == invalid_df.count()  # 4 lots of 4 invalid entries
    timestamp = rejects_rows[0][TIMESTAMP]
    for r in rejects_rows:
        # all should have the same timestamp
        assert isinstance(r[TIMESTAMP], datetime.datetime)
        assert r[TIMESTAMP] == timestamp
        assert r[RUN_ID] == pipeline_run.run_id
        assert r[PIPELINE] == pipeline_run.pipeline
        assert r[SOURCE] == pipeline_run.source_path
        # shit test
        assert r[ROW_ERRORS] in (
            ["parse_error"],
            ["missing_required: col2", "missing_required: col3", "missing_required: col4"],
            ["missing_required: col1", "missing_required: col2"],
            ["missing_required: col4", "missing_required: col5"],
            [
                "missing_required: col1",
                "missing_required: col2",
                "missing_required: col3",
                "missing_required: col4",
                "missing_required: col5",
            ],
        )
        # even worse tests
        assert isinstance(json.loads(r[PARSED_ROW]), dict)
        if r[ROW_ERRORS] != ["parse_error"]:
            assert r[RAW_ROW] is None
        else:
            assert len(r[RAW_ROW]) >= 1

    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.INFO
    assert (
        caplog.records[-1].message
        == f"{pipeline_run.pipeline} {pipeline_run.run_id}: invalid rows written to '{REJECTS}' audit table."
    )


@pytest.mark.requires_spark
def test_write_rejects_no_rejects(  # noqa: PLR0913
    spark: SparkSession,
    pipeline_run: PipelineRun,
    caplog: pytest.LogCaptureFixture,
    csv_schema: list[StructField],
    annotated_df_data: list[dict[str, Any]],
    annotated_df_schema: StructType,
) -> None:
    """Submit the annotated df to the write_rejects function for reject writing... but there are no rejects!"""
    create_table(spark, REJECTS)
    # filter out all invalid rows
    annotated_df = spark.createDataFrame(annotated_df_data, schema=annotated_df_schema)
    valid_df = annotated_df.filter(sf.size(ROW_ERRORS) == 0)
    assert valid_df.count() == n_rows_per_file
    write_rejects(
        pipeline_run,
        valid_df,
        csv_schema,
        INVALID_DATA_FIELD_NAME,
    )

    df = spark.table(f"{pipeline_run.namespace}.{REJECTS}")
    assert df.count() == 0

    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.INFO
    assert (
        caplog.records[-1].message
        == f"{pipeline_run.pipeline} {pipeline_run.run_id}: nothing to write to '{REJECTS}' audit table."
    )


@pytest.mark.requires_spark
def test_write_rejects_no_row_errors(
    empty_df: DataFrame,
    empty_df_schema: list[StructField],
    pipeline_run: PipelineRun,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Ensure that the writer behaves sensibly with an empty table as input."""
    err_msg = f"{pipeline_run.pipeline} {pipeline_run.run_id}: '{ROW_ERRORS}' column not present in dataframe; cannot record rejects."
    with pytest.raises(
        RuntimeError,
        match=err_msg,
    ):
        write_rejects(pipeline_run, empty_df, empty_df_schema, INVALID_DATA_FIELD_NAME)

    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.ERROR
    assert caplog.records[-1].message == err_msg


@pytest.mark.requires_spark
def test_write_rejects_empty_df(
    spark: SparkSession,
    empty_df_schema: list[StructField],
    pipeline_run: PipelineRun,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Ensure that the writer behaves sensibly with an empty table as input."""
    empty_df = spark.createDataFrame([], schema=StructType([*empty_df_schema, StructField(ROW_ERRORS, StringType())]))
    write_rejects(pipeline_run, empty_df, empty_df_schema, INVALID_DATA_FIELD_NAME)
    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.INFO
    assert (
        caplog.records[-1].message
        == f"{pipeline_run.pipeline} {pipeline_run.run_id}: nothing to write to '{REJECTS}' audit table."
    )
