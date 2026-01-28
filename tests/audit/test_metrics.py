"""Tests for the metrics audit table functions."""

import datetime
import logging
from typing import Any

import pytest
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import StructType

from cdm_data_loader_utils.audit.metrics import write_metrics
from cdm_data_loader_utils.audit.schema import (
    METRICS,
    N_INVALID,
    N_READ,
    N_VALID,
    PIPELINE,
    RUN_ID,
    SOURCE,
    VALIDATION_ERRORS,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from tests.audit.conftest import (
    DEFAULT_DATA,
    INIT_TIMESTAMP_FIELDS,
    check_saved_data,
    create_table,
)
from tests.conftest import NAMESPACE, PIPELINE_RUN, TEST_NS


@pytest.mark.requires_spark
@pytest.mark.parametrize("pipeline", [PIPELINE_RUN[PIPELINE], "pipe_2"])
@pytest.mark.parametrize("source", [PIPELINE_RUN[SOURCE], "/src/two"])
@pytest.mark.parametrize("add_default_data", [True, False])
def test_write_metrics(  # noqa: PLR0913
    spark: SparkSession,
    pipeline: str,
    source: str,
    add_default_data: bool,
    caplog: pytest.LogCaptureFixture,
    annotated_df_data: list[dict[str, Any]],
    annotated_df_schema: StructType,
    annotated_df_errors: set[str],
) -> None:
    """Calculate the metrics for a given set of results and add it to the metrics table."""
    create_table(spark, METRICS, add_default_data=add_default_data)
    is_saved_table = add_default_data and pipeline == PIPELINE_RUN[PIPELINE] and source == PIPELINE_RUN[SOURCE]
    pipeline_args = {RUN_ID: PIPELINE_RUN[RUN_ID], PIPELINE: pipeline, SOURCE: source}
    pipeline_run = PipelineRun(**{**pipeline_args, NAMESPACE: TEST_NS})

    annotated_df = spark.createDataFrame(annotated_df_data, schema=annotated_df_schema)
    metrics_out = {N_READ: 20, N_VALID: 4, N_INVALID: 16, VALIDATION_ERRORS: sorted(annotated_df_errors)}

    metrics = write_metrics(spark, annotated_df, pipeline_run)
    assert [metrics.asDict()] == [metrics_out]

    df = spark.table(f"{pipeline_run.namespace}.{METRICS}")
    if is_saved_table or not add_default_data:
        # if we've either updated an existing row or added a row afresh, expect one data point
        assert len(df.collect()) == 1
        check_saved_data(spark, METRICS, {**pipeline_args, **metrics_out})
    elif add_default_data:
        # expect two rows
        assert len(df.collect()) == 2
        row_data = [r.asDict() for r in df.collect()]
        for row in row_data:
            # filter out timestamp rows
            for f in INIT_TIMESTAMP_FIELDS.get(METRICS, []):
                assert isinstance(row[f], datetime.datetime)
                del row[f]
        assert DEFAULT_DATA[METRICS] in row_data
        assert {**pipeline_args, **metrics_out} in row_data

    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.INFO
    assert (
        caplog.records[-1].message
        == f"{pipeline_run.pipeline} {pipeline_run.run_id}: ingest metrics written to '{METRICS}' table."
    )


@pytest.mark.requires_spark
def test_write_metrics_empty_df(
    spark: SparkSession, pipeline_run: PipelineRun, caplog: pytest.LogCaptureFixture, empty_df: DataFrame
) -> None:
    """Ensure that the writer behaves sensibly with an empty table as input."""
    metrics = write_metrics(spark, empty_df, pipeline_run)
    assert metrics.asDict() == {N_READ: 0, N_INVALID: 0, N_VALID: 0, VALIDATION_ERRORS: []}
    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.INFO
    assert (
        caplog.records[-1].message
        == f"{pipeline_run.pipeline} {pipeline_run.run_id}: nothing to write to '{METRICS}' audit table."
    )
