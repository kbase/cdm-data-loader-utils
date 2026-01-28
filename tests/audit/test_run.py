"""Tests for the run audit table functions."""

import datetime
import logging

import pytest
from pyspark.sql import SparkSession

from cdm_data_loader_utils.audit.run import complete_run, fail_run, start_run
from cdm_data_loader_utils.audit.schema import (
    END_TIME,
    ERROR,
    PIPELINE,
    RECORDS_PROCESSED,
    RUN,
    RUN_ID,
    SOURCE,
    START_TIME,
    STATUS,
    STATUS_ERROR,
    STATUS_SUCCESS,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from tests.audit.conftest import DEFAULT_DATA, check_saved_data, create_table
from tests.conftest import PIPELINE_RUN, TEST_NS


def check_no_record_log_message(
    run: PipelineRun, table_name: str, status: str, caplog_message: logging.LogRecord
) -> None:
    """Checker for the log message if an update is attempted to an entry that doesn't exist.

    :param run: pipeline run
    :type run: PipelineRun
    :param table_name: table name
    :type table_name: str
    :param status: status
    :type status: str
    :param caplog_message: the log entry to test
    :type caplog_message:
    """
    assert caplog_message.levelno == logging.WARNING
    assert (
        caplog_message.message
        == f"{run.pipeline} {run.run_id}: cannot update '{table_name}' to status {status} because no record exists."
    )


# Tests for the run-related functions
@pytest.mark.requires_spark
def test_start_run(spark: SparkSession, pipeline_run: PipelineRun) -> None:
    """Test that starting a run creates a new entry in the RUN audit table."""
    create_table(spark, RUN)
    df = spark.table(f"{TEST_NS}.{RUN}")
    assert len(df.collect()) == 0

    start_run(spark, pipeline_run)
    check_saved_data(spark, RUN, DEFAULT_DATA[RUN])


@pytest.mark.requires_spark
@pytest.mark.parametrize("pipeline", [PIPELINE_RUN[PIPELINE], "pipe_2"])
@pytest.mark.parametrize("source", [PIPELINE_RUN[SOURCE], "/src/two"])
@pytest.mark.parametrize("add_default_data", [True, False])
def test_complete_run(
    spark: SparkSession,
    pipeline: str,
    source: str,
    add_default_data: bool,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test that completing a run updates the appropriate run table."""
    create_table(spark, RUN, add_default_data=add_default_data)
    is_saved_table = add_default_data and pipeline == PIPELINE_RUN[PIPELINE] and source == PIPELINE_RUN[SOURCE]
    pipeline_run = PipelineRun(PIPELINE_RUN[RUN_ID], pipeline, source, TEST_NS)
    n_records = 123_456_789_012_345

    complete_run(
        spark,
        pipeline_run,
        records_processed=n_records,
    )

    rows = [r.asDict() for r in spark.table(f"{TEST_NS}.{RUN}").collect()]
    if is_saved_table:
        assert len(rows) == 1
        row = rows[0]
        assert row[STATUS] == STATUS_SUCCESS
        assert row[RECORDS_PROCESSED] == n_records
        assert isinstance(row[END_TIME], datetime.datetime)
        assert row[END_TIME] > row[START_TIME]
        assert len(caplog.records) == 1

    else:
        if add_default_data:
            # only row will be the default data
            assert len(rows) == 1
            check_saved_data(spark, RUN, DEFAULT_DATA[RUN])
        else:
            # no data successfully added
            assert len(rows) == 0
        # two log messages: invalid row, completion
        assert len(caplog.records) == 2
        check_no_record_log_message(pipeline_run, RUN, STATUS_SUCCESS, caplog.records[0])

    assert caplog.records[-1].levelno == logging.INFO
    assert caplog.records[-1].message == f"{pipeline} {PIPELINE_RUN[RUN_ID]}: run completed"


@pytest.mark.requires_spark
@pytest.mark.parametrize("pipeline", [PIPELINE_RUN[PIPELINE], "pipe_2"])
@pytest.mark.parametrize("source", [PIPELINE_RUN[SOURCE], "/src/two"])
@pytest.mark.parametrize("add_default_data", [True, False])
def test_fail_run(
    spark: SparkSession,
    pipeline: str,
    source: str,
    add_default_data: bool,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test that failing a run does the appropriate fun stuff."""
    create_table(spark, RUN, add_default_data=add_default_data)
    is_saved_table = add_default_data and pipeline == PIPELINE_RUN[PIPELINE] and source == PIPELINE_RUN[SOURCE]
    pipeline_run = PipelineRun(PIPELINE_RUN[RUN_ID], pipeline, source, TEST_NS)
    error = RuntimeError("ZOMGWTELF! " * 200)

    fail_run(
        spark,
        pipeline_run,
        error=error,
    )
    rows = [r.asDict() for r in spark.table(f"{TEST_NS}.{RUN}").collect()]

    if is_saved_table:
        assert len(rows) == 1
        row = rows[0]
        assert row[RECORDS_PROCESSED] == 0
        assert isinstance(row[END_TIME], datetime.datetime)
        assert row[END_TIME] > row[START_TIME]

        assert row[STATUS] == STATUS_ERROR
        assert row[END_TIME] is not None
        assert len(row[ERROR]) == 1000
        assert len(caplog.records) == 1

    else:
        if add_default_data:
            # only row will be the default data
            assert len(rows) == 1
            check_saved_data(spark, RUN, DEFAULT_DATA[RUN])
        else:
            assert len(rows) == 0

        # two log messages: invalid row, completion
        assert len(caplog.records) == 2
        check_no_record_log_message(pipeline_run, RUN, STATUS_ERROR, caplog.records[0])

    assert caplog.records[-1].levelno == logging.ERROR
    assert caplog.records[-1].message.startswith(
        f"{pipeline} {PIPELINE_RUN[RUN_ID]}: run failed with RuntimeError('ZOMGWTELF! ZOMGWTELF! "
    )
