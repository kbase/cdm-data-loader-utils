"""Tests for the checkpoint audit table functions."""

import datetime
import logging

import pytest
from pyspark.sql import SparkSession

from cdm_data_loader_utils.audit.checkpoint import load_checkpoint, update_checkpoint_status, upsert_checkpoint
from cdm_data_loader_utils.audit.schema import (
    CHECKPOINT,
    LAST_ENTRY_ID,
    PIPELINE,
    RECORDS_PROCESSED,
    RUN_ID,
    SOURCE,
    STATUS,
    STATUS_ERROR,
    STATUS_RUNNING,
    UPDATED,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from tests.audit.conftest import DEFAULT_DATA, check_saved_data, create_table
from tests.conftest import PIPELINE_RUN, TEST_NS


@pytest.mark.requires_spark
def test_upsert_checkpoint_first_entry(spark: SparkSession, pipeline_run: PipelineRun) -> None:
    """Add an entry to the checkpoint table."""
    create_table(spark, CHECKPOINT)

    df = spark.table(f"{TEST_NS}.{CHECKPOINT}")
    assert len(df.collect()) == 0

    # insert an entry
    upsert_checkpoint(
        spark,
        pipeline_run,
        last_entry_id=DEFAULT_DATA[CHECKPOINT][LAST_ENTRY_ID],
        records_processed=DEFAULT_DATA[CHECKPOINT][RECORDS_PROCESSED],
    )
    check_saved_data(spark, CHECKPOINT, DEFAULT_DATA[CHECKPOINT])


@pytest.mark.requires_spark
@pytest.mark.parametrize("pipeline", [PIPELINE_RUN[PIPELINE], "some other pipeline"])
@pytest.mark.parametrize("source", [PIPELINE_RUN[SOURCE], "/some/random/dir"])
def test_upsert_checkpoint_update(spark: SparkSession, pipeline: str, source: str) -> None:
    """Add an entry to the checkpoint table and then update it."""
    df_rows = create_table(spark, CHECKPOINT, add_default_data=True)
    is_saved_table = pipeline == PIPELINE_RUN[PIPELINE] and source == PIPELINE_RUN[SOURCE]
    pipeline_run = PipelineRun(PIPELINE_RUN[RUN_ID], pipeline, source, TEST_NS)

    n_recs = 5
    upsert_checkpoint(
        spark,
        pipeline_run,
        last_entry_id="e2",
        records_processed=n_recs,
    )
    df = spark.table(f"{TEST_NS}.{CHECKPOINT}")
    new_df_data = [r.asDict() for r in df.collect()]

    if is_saved_table:
        # update to the current pipeline: checkpoint row should be updated
        assert len(new_df_data) == 1
        row = new_df_data[0]
    else:
        # different pipeline and source: new checkpoint row.
        assert len(new_df_data) == 2
        assert df_rows[0] in new_df_data
        row = next(r for r in new_df_data if r[LAST_ENTRY_ID] == "e2")

    # check the new data
    assert isinstance(row[UPDATED], datetime.datetime)
    assert row[UPDATED] > df_rows[0][UPDATED]
    assert row[LAST_ENTRY_ID] == "e2"
    assert row[RECORDS_PROCESSED] == n_recs
    assert row[STATUS] == STATUS_RUNNING


@pytest.mark.requires_spark
@pytest.mark.parametrize("pipeline", [PIPELINE_RUN[PIPELINE], "some other pipeline"])
@pytest.mark.parametrize("source", [PIPELINE_RUN[SOURCE], "/path/to/dir"])
@pytest.mark.parametrize("add_default_data", [True, False])
def test_load_checkpoint(spark: SparkSession, pipeline: str, source: str, add_default_data: bool) -> None:  # noqa: FBT001
    """Test loading a previously-saved checkpoint data."""
    is_saved_table = (
        add_default_data
        and pipeline == DEFAULT_DATA[CHECKPOINT][PIPELINE]
        and source == DEFAULT_DATA[CHECKPOINT][SOURCE]
    )
    create_table(spark, CHECKPOINT, add_default_data=add_default_data)
    pipeline_run = PipelineRun(PIPELINE_RUN[RUN_ID], pipeline, source, TEST_NS)

    last_entry = load_checkpoint(spark, pipeline_run)
    if is_saved_table:
        assert last_entry == DEFAULT_DATA[CHECKPOINT][LAST_ENTRY_ID]
    else:
        assert last_entry is None

    # add in a valid entry for the pipeline/data source of interest
    upsert_checkpoint(
        spark,
        pipeline_run,
        last_entry_id="some_entry",
        records_processed=500,
    )
    last_entry = load_checkpoint(spark, pipeline_run)
    assert last_entry == "some_entry"


@pytest.mark.requires_spark
@pytest.mark.parametrize("pipeline", [PIPELINE_RUN[PIPELINE], "pipe_2"])
@pytest.mark.parametrize("source", [PIPELINE_RUN[SOURCE], "/src/two"])
@pytest.mark.parametrize("add_default_data", [True, False])
def test_update_checkpoint_status(
    spark: SparkSession,
    pipeline: str,
    source: str,
    caplog: pytest.LogCaptureFixture,
    add_default_data: bool,  # noqa: FBT001
) -> None:
    """Test updating checkpoint status, with and without existing table data."""
    caplog.set_level(logging.INFO)
    is_saved_table = (
        add_default_data
        and pipeline == DEFAULT_DATA[CHECKPOINT][PIPELINE]
        and source == DEFAULT_DATA[CHECKPOINT][SOURCE]
    )
    create_table(spark, CHECKPOINT, add_default_data=add_default_data)
    pipeline_run = PipelineRun(PIPELINE_RUN[RUN_ID], pipeline, source, TEST_NS)

    update_checkpoint_status(
        spark,
        pipeline_run,
        status=STATUS_ERROR,
    )

    table_rows = spark.table(f"{TEST_NS}.{CHECKPOINT}").collect()

    # only if the pipeline run refers to the default data will the status have been updated
    if is_saved_table:
        assert len(table_rows) == 1
        assert table_rows[0][STATUS] == STATUS_ERROR
    elif add_default_data:
        # this row is from the default data; it won't have been updated
        assert len(table_rows) == 1
        assert table_rows[0][STATUS] == DEFAULT_DATA[CHECKPOINT][STATUS]
    else:
        assert len(table_rows) == 0

    # check log messages
    assert len(caplog.records) == 1
    if is_saved_table:
        assert caplog.records[-1].levelno == logging.INFO
        assert (
            caplog.records[-1].message
            == f"{pipeline} 1234-5678-90: checkpoint successfully updated to status {STATUS_ERROR}"
        )
    else:
        assert caplog.records[-1].levelno == logging.WARNING
        assert (
            caplog.records[-1].message
            == f"{pipeline} 1234-5678-90: cannot update 'checkpoint' to status ERROR because no record exists."
        )
