"""Tests for the checkpoint and run audit table functions."""

import datetime
from typing import Any

from pyspark.sql import SparkSession

from cdm_data_loaders.audit.schema import (
    AUDIT_SCHEMA,
    CHECKPOINT,
    END_TIME,
    ERROR,
    LAST_ENTRY_ID,
    METRICS,
    N_INVALID,
    N_READ,
    N_VALID,
    RECORDS_PROCESSED,
    REJECTS,
    RUN,
    START_TIME,
    STATUS,
    TIMESTAMP,
    UPDATED,
    VALIDATION_ERRORS,
)
from tests.conftest import PIPELINE_RUN, TEST_NS

DEFAULT_DATA = {
    CHECKPOINT: {
        **PIPELINE_RUN,
        STATUS: "RUNNING",
        RECORDS_PROCESSED: 1,
        LAST_ENTRY_ID: "e1",
        # UPDATED
    },
    RUN: {
        **PIPELINE_RUN,
        STATUS: "RUNNING",
        RECORDS_PROCESSED: 0,
        END_TIME: None,
        ERROR: None,
        # START_TIME
    },
    METRICS: {
        **PIPELINE_RUN,
        N_READ: 10,
        N_VALID: 8,
        N_INVALID: 2,
        VALIDATION_ERRORS: ["a mess o' trouble", "another fine mess"],
        # UPDATED
    },
    REJECTS: {
        **PIPELINE_RUN,
    },
}

INIT_TIMESTAMP_FIELDS = {CHECKPOINT: [UPDATED], RUN: [START_TIME], METRICS: [UPDATED], REJECTS: [TIMESTAMP]}
END_TIMESTAMP_FIELDS = {RUN: [END_TIME]}


def create_table(spark: SparkSession, table_name: str, add_default_data: bool = False) -> list[dict[str, Any]]:
    """Helper for creating audit tables for tests."""
    # ensure the db exists
    if table_name not in AUDIT_SCHEMA:
        msg = f"Invalid audit table: {table_name}"
        raise ValueError(msg)

    spark.sql(f"CREATE DATABASE IF NOT EXISTS {TEST_NS}")

    # add in default data, if there is any
    starter_data = []
    if add_default_data and table_name in (CHECKPOINT, METRICS, RUN):
        starter_dict = {**DEFAULT_DATA[table_name]}
        # add in any fields that are set when the table is added
        for f in INIT_TIMESTAMP_FIELDS.get(table_name, []):
            starter_dict[f] = datetime.datetime.now(tz=datetime.UTC)
        starter_data = [starter_dict]

    spark.createDataFrame(starter_data, AUDIT_SCHEMA[table_name]).write.format("delta").mode("error").saveAsTable(
        f"{TEST_NS}.{table_name}"
    )
    # validate
    df = spark.table(f"{TEST_NS}.{table_name}")
    if add_default_data:
        assert df.count() == 1
        check_saved_data(spark, table_name, DEFAULT_DATA[table_name])
    else:
        assert df.count() == 0

    return [r.asDict() for r in df.collect()]


def check_saved_data(spark: SparkSession, table_name: str, expected: dict[str, Any]) -> list[dict[str, Any]]:
    """Ensure the appropriate checkpoint data has been saved.

    :param spark: spark sesh
    :type spark: SparkSession
    :param table_name: the name of the table to check
    :type table_name: str
    :param expected: dictionary of expected data
    :type expected: dict[str, Any]
    :return: checkpoint data obj from the db
    :rtype: dict[str, Any]
    """
    df = spark.table(f"{TEST_NS}.{table_name}")
    assert len(df.collect()) == 1
    df_row = df.collect()[0].asDict()
    # check and then remove the timestamp fields
    for f in INIT_TIMESTAMP_FIELDS.get(table_name, []):
        assert isinstance(df_row[f], datetime.datetime)
        del df_row[f]
    # special case: end time may or may not be populated
    # if table_name == RUN and
    # the remaining fields should be checkable
    assert df_row == expected

    return [r.asDict() for r in df.collect()]
