"""Schema for audit tables."""

from pyspark.sql.types import ArrayType, IntegerType, LongType, StringType, StructField, StructType, TimestampType

from cdm_data_loader_utils.core.pipeline_run import PipelineRun

CHECKPOINT = "checkpoint"
METRICS = "metrics"
REJECTS = "rejects"
RUN = "run"

RUN_ID = "run_id"
PIPELINE = "pipeline"
SOURCE = "source_path"
NAMESPACE = "namespace"

END_TIME = "end_time"
ERROR = "error"
LAST_ENTRY_ID = "last_entry_id"
PARSED_ROW = "parsed_record"
RAW_ROW = "raw_record"
RECORDS_PROCESSED = "records_processed"
START_TIME = "start_time"
STATUS = "status"
TABLE = "table"
TIMESTAMP = "timestamp"
UPDATED = "updated"

STATUS_ERROR = "ERROR"
STATUS_RUNNING = "RUNNING"
STATUS_SUCCESS = "SUCCESS"

N_INVALID = "records_invalid"
N_READ = "records_read"
N_VALID = "records_valid"
ROW_ERRORS = "errors_in_record"
VALIDATION_ERRORS = "validation_errors"


FIELDS = {
    RUN_ID: StructField(RUN_ID, StringType(), nullable=False),
    PIPELINE: StructField(PIPELINE, StringType(), nullable=False),
    SOURCE: StructField(SOURCE, StringType(), nullable=False),
    STATUS: StructField(STATUS, StringType(), nullable=False),
    UPDATED: StructField(UPDATED, TimestampType(), nullable=False),
}

PIPELINE_FIELDS = [FIELDS[RUN_ID], FIELDS[PIPELINE], FIELDS[SOURCE]]

AUDIT_SCHEMA = {
    # represents the whole run
    RUN: StructType(
        [
            *PIPELINE_FIELDS,
            FIELDS[STATUS],
            StructField(RECORDS_PROCESSED, LongType(), nullable=True),
            StructField(START_TIME, TimestampType(), nullable=False),
            StructField(END_TIME, TimestampType(), nullable=True),
            StructField(ERROR, StringType(), nullable=True),
        ]
    ),
    # checkpoint during processing to allow parsing to be restarted after an error
    CHECKPOINT: StructType(
        [
            *PIPELINE_FIELDS,
            FIELDS[STATUS],
            StructField(RECORDS_PROCESSED, LongType(), nullable=True),
            StructField(LAST_ENTRY_ID, StringType(), nullable=True),
            FIELDS[UPDATED],
        ]
    ),
    # per-source metrics on data ingested
    METRICS: StructType(
        [
            *PIPELINE_FIELDS,
            StructField(N_READ, IntegerType(), nullable=False),
            StructField(N_VALID, IntegerType(), nullable=False),
            StructField(N_INVALID, IntegerType(), nullable=False),
            StructField(VALIDATION_ERRORS, ArrayType(StringType(), containsNull=False), nullable=False),
            FIELDS[UPDATED],
        ]
    ),
    # invalid data
    REJECTS: StructType(
        [
            *PIPELINE_FIELDS,
            StructField(RAW_ROW, StringType(), nullable=False),
            StructField(PARSED_ROW, StringType(), nullable=True),
            StructField(ROW_ERRORS, ArrayType(StringType()), nullable=False),
            StructField(TIMESTAMP, TimestampType(), nullable=False),
        ]
    ),
}


def current_run_expr(target: str = "t", source: str = "s") -> str:
    """SQL expression for matching the current run in any of the audit tables."""
    return " AND ".join([f"{target or 't'}.{param} = {source or 's'}.{param}" for param in [RUN_ID, SOURCE, PIPELINE]])


def match_run(run: PipelineRun) -> str:
    """SQL expression for matching the current run in any of the audit tables.

    :param run: current pipeline run
    :type run: PipelineRun
    :return: SQL expression string
    :rtype: str
    """
    return " AND ".join([f"{p} = '{getattr(run, p)}'" for p in [RUN_ID, PIPELINE, SOURCE]])
