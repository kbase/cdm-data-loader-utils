"""Run audit table functions: additions and updates to the run table, which tracks overall run status for a pipeline."""

from delta.tables import DeltaTable
from pyspark.sql import SparkSession
from pyspark.sql import functions as sf

from cdm_data_loader_utils.audit.schema import (
    AUDIT_SCHEMA,
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
    STATUS_RUNNING,
    STATUS_SUCCESS,
    match_run,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger


def _table_not_updated(delta: DeltaTable) -> bool:
    metrics = delta.history(1).select("operationMetrics").collect()[0][0]
    # check the number of updated rows, returning True if it is more than 0
    return int(metrics.get("numUpdatedRows", 0)) == 0


def start_run(spark: SparkSession, run: PipelineRun) -> None:
    """Write to the RUN table to indicate that a new ingestion run has started.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: pipeline run
    :type run: PipelineRun
    """
    df = spark.range(1).select(
        sf.lit(run.run_id).alias(RUN_ID),
        sf.lit(run.pipeline).alias(PIPELINE),
        sf.lit(run.source_path).alias(SOURCE),
        sf.lit(STATUS_RUNNING).alias(STATUS),
        sf.lit(sf.lit(0).cast("long")).alias(RECORDS_PROCESSED),
        sf.lit(sf.current_timestamp()).alias(START_TIME),
        sf.lit(sf.lit(None).cast("timestamp")).alias(END_TIME),
        sf.lit(None).cast("string").alias(ERROR),
    )
    df.write.format("delta").mode("append").saveAsTable(f"{run.namespace}.{RUN}")


def complete_run(spark: SparkSession, run: PipelineRun, records_processed: int) -> None:
    """Write to the RUN table to indicate that the ingestion run has completed.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: pipeline run
    :type run: PipelineRun
    :param records_processed: number of records parsed and saved to disk
    :type records_processed: int
    """
    delta: DeltaTable = DeltaTable.forName(spark, f"{run.namespace}.{RUN}")

    delta.update(
        match_run(run),
        {
            STATUS: sf.lit(STATUS_SUCCESS),
            END_TIME: sf.current_timestamp(),
            RECORDS_PROCESSED: sf.lit(records_processed).cast("long"),
        },
    )
    # check whether rows were updated by looking in the delta log
    if _table_not_updated(delta):
        get_cdm_logger().warning(
            "%s %s: cannot update '%s' to status %s because no record exists.",
            run.pipeline,
            run.run_id,
            RUN,
            STATUS_SUCCESS,
        )
    get_cdm_logger().info("%s %s: run completed", run.pipeline, run.run_id)


def fail_run(spark: SparkSession, run: PipelineRun, error: Exception) -> None:
    """Write to the RUN table to indicate that an error has occurred.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: pipeline run
    :type run: PipelineRun
    :param error: error object thrown during ingestion.
    :type error: Exception
    """
    delta = DeltaTable.forName(spark, f"{run.namespace}.{RUN}")
    delta.update(
        match_run(run),
        set={
            END_TIME: sf.current_timestamp(),
            STATUS: sf.lit(STATUS_ERROR),
            ERROR: sf.lit(str(error)[:1000]),
        },
    )
    # check whether rows were updated by looking in the delta log
    if _table_not_updated(delta):
        get_cdm_logger().warning(
            "%s %s: cannot update '%s' to status %s because no record exists.",
            run.pipeline,
            run.run_id,
            RUN,
            STATUS_ERROR,
        )
    get_cdm_logger().error("%s %s: run failed with %s", run.pipeline, run.run_id, error.__repr__())
