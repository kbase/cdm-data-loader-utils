"""Checkpoint audit table functions: adding and updating information on data import pipeline execution."""

from delta.tables import DeltaTable
from pyspark.sql import SparkSession
from pyspark.sql import functions as sf

from cdm_data_loader_utils.audit.schema import (
    AUDIT_SCHEMA,
    CHECKPOINT,
    LAST_ENTRY_ID,
    PIPELINE,
    RECORDS_PROCESSED,
    RUN_ID,
    SOURCE,
    STATUS,
    STATUS_RUNNING,
    UPDATED,
    current_run_expr,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger


# Checkpoint table-related functions
def upsert_checkpoint(
    spark: SparkSession,
    run: PipelineRun,
    last_entry_id: str,
    records_processed: int,
) -> None:
    """Add or update checkpoint records for the current pipeline ingest.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: pipeline run
    :type run: PipelineRun
    :param last_entry_id: ID of the last entry processed
    :type last_entry_id: str
    :param records_processed: number of entries processed
    :type records_processed: int
    """
    delta = DeltaTable.forName(spark, f"{run.namespace}.{CHECKPOINT}")

    df = spark.range(1).select(
        sf.lit(run.run_id).alias(RUN_ID),
        sf.lit(run.pipeline).alias(PIPELINE),
        sf.lit(run.source_path).alias(SOURCE),
        sf.lit(STATUS_RUNNING).alias(STATUS),
        sf.lit(records_processed).alias(RECORDS_PROCESSED),
        sf.lit(last_entry_id).alias(LAST_ENTRY_ID),
        sf.current_timestamp().alias(UPDATED),
    )
    updates = spark.createDataFrame(df.rdd, schema=AUDIT_SCHEMA[CHECKPOINT])

    (
        delta.alias("t")
        .merge(updates.alias("s"), current_run_expr())
        .whenMatchedUpdate(set={val: f"s.{val}" for val in [STATUS, RECORDS_PROCESSED, LAST_ENTRY_ID, UPDATED]})
        .whenNotMatchedInsertAll()
        .execute()
    )
    get_cdm_logger().info("%s %s: checkpoint created/updated", run.pipeline, run.run_id)


def update_checkpoint_status(spark: SparkSession, run: PipelineRun, status: str) -> None:
    """Update checkpoint status and timestamp.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: pipeline run
    :type run: PipelineRun
    :param status: pipeline status
    :type status: str
    """
    delta = DeltaTable.forName(spark, f"{run.namespace}.{CHECKPOINT}")
    delta.update(
        " AND ".join([f"{p} = '{getattr(run, p)}'" for p in [RUN_ID, PIPELINE, SOURCE]]),
        {STATUS: sf.lit(status), UPDATED: sf.current_timestamp()},
    )
    # check whether rows were updated by looking in the delta log
    # N.b. this may not work correctly if another process updates the table in the interim
    metrics = delta.history(1).select("operationMetrics").collect()[0][0]
    if int(metrics.get("numUpdatedRows", 0)) == 0:
        get_cdm_logger().warning(
            "%s %s: cannot update '%s' to status %s because no record exists.",
            run.pipeline,
            run.run_id,
            CHECKPOINT,
            status,
        )
    else:
        get_cdm_logger().info("%s %s: checkpoint successfully updated to status %s", run.pipeline, run.run_id, status)


def load_checkpoint(spark: SparkSession, run: PipelineRun) -> str | None:
    """Load any existing checkpoint data, filtered by current pipeline and data source path.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: pipeline run
    :type run: PipelineRun
    :return: either the last entry successfully saved and dumped to disk, or None
    :rtype: str | None
    """
    rows = (
        spark.table(f"{run.namespace}.{CHECKPOINT}")
        .filter(sf.col(RUN_ID) == run.run_id)
        .filter(sf.col(PIPELINE) == run.pipeline)
        .filter(sf.col(SOURCE) == run.source_path)
        .select(LAST_ENTRY_ID)
        .limit(1)
        .collect()
    )
    return rows[0][LAST_ENTRY_ID] if rows else None
