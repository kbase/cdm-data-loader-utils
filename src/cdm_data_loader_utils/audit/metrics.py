"""Audit table for recording metrics."""

from delta.tables import DeltaTable
from pyspark.sql import DataFrame, Row, SparkSession
from pyspark.sql import functions as sf

from cdm_data_loader_utils.audit.schema import (
    AUDIT_SCHEMA,
    METRICS,
    N_INVALID,
    N_READ,
    N_VALID,
    PIPELINE,
    ROW_ERRORS,
    RUN_ID,
    SOURCE,
    UPDATED,
    VALIDATION_ERRORS,
    current_run_expr,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger

logger = get_cdm_logger()


def write_metrics(
    spark: SparkSession,
    annotated_df: DataFrame,
    run: PipelineRun,
) -> Row:
    """Write metrics for the current batch of imports to disk.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: current pipeline run
    :type run: PipelineRun
    :param records_read: total number of records read
    :type records_read: int
    :param records_valid: number of valid records
    :type records_valid: int
    :param records_invalid: number of invalid records
    :type records_invalid: int
    :param validation_errors: list of validation errors encountered in the batch
    :type validation_errors: list[str]
    :return: row of a dataframe with validation metrics
    :rtype: Row
    """
    if annotated_df.count() == 0:
        # nothing to do here
        logger.info("%s %s: nothing to write to '%s' audit table.", run.pipeline, run.run_id, METRICS)
        return Row(records_read=0, records_valid=0, records_invalid=0, validation_errors=[])

    invalid_df = annotated_df.filter(sf.size(ROW_ERRORS) > 0)

    validation_errors = sorted(
        [r.reason for r in invalid_df.select(sf.explode(ROW_ERRORS).alias("reason")).distinct().collect()]
    )

    metrics_df = annotated_df.agg(
        sf.count("*").alias(N_READ),
        sf.sum(sf.when(sf.size(ROW_ERRORS) == 0, 1).otherwise(0)).alias(N_VALID),
        sf.sum(sf.when(sf.size(ROW_ERRORS) > 0, 1).otherwise(0)).alias(N_INVALID),
        sf.lit(validation_errors).alias(VALIDATION_ERRORS),
    )
    metrics = metrics_df.collect()[0]

    df = spark.range(1).select(
        sf.lit(run.run_id).alias(RUN_ID),
        sf.lit(run.pipeline).alias(PIPELINE),
        sf.lit(run.source_path).alias(SOURCE),
        sf.lit(metrics.records_read).alias(N_READ),
        sf.lit(metrics.records_valid).alias(N_VALID),
        sf.lit(metrics.records_invalid).alias(N_INVALID),
        sf.lit(metrics.validation_errors).alias(VALIDATION_ERRORS),
        sf.current_timestamp().alias(UPDATED),
    )

    target = DeltaTable.forName(
        spark,
        f"{run.namespace}.{METRICS}",
    )

    (
        target.alias("t")
        .merge(
            df.alias("s"),
            current_run_expr(),
        )
        .whenMatchedUpdate(set={k: f"s.{k}" for k in [N_READ, N_VALID, N_INVALID, VALIDATION_ERRORS, UPDATED]})
        .whenNotMatchedInsertAll()
        .execute()
    )
    get_cdm_logger().info("%s %s: ingest metrics written to '%s' table.", run.pipeline, run.run_id, METRICS)

    return metrics
