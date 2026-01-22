"""Audit table for recording data rejected as invalid during ingest."""

import pyspark.sql.functions as sf
from pyspark.sql import DataFrame
from pyspark.sql.types import StructField

from cdm_data_loader_utils.audit.schema import (
    AUDIT_SCHEMA,
    PARSED_ROW,
    PIPELINE,
    RAW_ROW,
    REJECTS,
    ROW_ERRORS,
    RUN_ID,
    SOURCE,
    TIMESTAMP,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger

logger = get_cdm_logger()


def write_rejects(
    run: PipelineRun,
    annotated_df: DataFrame,
    schema_fields: list[StructField],
    invalid_col: str,
) -> None:
    """Write rejected data to the rejects audit table.

    This should plug in directly to readers like the spark CSV reader, which put any non-compliant data into
    a single column when run in PERMISSIVE mode (default for the cdm_data_loader_utils readers).

    It is expected that the dataframe will contain a column called ROW_ERRORS, which contains a list of strings
    describing the errors found in the rows.

    :param run: current pipeline run
    :type run: PipelineRun
    :param annotated_df: dataframe with errors to be written out
    :type annotated_df: DataFrame
    :param schema: schema of the dataframe -- s1hould be the unamended version without a column for invalid data
    :type schema: StructType
    :param invalid_col: name of the column with the invalid data in it
    :type invalid_col: str
    """
    if ROW_ERRORS not in annotated_df.columns:
        err_msg = f"{run.pipeline} {run.run_id}: '{ROW_ERRORS}' column not present in dataframe; cannot record rejects."
        logger.error(err_msg)
        raise RuntimeError(err_msg)

    if annotated_df.count() == 0:
        logger.info("%s %s: nothing to write to '%s' audit table.", run.pipeline, run.run_id, REJECTS)
        return

    # add in a dummy column so that spark doesn't optimise away everything except the error col
    invalid_df: DataFrame = annotated_df.withColumn("_dummy", sf.lit(1)).filter(sf.size(ROW_ERRORS) > 0)
    if not invalid_df.select("_dummy").head(1):
        # nothing to do here
        logger.info("%s %s: nothing to write to '%s' audit table.", run.pipeline, run.run_id, REJECTS)
        return

    data_fields = [f.name for f in schema_fields]
    # drop the dummy column
    invalid_df = invalid_df.drop("_dummy")
    rejects_df = invalid_df.select(
        sf.lit(run.run_id).alias(RUN_ID),
        sf.lit(run.pipeline).alias(PIPELINE),
        sf.lit(run.source_path).alias(SOURCE),
        sf.col(invalid_col).alias(RAW_ROW),
        sf.to_json(sf.struct(*[sf.col(c) for c in data_fields])).alias(PARSED_ROW),
        sf.col(ROW_ERRORS),
        sf.current_timestamp().alias(TIMESTAMP),
    )
    # ensure that it conforms to the schema
    rejects_df = rejects_df.select(
        *[sf.col(f.name).cast(f.dataType).alias(f.name) for f in AUDIT_SCHEMA[REJECTS].fields]
    )
    # write to disk
    rejects_df.write.format("delta").mode("append").saveAsTable(f"{run.namespace}.{REJECTS}")

    get_cdm_logger().info("%s %s: invalid rows written to '%s' audit table.", run.pipeline, run.run_id, REJECTS)
