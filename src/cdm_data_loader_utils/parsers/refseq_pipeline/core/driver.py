from typing import List, Dict, Any
from pyspark.sql import SparkSession
from refseq_pipeline.core.cdm_parse import parse_reports
from refseq_pipeline.core.spark_delta import write_delta_table, cleanup_after_write


def process_and_write_reports(
    spark: SparkSession,
    reports: List[Dict[str, Any]],
    database: str,
    table: str = "assembly_stats",
    mode: str = "append",
    data_dir: str | None = None,
    prefer_spark: bool = True,
    optimize: bool = False,
    vacuum: bool = False
):
    """
    Parse a list of assembly reports and write them to Delta Lake.

    Args:
        spark: SparkSession
        reports: list of raw JSON records
        database: target Delta database
        table: target Delta table name
        mode: 'append' or 'overwrite'
        data_dir: external Delta path (optional)
        prefer_spark: if True, use Spark-native pipeline
        optimize: whether to OPTIMIZE ZORDER the table
        vacuum: whether to VACUUM old files
    """
    if not reports:
        print("[driver] No reports to process.")
        return

    df = parse_reports(reports, return_spark=prefer_spark, spark=spark)

    write_delta_table(
        sdf=df,
        spark=spark,
        database=database,
        table=table,
        mode=mode,
        data_dir=data_dir,
    )

    cleanup_after_write(
        spark=spark,
        database=database,
        table=table,
        do_optimize=optimize,
        do_vacuum=vacuum,
    )
