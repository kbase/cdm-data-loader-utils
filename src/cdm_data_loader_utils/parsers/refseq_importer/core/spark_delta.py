import os
import shutil
from typing import Optional
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.types import StructType, StructField, StringType


def build_spark(database: str) -> SparkSession:
    """
    Initialize a Spark session with Delta Lake support and create the specified database if it doesn't exist.
    """
    builder = (
        SparkSession.builder.appName("NCBI Datasets -> CDM")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )
    spark = configure_spark_with_delta_pip(builder).getOrCreate()

    # Create the database namespace if it doesn't exist
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {database}")
    return spark


def write_delta(
    spark: SparkSession, df: DataFrame, database: str, table: str, mode: str = "append", data_dir: Optional[str] = None
) -> None:
    """
    Write Spark DataFrame to Delta.
    If data_dir is provided, writes to external LOCATION {data_dir}/{database}/{table}.
    Otherwise writes to managed table.

    """

    if df is None or df.rdd.isEmpty():
        print(f"No data to write to {database}.{table}")
        return

    print(f"Writing table={table}, rows={df.count()}")
    df.printSchema()
    df.show(10, truncate=False)

    # Special schema case
    if table == "contig_collection":
        schema = StructType(
            [
                StructField("collection_id", StringType(), True),
                StructField("contig_collection_type", StringType(), True),
                StructField("ncbi_taxon_id", StringType(), True),
                StructField("gtdb_taxon_id", StringType(), True),
            ]
        )
        df = spark.createDataFrame(df.rdd, schema=schema)

    writer = df.write.format("delta").mode(mode)
    writer = writer.option("mergeSchema", "true") if mode == "append" else writer.option("overwriteSchema", "true")

    full_table = f"{database}.{table}"

    if data_dir:
        target_path = os.path.abspath(os.path.join(data_dir, database, table))
        os.makedirs(target_path, exist_ok=True)

        delta_log = os.path.join(target_path, "_delta_log")
        if mode == "overwrite" and os.path.exists(target_path) and not os.path.exists(delta_log):
            print(f"[WARN] Non-Delta data at {target_path}")
            shutil.rmtree(target_path, ignore_errors=True)
            os.makedirs(target_path, exist_ok=True)

        writer.save(target_path)

        # Register/create an external table using LOCATION
        spark.sql(f"""
            CREATE TABLE IF NOT EXISTS {full_table}
            USING DELTA
            LOCATION '{target_path}'
        """)

        print(f"Saved table {full_table} (rows={df.count()}) -> {target_path}")

    else:
        writer.saveAsTable(full_table)
        print(f"Saved managed table {full_table} (rows={df.count()})")


def preview_or_skip(spark: SparkSession, database: str, table: str, limit: int = 20) -> None:
    """
    Preview table if it exists.
    """

    full_table = f"{database}.{table}"
    if spark.catalog.tableExists(full_table):
        print(f"Preview for {full_table}:")
        spark.sql(f"SELECT * FROM {full_table} LIMIT {limit}").show(truncate=False)
    else:
        print(f"Table {full_table} not found. Skipping preview.")
