"""Utilities for running tests."""

from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.types import StructType


def create_empty_delta_table(
    spark: SparkSession,
    db: str,
    table: str,
    schema: StructType,
) -> None:
    """Create an empty delta table, initialising the db namespace first."""
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")
    df = spark.createDataFrame([], schema)
    df.write.format("delta").mode("error").saveAsTable(f"{db}.{table}")
    assert df.count() == 0
