"""Utilities for running tests."""

from pyspark.sql import DataFrame, SparkSession
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


def assertDataFrameEqual(result_df: DataFrame, expected_df: DataFrame) -> None:
    """Workaround for assertDataFrameEqual from pyspark.testing being broken by pandas 3.0.

    :param result_df: result dataframe
    :type result_df: DataFrame
    :param expected_df: expected dataframe
    :type expected_df: DataFrame
    """
    results_dict = [r.asDict() for r in result_df]
    expected_dict = [r.asDict() for r in expected_df]
    assert len(results_dict) == len(expected_dict)
    for row in results_dict:
        assert row in expected_dict
    for row in expected_dict:
        assert row in results_dict
