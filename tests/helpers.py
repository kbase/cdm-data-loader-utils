"""Utilities for running tests."""

from pyspark.sql import Row, SparkSession
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


def assertDataFrameEqual(result_rows: list[Row], expected_rows: list[Row]) -> None:  # noqa: N802
    """Workaround for assertDataFrameEqual from pyspark.testing being broken by pandas 3.0.

    :param result_df: list of dataframe rows, as returned by df.collect()
    :type result_df: list[Row]
    :param expected_df: expected dataframe rows, as returned by df.collect()
    :type expected_df: list[Row]
    """
    results_dict = [r.asDict() for r in result_rows]
    expected_dict = [r.asDict() for r in expected_rows]
    assert len(results_dict) == len(expected_dict)
    for row in results_dict:
        assert row in expected_dict
    for row in expected_dict:
        assert row in results_dict
