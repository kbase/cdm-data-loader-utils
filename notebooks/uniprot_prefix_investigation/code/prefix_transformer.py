"""
Prefix Transformer

Reads identifier parquet,
normalizes db prefixes using prefix_normalizer,
writes transformed parquet.

Usage:
PYTHONPATH=src python src/cdm_data_loader_utils/parsers/prefix_transformer.py \
    --input dataset/input.parquet \
    --output dataset_normalized \
    --registry registry.json \
    --strict
"""

from __future__ import annotations

import argparse
import json
from typing import Any

from pyspark.sql import SparkSession
from pyspark.sql.functions import udf, col
from pyspark.sql.types import (
    StructType,
    StructField,
    StringType,
    BooleanType,
)

from cdm_data_loader_utils.parsers.prefix_normalizer import normalize_prefix


# ----------------------------------
# CLI
# ----------------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prefix Transformer")

    parser.add_argument("--input", required=True, help="Input parquet path")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--registry", required=True, help="BioRegistry JSON file")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on unclassified prefixes",
    )

    return parser.parse_args()


# ----------------------------------
# Main
# ----------------------------------
def main() -> None:
    args = parse_args()

    spark = SparkSession.builder.appName("PrefixTransformer").getOrCreate()
    spark.sparkContext.setLogLevel("WARN")

    # ----------------------------------
    # Load BioRegistry
    # ----------------------------------
    with open(args.registry) as f:
        registry = json.load(f)

    registry_set = {k.lower() for k in registry.keys()}

    # broadcast registry to executors
    registry_bc = spark.sparkContext.broadcast(registry_set)

    # ----------------------------------
    # Define UDF
    # ----------------------------------
    def normalize_wrapper(db_value: str | None) -> tuple[Any, Any, Any, Any]:
        if db_value is None:
            return (None, None, "null", False)

        result = normalize_prefix(db_value, registry_bc.value)

        category = result["category"]
        normalized = result["normalized"]
        is_gap = result["is_registry_gap"]

        if category == "annotation":
            return (db_value, None, "annotation", False)

        if args.strict and category not in {
            "exact",
            "map",
            "synonym",
            "registry_gap",
        }:
            raise ValueError(f"Unclassified prefix detected: {db_value}")

        return (db_value, normalized, category, is_gap)

    schema = StructType(
        [
            StructField("db_original", StringType(), True),
            StructField("db_normalized", StringType(), True),
            StructField("prefix_category", StringType(), True),
            StructField("is_registry_gap", BooleanType(), True),
        ]
    )

    normalize_udf = udf(normalize_wrapper, schema)

    # ----------------------------------
    # Read Input
    # ----------------------------------
    df = spark.read.parquet(args.input)

    # ----------------------------------
    # Apply Transformation
    # ----------------------------------
    df_transformed = df.withColumn("prefix_struct", normalize_udf(col("db")))

    df_transformed = df_transformed.select(
        "*",
        col("prefix_struct.db_original"),
        col("prefix_struct.db_normalized"),
        col("prefix_struct.prefix_category"),
        col("prefix_struct.is_registry_gap"),
    ).drop("prefix_struct")

    # optional: remove annotation rows
    df_transformed = df_transformed.filter(col("prefix_category") != "annotation")

    # ----------------------------------
    # Write Output
    # ----------------------------------
    df_transformed.write.mode("overwrite").parquet(args.output)

    print("Transformation complete.")
    spark.stop()


if __name__ == "__main__":
    main()
