"""
PySpark-based normalization pipeline for GO Gene Association Files (GAF).

This script processes a raw GAF-like annotation CSV (e.g., `annotations_data100.csv`)
and produces a normalized output with a schema consistent with the previous
Pandas-based implementation.

Two output formats are supported:

1) CSV Output:
    python3 association_update.py \
        --input annotations_data100.csv \
        --output normalized_annotation_update.csv

2) Delta Lake Output:
    python3 association_update.py \
        --input annotations_data100.csv \
        --output ./delta_output

Result:
- A Parquet-backed Delta table containing normalized annotations.
- Schema conforms to the CDM-style structured annotation model.
"""

import logging
import os
import sys
import urllib.request

import click
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.functions import (
    col,
    concat,
    concat_ws,
    explode,
    lit,
    regexp_replace,
    split,
    to_date,
    trim,
    upper,
    when,
)
from pyspark.sql.types import StringType

# ---------------------- Logging Setup ----------------------
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# --- Constants ---
SUBJECT = "subject"
PREDICATE = "predicate"
OBJECT = "object"
PUBLICATIONS = "publications"
EVIDENCE_CODE = "Evidence_Code"
SUPPORTING_OBJECTS = "supporting_objects"
ANNOTATION_DATE = "annotation_date"
PRIMARY_KNOWLEDGE_SOURCE = "primary_knowledge_source"
AGGREGATOR = "aggregator"
PROTOCOL_ID = "protocol_id"
NEGATED = "negated"
EVIDENCE_TYPE = "evidence_type"

# GAF Field Names
DB = "DB"
DB_OBJ_ID = "DB_Object_ID"
QUALIFIER = "Qualifier"
GO_ID = "GO_ID"
DB_REF = "DB_Reference"
WITH_FROM = "With_From"
DATE = "Date"
ASSIGNED_BY = "Assigned_By"

# ECO Mapping
ECO_MAPPING_URL = "http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt"

ALLOWED_PREDICATES = [
    "enables",
    "contributes_to",
    "acts_upstream_of_or_within",
    "involved_in",
    "acts_upstream_of",
    "acts_upstream_of_positive_effect",
    "acts_upstream_of_negative_effect",
    "acts_upstream_of_or_within_negative_effect",
    "acts_upstream_of_or_within_positive_effect",
    "located_in",
    "part_of",
    "is_active_in",
    "colocalizes_with",
]


def get_spark():
    """Initialize and return a Spark session configured for Delta Lake."""
    builder = (
        SparkSession.builder.appName("GO-GAF-Spark-Parser")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
        .config("spark.sql.shuffle.partitions", "200")
    )
    return configure_spark_with_delta_pip(builder).getOrCreate()


def load_annotation(spark, input_path):
    """Load and preprocess raw annotation CSV."""
    df = spark.read.csv(input_path, header=True, inferSchema=True)

    df = df.select(DB, DB_OBJ_ID, QUALIFIER, GO_ID, DB_REF, EVIDENCE_CODE, WITH_FROM, DATE, ASSIGNED_BY)

    df = (
        df.withColumn(PREDICATE, col(QUALIFIER))
        .withColumn(OBJECT, col(GO_ID))
        .withColumn(PUBLICATIONS, split(trim(when(col(DB_REF).isNotNull(), col(DB_REF)).otherwise(lit(""))), "\\|"))
        .withColumn(SUPPORTING_OBJECTS, split(trim(col(WITH_FROM)), "\\|"))
        .withColumn(ANNOTATION_DATE, col(DATE))
        .withColumn(PRIMARY_KNOWLEDGE_SOURCE, col(ASSIGNED_BY))
    )

    return df


def normalize_dates(df):
    """Normalize annotation dates to yyyy-MM-dd format if 8-digit string."""
    df = df.withColumn(
        ANNOTATION_DATE, when(col(ANNOTATION_DATE).rlike("^[0-9]{8}$"), to_date(col(ANNOTATION_DATE), "yyyyMMdd"))
    )
    return df


def process_predicates(df):
    """Validate and clean predicate values (e.g., remove NOT| prefix)."""
    df = df.withColumn(NEGATED, col(PREDICATE).startswith("NOT|")).withColumn(
        PREDICATE, regexp_replace(col(PREDICATE), "^NOT\\|", "")
    )

    invalid = df.filter(~col(PREDICATE).isin(ALLOWED_PREDICATES))
    if invalid.count() > 0:
        invalid_values = [r[PREDICATE] for r in invalid.select(PREDICATE).distinct().collect()]
        raise ValueError(f"Invalid predicate found {invalid_values}")
    return df


def add_metadata(df):
    """Add aggregator, protocol ID, and subject URI."""
    return (
        df.withColumn(AGGREGATOR, lit("UniProt"))
        .withColumn(PROTOCOL_ID, lit(None).cast(StringType()))
        .withColumn(SUBJECT, concat(col(DB).cast("string"), lit(":"), col(DB_OBJ_ID).cast("string")))
    )


def load_eco_mapping(spark, local_path="gaf-eco-mapping.txt"):
    """Download and load ECO evidence mapping table."""
    if not os.path.exists(local_path):
        print(f"Downloading ECO mapping file to: {local_path}")
        urllib.request.urlretrieve(ECO_MAPPING_URL, local_path)

    df = spark.read.csv(local_path, sep="\t", comment="#", header=False)
    return df.toDF(EVIDENCE_CODE, DB_REF, EVIDENCE_TYPE)


def merge_evidence(df, eco):
    """Join annotation DataFrame with ECO evidence mapping."""
    df = (
        df.withColumn(PUBLICATIONS, explode(col(PUBLICATIONS)))
        .filter(col(PUBLICATIONS).isNotNull() & (col(PUBLICATIONS) != ""))
        .withColumn(PUBLICATIONS, upper(trim(col(PUBLICATIONS))))
    )

    eco = eco.withColumn(DB_REF, upper(trim(col(DB_REF)))).withColumn(EVIDENCE_CODE, upper(trim(col(EVIDENCE_CODE))))

    merged = (
        df.alias("df")
        .join(
            eco.alias("eco"),
            on=(col("df." + EVIDENCE_CODE) == col("eco." + EVIDENCE_CODE))
            & (col("df." + PUBLICATIONS) == col("eco." + DB_REF)),
            how="left",
        )
        .drop(col("eco." + DB_REF))
        .drop(col("eco." + EVIDENCE_CODE))
    )

    fallback = (
        eco.filter(col(DB_REF) == "DEFAULT")
        .select(EVIDENCE_CODE, EVIDENCE_TYPE)
        .withColumnRenamed(EVIDENCE_TYPE, "fallback")
    )

    merged = (
        merged.join(fallback, on=EVIDENCE_CODE, how="left")
        .withColumn(EVIDENCE_TYPE, when(col(EVIDENCE_TYPE).isNull(), col("fallback")).otherwise(col(EVIDENCE_TYPE)))
        .drop("fallback")
    )

    return merged


def reorder_columns(df):
    """Ensure correct column order and clean up types."""
    df = (
        df.withColumn(PUBLICATIONS, concat_ws("|", col(PUBLICATIONS)))
        .withColumn(SUPPORTING_OBJECTS, concat_ws("|", col(SUPPORTING_OBJECTS)))
        .withColumn(SUPPORTING_OBJECTS, when(col(SUPPORTING_OBJECTS) == "", None).otherwise(col(SUPPORTING_OBJECTS)))
        .withColumn(NEGATED, col(NEGATED).cast("boolean").cast("string"))
    )

    final_cols = [
        OBJECT,
        DB,
        ANNOTATION_DATE,
        PREDICATE,
        EVIDENCE_CODE,
        PUBLICATIONS,
        DB_OBJ_ID,
        PRIMARY_KNOWLEDGE_SOURCE,
        SUPPORTING_OBJECTS,
        AGGREGATOR,
        PROTOCOL_ID,
        NEGATED,
        SUBJECT,
        EVIDENCE_TYPE,
    ]
    return df.select([col(c) for c in final_cols])


def write_output(df, output_path, mode="overwrite"):
    df.write.format("delta").mode(mode).save(output_path)


def register_table(spark, output_path, table_name="normalized_annotation", permanent=True):
    if permanent:
        logger.info(f"Registering Delta table as permanent table {table_name}")

        spark.sql(f"""
            CREATE TABLE IF NOT EXISTS {table_name}
            USING DELTA
            LOCATION '{output_path}'
        """)

    else:
        logger.info(f"Registering Delta table as temporary view: {table_name}")
        df = spark.read.format("delta").load(output_path)
        df.createOrReplaceTempView(table_name)


def run(
    input_path,
    output_path,
    register=False,
    table_name="normalized_annotation",
    permanent=True,
    dry_run=False,
    mode="overwrite",
):
    spark = None
    try:
        spark = get_spark()
        logger.info("Starting annotation pipeline")

        eco = load_eco_mapping(spark)
        df = load_annotation(spark, input_path)
        df = normalize_dates(df)
        df = process_predicates(df)
        df = add_metadata(df)
        df = merge_evidence(df, eco)
        df = reorder_columns(df)

        if dry_run:
            logger.info("showing top 5 rows")
            df.show(5, truncate=False)
        else:
            write_output(df, output_path, mode=mode)
            logger.info(f"Data written to {output_path}")
            if register:
                register_table(spark, output_path, table_name=table_name, permanent=permanent)

    except Exception as e:
        logger.error(f"Pipeline failed: {e}")
        sys.exit(1)
    finally:
        if spark:
            spark.stop()


@click.command()
@click.option("--input", "-i", required=True, help="Path to input CSV file")
@click.option("--output", "-o", required=True, help="Target Delta table output directory")
@click.option("--register", is_flag=True, help="Register the output as Spark SQL table")
@click.option("--table-name", default="normalized_annotation", help="SQL table name to register")
@click.option("--temp", is_flag=True, help="Register as temporary view (default is permanent)")
@click.option(
    "--mode", default="overwrite", type=click.Choice(["overwrite", "append", "ignore"]), help="Delta write mode"
)
@click.option("--dry-run", is_flag=True, help="Dry run without writing output")
def main(input, output, register, table_name, temp, mode, dry_run):
    if not os.path.isfile(input):
        logger.error(f"Input file does not exist: {input}")
        sys.exit(1)
    run(input, output, register=register, table_name=table_name, permanent=not temp, dry_run=dry_run, mode=mode)


if __name__ == "__main__":
    main()
