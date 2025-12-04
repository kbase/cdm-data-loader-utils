import pytest
from pyspark.sql import SparkSession
from pyspark.sql import Row
import shutil
import os
from pyspark.sql.types import StructType, StructField, StringType
from refseq_importer.core.spark_delta import (
    build_spark,
    write_delta,
    preview_or_skip,
)



# =============================================================
# Spark fixture
# =============================================================
@pytest.fixture(scope="session")
def spark():
    spark = (
        SparkSession.builder
        .master("local[1]")
        .appName("spark-delta-test")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
        .getOrCreate()
    )
    yield spark
    spark.stop()


# =============================================================
# build_spark
# =============================================================

def test_build_spark_creates_database(tmp_path):
    db = "testdb"
    spark = build_spark(db)
    dbs = [d.name for d in spark.catalog.listDatabases()]
    assert db in dbs


# =============================================================
# write_delta (managed table)
# =============================================================

def test_write_delta_managed_table(spark):
    db = "writetest"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")

    df = spark.createDataFrame([Row(a="X", b="Y")])

    write_delta(
        spark=spark,
        df=df,
        database=db,
        table="example",
        mode="overwrite",
        data_dir=None,
    )

    # Table should exist
    assert spark.catalog.tableExists(f"{db}.example")

    # Data should exist
    rows = spark.sql(f"SELECT a, b FROM {db}.example").collect()
    assert rows[0]["a"] == "X"
    assert rows[0]["b"] == "Y"


# =============================================================
# write_delta with external LOCATION
# =============================================================

def test_write_delta_external_location(spark, tmp_path):
    db = "externaldb"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")

    df = spark.createDataFrame([Row(id="1", value="A")])

    write_delta(
        spark=spark,
        df=df,
        database=db,
        table="exttable",
        mode="overwrite",
        data_dir=str(tmp_path),
    )

    # Table should be registered
    assert spark.catalog.tableExists(f"{db}.exttable")

    # Data should be readable
    rows = spark.sql(f"SELECT * FROM {db}.exttable").collect()
    assert rows[0]["id"] == "1"
    assert rows[0]["value"] == "A"


# =============================================================
# write_delta special schema: contig_collection
# =============================================================

def test_write_delta_contig_collection_schema(spark):
    db = "cdmdb"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")

    schema = StructType([
        StructField("collection_id", StringType(), True),
        StructField("contig_collection_type", StringType(), True),
        StructField("ncbi_taxon_id", StringType(), True),
        StructField("gtdb_taxon_id", StringType(), True),
    ])

    df = spark.createDataFrame(
        [
            ("C1", "isolate", "NCBITaxon:123", None)
        ],
        schema=schema,
    )

    write_delta(
        spark=spark,
        df=df,
        database=db,
        table="contig_collection",
        mode="overwrite",
        data_dir=None,
    )

    result = spark.sql(f"SELECT * FROM {db}.contig_collection").collect()[0]

    assert result["collection_id"] == "C1"
    assert result["contig_collection_type"] == "isolate"
    assert result["ncbi_taxon_id"] == "NCBITaxon:123"
    assert result["gtdb_taxon_id"] is None

# =============================================================
# write_delta skip when empty
# =============================================================

def test_write_delta_empty_df(spark, capsys):
    db = "emptydb"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")

    # Create empty df
    df = spark.createDataFrame([], schema="a string")

    write_delta(
        spark=spark,
        df=df,
        database=db,
        table="emptytable",
        mode="overwrite",
        data_dir=None,
    )

    captured = capsys.readouterr().out
    assert "No data to write" in captured


# =============================================================
# preview_or_skip
# =============================================================

def test_preview_or_skip_existing(spark, capsys):
    db = "previewdb"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")

    # Drop table
    spark.sql(f"DROP TABLE IF EXISTS {db}.t1")

    # Delete physical directory to avoid LOCATION_ALREADY_EXISTS
    warehouse_dir = os.path.abspath("spark-warehouse/previewdb.db/t1")
    shutil.rmtree(warehouse_dir, ignore_errors=True)

    # Create table again
    spark.sql(f"CREATE TABLE {db}.t1 (x STRING)")

    # Insert sample row
    spark.sql(f"INSERT INTO {db}.t1 VALUES ('hello')")

    preview_or_skip(spark, db, "t1")

    captured = capsys.readouterr().out
    assert "hello" in captured


def test_preview_or_skip_missing(spark, capsys):
    db = "missingdb"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")

    preview_or_skip(spark, db, "t9999")

    out = capsys.readouterr().out
    assert "Skipping preview" in out

