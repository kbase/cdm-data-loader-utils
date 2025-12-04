import pytest
from pyspark.sql import Row, SparkSession
from pyspark.sql.types import (
    StructType, StructField, StringType
)

from refseq_importer.core.tables_finalize import (
    list_of_dicts_to_spark,
    finalize_tables
)


# -------------------------------------------------------------------
# Spark fixture
# -------------------------------------------------------------------
@pytest.fixture(scope="session")
def spark():
    spark = (
        SparkSession.builder
        .master("local[*]")
        .appName("test-tables-finalize")
        .getOrCreate()
    )
    yield spark
    spark.stop()


# -------------------------------------------------------------------
# Test list_of_dicts_to_spark
# -------------------------------------------------------------------
def test_list_of_dicts_to_spark(spark):
    schema = StructType([
        StructField("a", StringType(), True),
        StructField("b", StringType(), True),
    ])

    rows = [{"a": "1", "b": "x"}, {"a": "2", "b": "y"}]
    df = list_of_dicts_to_spark(spark, rows, schema)

    out = {(r.a, r.b) for r in df.collect()}
    assert out == {("1", "x"), ("2", "y")}


# -------------------------------------------------------------------
# Test finalize_tables end-to-end
# -------------------------------------------------------------------
def test_finalize_tables_basic(spark):

    # ---------- entity ----------
    e_schema = StructType([
        StructField("entity_id", StringType(), True),
        StructField("entity_type", StringType(), True),
        StructField("data_source", StringType(), True),
        StructField("created", StringType(), True),
        StructField("updated", StringType(), True),
    ])

    e1 = spark.createDataFrame(
        [Row(entity_id="E1", entity_type="genome",
             data_source="RefSeq", created="2020", updated="2021")],
        schema=e_schema
    )
    e2 = spark.createDataFrame(
        [Row(entity_id="E2", entity_type="genome",
             data_source="RefSeq", created="2020", updated="2021")],
        schema=e_schema
    )

    # ---------- contig_collection (schema REQUIRED due to None!) ----------
    coll_schema = StructType([
        StructField("collection_id", StringType(), True),
        StructField("contig_collection_type", StringType(), True),
        StructField("ncbi_taxon_id", StringType(), True),
        StructField("gtdb_taxon_id", StringType(), True),
    ])

    c1 = spark.createDataFrame(
        [
            Row(
                collection_id="E1",
                contig_collection_type="isolate",
                ncbi_taxon_id="NCBITaxon:1",
                gtdb_taxon_id=None,
            )
        ],
        schema=coll_schema
    )

    c2 = spark.createDataFrame(
        [
            Row(
                collection_id="E2",
                contig_collection_type="isolate",
                ncbi_taxon_id="NCBITaxon:2",
                gtdb_taxon_id=None,
            )
        ],
        schema=coll_schema
    )

    # ---------- name ----------
    names = [
        {"entity_id": "E1", "name": "A", "description": "d1", "source": "RefSeq"},
        {"entity_id": "E2", "name": "B", "description": "d2", "source": "RefSeq"},
    ]

    # ---------- identifier ----------
    identifiers = [
        {"entity_id": "E1", "identifier": "BioSample:1", "source": "RefSeq", "description": "bs"},
        {"entity_id": "E2", "identifier": "BioSample:2", "source": "RefSeq", "description": "bs"},
    ]

    df_entity, df_coll, df_name, df_ident = finalize_tables(
        spark,
        [e1, e2],
        [c1, c2],
        names,
        identifiers
    )

    # ---------- Assertions ----------
    assert df_entity.count() == 2
    assert df_coll.count() == 2
    assert df_name.count() == 2
    assert df_ident.count() == 2

    assert set(r.entity_id for r in df_entity.collect()) == {"E1", "E2"}
    assert set(r.collection_id for r in df_coll.collect()) == {"E1", "E2"}
    assert set(r.name for r in df_name.collect()) == {"A", "B"}
    assert set(r.identifier for r in df_ident.collect()) == {"BioSample:1", "BioSample:2"}

    