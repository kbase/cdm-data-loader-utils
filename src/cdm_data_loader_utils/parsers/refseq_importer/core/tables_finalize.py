from typing import List, Dict, Any
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql import Row
from pyspark.sql.types import StructType, StructField, StringType
from .spark_delta import write_delta, preview_or_skip


# -------------------------------------------------------------------
# Convert list[dict] → Spark DataFrame with explicit schema
# -------------------------------------------------------------------
def list_of_dicts_to_spark(
    spark: SparkSession,
    rows: List[Dict[str, Any]],
    schema: StructType
) -> DataFrame:
    if not rows:
        return spark.createDataFrame([], schema=schema)
    spark_rows = [Row(**r) for r in rows]
    return spark.createDataFrame(spark_rows, schema=schema)


# -------------------------------------------------------------------
# finalize_tables: convert 4 partial lists into Spark DFs
# -------------------------------------------------------------------
def finalize_tables(
    spark: SparkSession,
    entities: List[DataFrame],
    collections: List[DataFrame],
    names: List[Dict[str, Any]],
    identifiers: List[Dict[str, Any]]
):
    """
    Combine Spark partial entity/collection DFs,
    and convert name + identifier lists into Spark DFs.
    """

    # entity table
    if entities:
        df_entity = entities[0]
        for df in entities[1:]:
            df_entity = df_entity.unionByName(df, allowMissingColumns=True)
    else:
        df_entity = spark.createDataFrame([], schema="""
            entity_id STRING,
            entity_type STRING,
            data_source STRING,
            created STRING,
            updated STRING
        """)

    # contig_collection table
    if collections:
        df_coll = collections[0]
        for df in collections[1:]:
            df_coll = df_coll.unionByName(df, allowMissingColumns=True)
    else:
        df_coll = spark.createDataFrame([], schema="""
            collection_id STRING,
            contig_collection_type STRING,
            ncbi_taxon_id STRING,
            gtdb_taxon_id STRING
        """)

    # name table
    name_schema = StructType([
        StructField("entity_id", StringType(), False),
        StructField("name", StringType(), False),
        StructField("description", StringType(), True),
        StructField("source", StringType(), True),
    ])
    df_name = list_of_dicts_to_spark(spark, names, name_schema)

    # identifier table
    ident_schema = StructType([
        StructField("entity_id", StringType(), False),
        StructField("identifier", StringType(), False),
        StructField("source", StringType(), True),
        StructField("description", StringType(), True),
    ])
    df_ident = list_of_dicts_to_spark(spark, identifiers, ident_schema)

    return df_entity, df_coll, df_name, df_ident



# -------------------------------------------------------------------
# Write all 4 tables 
# -------------------------------------------------------------------
def write_and_preview(
    spark: SparkSession,
    database: str,
    mode: str,
    df_entity: DataFrame,
    df_coll: DataFrame,
    df_name: DataFrame,
    df_ident: DataFrame,
    data_dir: str | None = None,
):
    write_delta(spark, df_entity, database, "entity", mode, data_dir=data_dir)
    write_delta(spark, df_coll, database, "contig_collection", mode, data_dir=data_dir)
    write_delta(spark, df_name, database, "name", mode, data_dir=data_dir)
    write_delta(spark, df_ident, database, "identifier", mode, data_dir=data_dir)

    print("\nDelta tables written:")
    for tbl in ["datasource", "entity", "contig_collection", "name", "identifier"]:
        preview_or_skip(spark, database, tbl)

        