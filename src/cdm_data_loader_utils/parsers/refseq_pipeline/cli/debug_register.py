from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from refseq_pipeline.core.spark_delta import read_delta_table, register_table

"""
PYTHONPATH=/global_share/alinawang/cdm-data-loader-utils/src/parsers python debug_register.py
"""

database = "refseq_api"
table = "assembly_hashes"
delta_path = "/global_share/alinawang/cdm-data-loader-utils/delta_data/refseq/refseq_api/assembly_hashes"

builder = (
    SparkSession.builder.appName("Delta Table Inspector")
    .master("local[*]")
    .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
    .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
)
spark = configure_spark_with_delta_pip(builder).getOrCreate()


register_table(spark, database, table, delta_path)
df = read_delta_table(spark, database, table)

print("[SCHEMA]")
df.printSchema()

print("\n[Show dataset]")
df.show(10, truncate=False)

print("\n[count]")
print(df.count())

if "updated" in df.columns:
    print("\n[timestamp range]")
    df.selectExpr("min(updated)", "max(updated)").show()
