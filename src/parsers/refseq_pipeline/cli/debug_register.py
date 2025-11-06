import os
from pyspark.sql import SparkSession
from delta import configure_spark_with_delta_pip
from refseq_pipeline.core.spark_delta import register_table, read_delta_table


###### Test script to inspect a Delta table ######

def main():
    # set parameters
    database = "refseq_api"
    table = "assembly_hashes"
    delta_path = os.path.abspath("delta_data/refseq/refseq_api/assembly_hashes")

    # ===== Spark Initialization =====
    builder = (
        SparkSession.builder
        .appName("Delta Table Inspector")
        .master("local[*]")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )

    spark = configure_spark_with_delta_pip(builder).getOrCreate()

    # ===== Register Delta Table =====
    register_table(spark, database, table, delta_path)

    # ===== Read Delta Table =====
    df = read_delta_table(spark, database, table)

    print("\n[SCHEMA]")
    df.printSchema()

    print("\n[Show dataset]")
    df.show(10, truncate=False)

    print("\n[count]")
    print(df.count())

    if "updated" in df.columns:
        print("\n[timestamp range]")
        df.selectExpr("min(updated)", "max(updated)").show()

    spark.stop()


if __name__ == "__main__":
    main()


"""

python -m refseq_pipeline.cli.debug_register

OUPTUT 

[SCHEMA]
root
 |-- accession: string (nullable = true)
 |-- ftp_path: string (nullable = true)
 |-- kind: string (nullable = true)
 |-- content_sha256: string (nullable = true)
 |-- retrieved_at: string (nullable = true)
 |-- tag: string (nullable = true)

 
+---------------+----------------------------------------------------------------------------------------------------+----+----------------------------------------------------------------+-------------------------+--------+
|accession      |ftp_path                                                                                            |kind|content_sha256                                                  |retrieved_at             |tag     |
+---------------+----------------------------------------------------------------------------------------------------+----+----------------------------------------------------------------+-------------------------+--------+
|GCF_000158355.2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/355/GCF_000158355.2_Citrobacter_sp_30_2_V2     |md5 |98e87f1cdc11052fe4f2d0d196c9e50bd783ca4772ce3addf44acf1b796abea8|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158375.2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/375/GCF_000158375.2_Clostridium_sp_7_2_43FAA_V2|md5 |d932c7ac049ad685fd444fccd5d24b450dd0147eef8c78ade1477a1a766146d4|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158395.1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/395/GCF_000158395.1_ASM15839v1                 |md5 |70cbd156a704d4453149d1c8f7e55074c1bd0340e0bfad887cb929467ec44acc|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158415.2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/415/GCF_000158415.2_Escherichia_sp_4_1_40B_V2  |md5 |9c70abfac82c8856f5a1a2b574f0e8b7a2f598ae9bf4888998d341cc29616499|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158435.2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/435/GCF_000158435.2_Heli_bilis_ATCC43879_V2    |md5 |826e682e6d2fc643ee70e667f62b9a9884d9a82482ce2777ad29b50d393d3dbe|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158455.1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/455/GCF_000158455.1_ASM15845v1                 |md5 |46624e180f51c2862cc8798372753146355620b30815b3663c757147ff1634ce|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158475.2|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/475/GCF_000158475.2_Oxal_for_HOxBLS_2_V2       |md5 |3b5ebd60bf348fd62e72b725f8f9c4abf1ac7485eebd4935ff8b55d8c9559bb5|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158495.1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/495/GCF_000158495.1_ASM15849v1                 |md5 |0a173b7e5469ab2734c8873c0af706a48b33d8d37cf900800bb16a73618242c1|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158575.1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/575/GCF_000158575.1_ASM15857v1                 |md5 |3e3d2ec8e45aea715b1fb6895dcadba5980ee90477c1c7e1f4c5f91d28f923d5|2025-10-14T21:13:17+00:00|20251014|
|GCF_000158595.1|https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/158/595/GCF_000158595.1_ASM15859v1                 |md5 |f1ef45b1527c7d8d0369e0d44f9f96a6ec92dcb3f162fc731e159e34e1317108|2025-10-14T21:13:17+00:00|20251014|
+---------------+----------------------------------------------------------------------------------------------------+----+----------------------------------------------------------------+-------------------------+--------+

[count]
940431

"""

