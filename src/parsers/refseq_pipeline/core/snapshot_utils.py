from pyspark.sql import SparkSession, functions as F, DataFrame


def get_changed_accessions(spark: SparkSession, delta_path: str, old_tag: str, new_tag: str) -> list[str]:
    df = spark.read.format("delta").load(delta_path)
    df_old = df.filter(F.col("tag") == old_tag).select("accession", "content_sha256")
    df_new = df.filter(F.col("tag") == new_tag).select("accession", "content_sha256")

    joined = df_old.alias("a").join(df_new.alias("b"), on="accession", how="inner")
    changed = joined.filter(F.col("a.content_sha256") != F.col("b.content_sha256"))

    return [row["accession"] for row in changed.select("accession").collect()]


def get_new_accessions(spark: SparkSession, delta_path: str, old_tag: str, new_tag: str) -> list[str]:
    df = spark.read.format("delta").load(delta_path)
    df_old = df.filter(F.col("tag") == old_tag).select("accession").distinct()
    df_new = df.filter(F.col("tag") == new_tag).select("accession").distinct()

    new_only = df_new.join(df_old, on="accession", how="left_anti")
    return [row["accession"] for row in new_only.select("accession").collect()]


def get_removed_accessions(spark: SparkSession, delta_path: str, old_tag: str, new_tag: str) -> list[str]:
    df = spark.read.format("delta").load(delta_path)
    df_old = df.filter(F.col("tag") == old_tag).select("accession").distinct()
    df_new = df.filter(F.col("tag") == new_tag).select("accession").distinct()

    removed_only = df_old.join(df_new, on="accession", how="left_anti")
    return [row["accession"] for row in removed_only.select("accession").collect()]

# ========== For CLI detect_updates.py ==========

def detect_updated_or_new_hashes_from_path(
    spark: SparkSession,
    delta_path: str,
    old_tag: str,
    new_tag: str) -> DataFrame:
    """
    Return a unified DataFrame of changes between two snapshots (new or updated accessions),
    using Delta table path instead of database.table.
    """
    df = spark.read.format("delta").load(delta_path).select("accession", "content_sha256", "tag")

    df_old = df.filter(F.col("tag") == old_tag).selectExpr(
        "accession as acc", "content_sha256 as old_sha256"
    )
    df_new = df.filter(F.col("tag") == new_tag).selectExpr(
        "accession as acc", "content_sha256 as new_sha256"
    )

    joined = df_new.join(df_old, on="acc", how="outer")

    return (
        joined.withColumn(
            "change_type",
            F.when(F.col("old_sha256").isNull(), F.lit("new"))
             .when(F.col("new_sha256").isNull(), F.lit("deleted"))
             .when(F.col("old_sha256") != F.col("new_sha256"), F.lit("updated"))
             .otherwise(F.lit("unchanged"))
        )
        .filter("change_type in ('new', 'updated')")
        .selectExpr("acc as accession", "change_type", "old_sha256", "new_sha256")
    )

