import os
import shutil
from pathlib import Path
from typing import Iterable
from delta.tables import DeltaTable
from pyspark.sql import SparkSession, functions as F, Window


def build_spark(database: str, master: str = "local[*]") -> SparkSession:
    """
    Initialize SparkSession with Delta Lake support and create the database if not exists.
    """
    from delta import configure_spark_with_delta_pip

    builder = (
        SparkSession.builder.appName("RefSeq Pipeline")
        .master(master)
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )

    spark = configure_spark_with_delta_pip(builder).getOrCreate()
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {database}")
    return spark


def write_delta_table(
    sdf, spark: SparkSession, database: str, table: str, mode: str = "append", data_dir: str | None = None
):
    """
    Write a Spark DataFrame to either external or managed Delta table.
    Automatically cleans up non-Delta directories.
    """
    if sdf.rdd.isEmpty():
        print(f"[delta] empty dataframe, skip write {database}.{table}")
        return

    # ===== External Delta Table Mode =====
    if data_dir:
        db = database.replace(".db", "")
        target_path = os.path.abspath(os.path.join(data_dir, db, table))
        Path(target_path).mkdir(parents=True, exist_ok=True)

        is_delta = DeltaTable.isDeltaTable(spark, target_path)
        spark.sql(f"CREATE DATABASE IF NOT EXISTS `{db}`")
        full_table = f"`{db}`.`{table}`"

        try:
            if not is_delta:
                print(f"[delta] Target path {target_path} is not Delta table. Cleaning before write...")
                shutil.rmtree(target_path, ignore_errors=True)
                Path(target_path).mkdir(parents=True, exist_ok=True)
                mode = "overwrite"

            writer = sdf.write.format("delta").mode(mode)
            if mode == "append":
                writer = writer.option("mergeSchema", "true")
            else:
                writer = writer.option("overwriteSchema", "true")

            print(f"[delta] DataFrame schema: {sdf.schema.simpleString()}")
            print(f"[delta] Row count to write: {sdf.count()}")

            writer.save(target_path)
            spark.sql(f"CREATE TABLE IF NOT EXISTS {full_table} USING DELTA LOCATION '{target_path}'")
            print(f"[delta] saved {database}.{table} to {target_path} (mode={mode})")

        except Exception as e:
            print(f"[delta] ERROR while writing {full_table}: {e}")
            raise

    # ===== Managed Delta Table Mode =====
    else:
        full = f"{database}.{table}"
        try:
            writer = sdf.write.format("delta").mode(mode)

            if mode == "append":
                writer = writer.option("mergeSchema", "true")
            else:
                writer = writer.option("overwriteSchema", "true")

            print(f"[delta] DataFrame schema: {sdf.schema.simpleString()}")
            print(f"[delta] Row count to write: {sdf.count()}")

            writer.saveAsTable(full)
            print(f"[delta] saved {full} (mode={mode})")
        except Exception as e:
            print(f"[delta] ERROR while writing managed table {full}: {e}")
            raise


def cleanup_after_write(
    spark: SparkSession,
    database: str,
    table: str,
    *,
    bad_checkm_patterns: Iterable[str] = ("net.razorvine.pickle.objects.%",),
    prefer_cols=("checkm_version",),
    keep_latest_by=None,
    do_optimize=False,
    do_vacuum=False,
    vacuum_hours=168,
    zorder_by=("cdm_id",),
):
    """
    Clean and deduplicate a Delta table by CDM ID after write.
    Optionally OPTIMIZE and VACUUM the table.
    """
    full_sql = f"`{database}`.`{table}`"
    full = f"{database}.{table}"

    if not spark.catalog.tableExists(full):
        print(f"[cleanup] table {full} not found; skip.")
        return

    src = spark.table(full_sql)
    if src.rdd.isEmpty():
        print(f"[cleanup] table {full} is empty; skip.")
        return

    # --- Clean invalid checkm_version patterns ---
    cv = F.col("checkm_version")
    for pat in bad_checkm_patterns:
        src = src.withColumn("checkm_version", F.when(cv.like(pat), F.lit(None)).otherwise(F.col("checkm_version")))

    # --- Deduplication rule ---
    weight_exprs = [F.when(F.col(c).isNull(), 1).otherwise(0) for c in prefer_cols]
    total_weight = sum(weight_exprs)

    order_cols = [total_weight.asc()]
    if keep_latest_by and keep_latest_by in src.columns:
        order_cols.append(F.col(keep_latest_by).desc())
    order_cols.append(F.col("cdm_id").asc())

    w = Window.partitionBy("cdm_id").orderBy(*order_cols)
    cleaned = src.withColumn("rn", F.row_number().over(w)).filter("rn = 1").drop("rn")

    # --- Overwrite the original table location ---
    (cleaned.write.format("delta").mode("overwrite").option("overwriteSchema", "true").saveAsTable(full))

    # --- OPTIMIZE ---
    if do_optimize:
        try:
            z = ", ".join([f"`{c}`" for c in zorder_by]) if zorder_by else ""
            spark.sql(f"OPTIMIZE {full_sql} {'ZORDER BY (' + z + ')' if z else ''}")
        except Exception as e:
            print("OPTIMIZE skipped:", e)

    # --- VACUUM ---
    if do_vacuum:
        try:
            spark.sql(f"VACUUM {full_sql} RETAIN {vacuum_hours} HOURS")
        except Exception as e:
            print("VACUUM skipped:", e)


def register_table(spark: SparkSession, database: str, table: str, path: str):
    """
    Explicitly register an external Delta table from a path.
    This is useful when comparing snapshots or accessing by SQL.
    """
    spark.sql(f"CREATE DATABASE IF NOT EXISTS `{database}`")
    full_table = f"`{database}`.`{table}`"
    if not spark.catalog.tableExists(full_table):
        spark.sql(f"""
            CREATE TABLE {full_table}
            USING DELTA
            LOCATION '{os.path.abspath(path)}'
        """)
        print(f"Registered table {full_table} at {path}")
    else:
        print(f"Table {full_table} already registered.")


def read_delta_table(spark: SparkSession, database: str, table: str):
    """
    Load an existing Delta table from metastore.
    """
    full = f"{database}.{table}"
    if spark.catalog.tableExists(full):
        return spark.table(full)
    else:
        raise ValueError(f"Table {full} not found in catalog.")
