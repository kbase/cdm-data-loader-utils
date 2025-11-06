import json
import click
from pyspark.sql import SparkSession, functions as F, types as T
from refseq_pipeline.core.spark_delta import build_spark, write_delta_table, cleanup_after_write
from refseq_pipeline.core.datasets_api import fetch_reports_by_taxon
from refseq_pipeline.core.cdm_parse import parse_reports


# ----------  CDM Schema ----------
CDM_SCHEMA = T.StructType([
    T.StructField("cdm_id", T.StringType(), nullable=False),
    T.StructField("taxid", T.StringType(), nullable=True),
    T.StructField("n_contigs", T.LongType(), nullable=True),
    T.StructField("contig_n50", T.LongType(), nullable=True),
    T.StructField("contig_l50", T.LongType(), nullable=True),
    T.StructField("n_scaffolds", T.LongType(), nullable=True),
    T.StructField("scaffold_n50", T.LongType(), nullable=True),
    T.StructField("scaffold_l50", T.LongType(), nullable=True),
    T.StructField("n_component_sequences", T.LongType(), nullable=True),
    T.StructField("gc_percent", T.DoubleType(), nullable=True),
    T.StructField("n_chromosomes", T.DoubleType(), nullable=True),
    T.StructField("contig_bp", T.LongType(), nullable=True),
    T.StructField("checkm_completeness", T.DoubleType(), nullable=True),
    T.StructField("checkm_contamination", T.DoubleType(), nullable=True),
    T.StructField("checkm_version", T.StringType(), nullable=True)
])


def cast_df_to_schema(df, schema: T.StructType):
    """Ensure dataframe has all columns in schema, in correct order and type."""
    for field in schema:
        if field.name not in df.columns:
            df = df.withColumn(field.name, F.lit(None).cast(field.dataType))
        else:
            df = df.withColumn(field.name, df[field.name].cast(field.dataType))
    return df.select([field.name for field in schema])


@click.command()
@click.option("--database", required=True, help="Target Spark database name.")
@click.option("--table", default="assembly_stats", show_default=True, help="Target table name.")
@click.option("--taxids-json", required=True, type=click.Path(exists=True),
              help="Path to JSON file containing list of taxids.")
@click.option("--mode", type=click.Choice(["append", "overwrite"]), default="append",
              show_default=True, help="Delta write mode.")
@click.option("--prefer-spark", is_flag=True,
              help="Prefer Spark parsing instead of pandas → Spark conversion.")


def main(database, table, taxids_json, mode, prefer_spark):
    """
    Fetch assembly reports from NCBI Datasets API by taxids and write to a Spark Delta table.
    Each entry will include its corresponding taxid to avoid NULL records.
    """
    with open(taxids_json) as f:
        taxids = json.load(f)
    if not taxids:
        print("[fetch] No taxids found. Exiting.")
        return

    spark: SparkSession = build_spark(database)
    num_success = 0
    num_failed = 0
    DATA_DIR = "/Users/yuewang/Documents/RefSeq/delta_data"

    for i, tx in enumerate(taxids, 1):
        print(f"[fetch] Processing taxid {i}/{len(taxids)}: {tx}")
        try:
            reports = list(fetch_reports_by_taxon(taxon=tx))
            if not reports:
                print(f"[fetch] taxid={tx}: no reports found.")
                continue

            df = parse_reports(reports, return_spark=prefer_spark, spark=spark)
            if df is None or (not prefer_spark and df.empty):
                print(f"[fetch] taxid={tx}: empty parsed dataframe.")
                continue

            if prefer_spark:
                sdf = df.withColumn("taxid", F.lit(str(tx)))
            else:
                df["taxid"] = str(tx)
                sdf = spark.createDataFrame(df)

            # ---------- Schema Alignment ----------
            sdf = cast_df_to_schema(sdf, CDM_SCHEMA)

            # ---------- DEBUG ----------
            print("[debug] Preview of full dataframe before writing:")
            sdf.select(*CDM_SCHEMA.fieldNames()).show(5, truncate=False)

            write_delta_table(
                sdf,
                spark,
                database,
                table,
                mode=mode,
                data_dir=DATA_DIR
            )
            num_success += 1

        except Exception as e:
            print(f"[fetch] taxid={tx} failed: {e}")
            num_failed += 1

    cleanup_after_write(spark, database, table)
    print(f"[fetch] Done. Table: {database}.{table}")
    print(f"[fetch] Summary: {num_success} succeeded, {num_failed} failed.")

    # ---------- Final dataframe ----------
    final_df = spark.read.format("delta").load(f"{DATA_DIR}/{database}/{table}")
    print("\n[verify] Example rows after write:")
    final_df.select(*CDM_SCHEMA.fieldNames()).show(10, truncate=False)
    print("[verify] TaxID NULL count =", final_df.filter(F.col("taxid").isNull()).count())


if __name__ == "__main__":
    main()



"""

python -m refseq_pipeline.cli.fetch_taxon_reports \
  --database refseq_api \
  --table assembly_stats \
  --taxids-json compare_snapshot_data/20250930_vs_20251014/changed_taxids.json \
  --mode overwrite \
  --prefer-spark
    
"""

