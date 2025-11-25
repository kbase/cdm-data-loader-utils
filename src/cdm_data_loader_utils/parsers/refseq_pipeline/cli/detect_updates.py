import os
import click
from pathlib import Path
from refseq_pipeline.core.spark_delta import build_spark
from refseq_pipeline.core.snapshot_utils import detect_updated_or_new_hashes_from_path

"""
We implemented a snapshot comparison tool for RefSeq assemblies. 
By comparing SHA256 hashes from two Delta table tags 
We identified several genomes that were either updated or newly added. 
This diff was saved as a CSV file containing accession IDs, change types, and old/new checksums for downstream processing or API calls.


Usage: 

  PYTHONPATH=src/parsers \
  python src/parsers/refseq_pipeline/cli/detect_updates.py \
  --database refseq_api \
  --table assembly_hashes \
  --old-tag 20251014 \
  --new-tag 20251120 \
  --output compare_snapshot_data

"""


@click.command()
@click.option("--database", required=True, help="Spark SQL database name.")
@click.option("--table", default="assembly_hashes", show_default=True, help="Delta table name.")
@click.option("--old-tag", required=True, help="Previous snapshot tag to compare.")
@click.option("--new-tag", required=True, help="Current snapshot tag to compare.")
@click.option("--output", default=None, help="CSV output directory.")


def main(database, table, old_tag, new_tag, output):
    """
    Compare two Delta snapshot tags (ex: 20200930 VS 20251001)
    Detect genomes that are new or updated between the two versions
    """

    print(f"\n[detect] Comparing snapshot tags: {old_tag} → {new_tag}")

    # --- Initialize Spark ---
    spark = build_spark(database)

    # --- Preview available tags ---
    try:
        print(f"[detect] Available tags in {database}.{table}:")
        spark.sql(f"SELECT DISTINCT tag FROM {database}.{table} ORDER BY tag").show(truncate=False)
    except Exception as e:
        print(f"[WARN] Could not preview registered table tags: {e}")


    # --- Build Delta path ---
    #project_root = Path(__file__).resolve().parents[2]  
    #delta_path = project_root / "delta_data" / "refseq" / database / table

    project_root = Path("/global_share/alinawang/cdm-data-loader-utils")
    delta_path = project_root / "delta_data" / "refseq" / "refseq_api" / "assembly_hashes"
    
    print(f"[detect] Using Delta path: {delta_path}")

    if not delta_path.exists():
        print(f"[ERROR] Delta path not found: {delta_path}")
        raise SystemExit(1)
    
    # --- Debugging block to inspect Delta contents ---
    try:
        print("\n[debug] Loading Delta table from path...")
        df_all = spark.read.format("delta").load(str(delta_path))

        print("[debug] Schema:")
        df_all.printSchema()

        print("[debug] Columns:")
        print(df_all.columns)

        print("[debug] Distinct tags (up to 100):")
        df_all.select("tag").distinct().show(100, truncate=False)

        print(f"[debug] Count where tag == '{new_tag}':")
        print(df_all.filter(df_all["tag"] == new_tag).count())

        print(f"[debug] Count where tag == '{old_tag}':")
        print(df_all.filter(df_all["tag"] == old_tag).count())
    except Exception as e:
        print(f"[WARN] Debug block failed: {e}")


    # --- Run comparison ---
    df = detect_updated_or_new_hashes_from_path(spark, str(delta_path), old_tag, new_tag)

    if df is None or df.rdd.isEmpty():
        print(f"[detect] No updated or new genomes found between {old_tag} and {new_tag}.")
        return

    count = df.count()
    print(f"[detect] Found {count} updated/new genomes between {old_tag} and {new_tag}.")
    df.show(20, truncate=False)

    # --- Output to CSV ---
    if output:
        os.makedirs(output, exist_ok=True)
        df.coalesce(1).write.option("header", True).csv(output, mode="overwrite")
        print(f"[detect] CSV diff written to folder: {output}\n")
    else:
        print("[detect] No output path specified; skipping CSV write.\n")


if __name__ == "__main__":
    main()


