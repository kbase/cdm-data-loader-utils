"""

PYTHONPATH=src/parsers python src/parsers/refseq_pipeline/cli/compare_snapshots.py \
  --database refseq_api \
  --table assembly_hashes \
  --old-tag 20250930 \
  --new-tag 20251120 \
  --output-dir compare_snapshot_data/20250930_vs_20251120

"""

import os
import click
from pathlib import Path
from pyspark.sql import SparkSession
from delta import configure_spark_with_delta_pip
from refseq_pipeline.core.snapshot_utils import detect_updated_or_new_hashes_from_path


def build_spark_session(app_name="Compare Snapshot Hashes") -> SparkSession:
    builder = (
        SparkSession.builder.appName(app_name)
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )
    return configure_spark_with_delta_pip(builder).getOrCreate()


def run_compare_snapshots(spark, delta_path: Path, old_tag: str, new_tag: str):
    """
    Return a DataFrame with new or updated records.
    """
    df_all = spark.read.format("delta").load(str(delta_path))
    available_tags = [row["tag"] for row in df_all.select("tag").distinct().collect()]
    print("[Available Tags]")
    for t in available_tags:
        print(" -", t)

    if old_tag not in available_tags or new_tag not in available_tags:
        raise ValueError(f"Tags {old_tag} or {new_tag} not found in Delta table.")

    # Tag count info
    for tag in [old_tag, new_tag]:
        count = df_all.filter(f"tag = '{tag}'").count()
        print(f"[INFO] Count of entries with tag '{tag}': {count}")

    # sha256 preview
    df_all.select("tag", "accession", "content_sha256").orderBy("tag").show(20, truncate=False)

    # Null check
    df_nulls = df_all.filter("content_sha256 IS NULL")
    null_count = df_nulls.count()
    if null_count > 0:
        print(f"[WARNING] Found {null_count} rows with null sha256!")
        df_nulls.select("tag", "accession").show()

    # Run detection logic
    df_diff = detect_updated_or_new_hashes_from_path(spark, str(delta_path), old_tag, new_tag)

    print(f"[compare] Diff result count: {df_diff.count()}")
    df_diff.show(20, truncate=False)

    return df_diff


@click.command()
@click.option("--database", required=True, help="Delta database name (refseq_api)")
@click.option("--table", required=True, help="Delta table name for hashes (assembly_hashes)")
@click.option("--old-tag", required=True, help="Snapshot tag for the previous state (20250930)")
@click.option("--new-tag", required=True, help="Snapshot tag for the new state (20251001)")
@click.option("--output-dir", default=None, help="Directory to save result CSV and TSV files")
def main(database, table, old_tag, new_tag, output_dir):
    spark = build_spark_session()

    project_root = Path("/global_share/alinawang/cdm-data-loader-utils")
    delta_path = project_root / "delta_data" / "refseq" / "refseq_api" / "assembly_hashes"

    print(f"[compare] Using Delta path: {delta_path}")
    if not delta_path.exists():
        raise FileNotFoundError(f"[ERROR] Delta path not found: {delta_path}")

    if output_dir is None:
        output_dir = PROJECT_ROOT / "compare_snapshot_data" / f"{old_tag}_vs_{new_tag}"
    else:
        output_dir = Path(output_dir)

    os.makedirs(output_dir, exist_ok=True)
    print(f"[compare] Output will be saved to: {output_dir}")

    try:
        df_diff = run_compare_snapshots(spark, delta_path, old_tag, new_tag)
    except Exception as e:
        print("[ERROR] Failed during snapshot comparison:", e)
        raise

    if df_diff.rdd.isEmpty():
        print("[compare] No differences detected between snapshots.")
        return

    # Write diff summary CSV
    summary_csv = output_dir / "diff_summary.csv"
    df_diff.coalesce(1).write.option("header", True).mode("overwrite").csv(str(summary_csv))
    print(f"[compare] Summary CSV written: {summary_csv}")

    # Write new/updated accessions as TSV
    for change_type in ["new", "updated"]:
        df_subset = df_diff.filter(f"change_type = '{change_type}'").select("accession")
        tsv_path = output_dir / f"{change_type}_accessions.tsv"
        with open(tsv_path, "w") as f:
            for row in df_subset.collect():
                f.write(f"{row['accession']}\n")
        print(f"[compare] {change_type} accession list saved: {tsv_path}")

    print(f"[compare] Done. {df_diff.count()} total differences written.")


if __name__ == "__main__":
    main()
