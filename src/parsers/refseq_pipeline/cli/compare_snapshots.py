import os
import click
from pathlib import Path
from pyspark.sql import SparkSession
from delta import configure_spark_with_delta_pip
from refseq_pipeline.core.snapshot_utils import detect_updated_or_new_hashes_from_path


def build_spark_session(app_name="Compare Snapshot Hashes") -> SparkSession:
    builder = (
        SparkSession.builder
        .appName(app_name)
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )
    return configure_spark_with_delta_pip(builder).getOrCreate()


@click.command()
@click.option("--database", required=True, help="Delta database name (e.g. refseq_api)")
@click.option("--table", required=True, help="Delta table name for hashes (e.g. assembly_hashes)")
@click.option("--old-tag", required=True, help="Snapshot tag for the previous state (e.g. 20250930)")
@click.option("--new-tag", required=True, help="Snapshot tag for the new state (e.g. 20251001)")
#@click.option("--output-dir", required=True, help="Directory to save result CSV and TSV files")
@click.option("--output-dir", default=None, help="Directory to save result CSV and TSV files")


def run_compare_snapshots(spark, delta_path: Path, old_tag: str, new_tag: str):
    """
    Core logic extracted for testing.
    Returns a Spark DataFrame with changed rows.
    """
    from refseq_pipeline.core.snapshot_utils import detect_updated_or_new_hashes_from_path

    df_all = spark.read.format("delta").load(str(delta_path))
    available_tags = [row["tag"] for row in df_all.select("tag").distinct().collect()]

    if old_tag not in available_tags:
        raise ValueError(f"[error] Old tag '{old_tag}' not found in table.")
    if new_tag not in available_tags:
        raise ValueError(f"[error] New tag '{new_tag}' not found in table.")

    return detect_updated_or_new_hashes_from_path(spark, str(delta_path), old_tag, new_tag)


def main(database, table, old_tag, new_tag, output_dir):
    """
    Compare two Delta snapshot tags and export a unified CSV of all changes
    (new, updated), and optional TSVs for specific categories.
    """
    spark = build_spark_session()

    # Resolve Delta path
    PROJECT_ROOT = Path(__file__).resolve().parents[2]
    delta_path = PROJECT_ROOT / "delta_data" / "refseq" / database / table
    print(f"[compare] Using Delta path: {delta_path}")

    if output_dir is None:
        output_dir = PROJECT_ROOT / "compare_snapshot_data" / f"{old_tag}_vs_{new_tag}"
    else:
        output_dir = Path(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    print(f"[compare] Output directory: {output_dir}")

    # Load full table and inspect available tags
    df_all = spark.read.format("delta").load(str(delta_path))
    available_tags = [row["tag"] for row in df_all.select("tag").distinct().collect()]
    print(f"[compare] Available tags: {available_tags}")

    if old_tag not in available_tags:
        raise ValueError(f"[error] Old tag '{old_tag}' not found in table. Aborting.")
    if new_tag not in available_tags:
        raise ValueError(f"[error] New tag '{new_tag}' not found in table. Aborting.")

    # Filter count for double-check logging
    old_count = df_all.filter(f"tag = '{old_tag}'").count()
    new_count = df_all.filter(f"tag = '{new_tag}'").count()
    print(f"[DEBUG] Count where tag == '{old_tag}': {old_count}")
    print(f"[DEBUG] Count where tag == '{new_tag}': {new_count}")

    # Run comparison
    df_diff = detect_updated_or_new_hashes_from_path(spark, str(delta_path), old_tag, new_tag)

    if df_diff.rdd.isEmpty():
        print("[compare] No differences found between tags.")
        return

    # Create output directory
    #os.makedirs(output_dir, exist_ok=True)
    # compare_snapshot_data
    csv_path = output_dir / "diff_summary.csv"

    # Save full diff summary as CSV
    #df_diff.coalesce(1).write.option("header", True).csv(os.path.join(output_dir, "diff_summary.csv"))
    df_diff.coalesce(1).write.option("header", True).csv(str(csv_path))
    print(f"[compare] Summary CSV written to: {csv_path}")

    # Write category-specific TSVs
    for change_type in ["new", "updated"]:
        subset = df_diff.filter(f"change_type = '{change_type}'").select("accession")
        #path = os.path.join(output_dir, f"{change_type}_accessions.tsv")
        path = output_dir / f"{change_type}_accessions.tsv"
        with open(path, "w") as f:
            for row in subset.collect():
                f.write(f"{row['accession']}\n")
        print(f"[compare] {change_type} accessions saved to: {path}")
    print(f"[compare] Summary CSV + {df_diff.count()} changes written to: {output_dir}")


if __name__ == "__main__":
    main()


"""
Usage:

python -m refseq_pipeline.cli.compare_snapshots \
  --database refseq_api \
  --table assembly_hashes \
  --old-tag 20250930 \
  --new-tag 20251014


"""


