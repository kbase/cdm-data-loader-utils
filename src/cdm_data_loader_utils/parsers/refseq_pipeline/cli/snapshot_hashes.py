import time
import click
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from pyspark.sql import functions as F
from refseq_pipeline.core.refseq_io import parse_assembly_summary
from refseq_pipeline.core.spark_delta import build_spark, write_delta_table
from refseq_pipeline.core.hashes_snapshot import snapshot_hashes_for_accessions


# ---------------------------
# Helper functions
# ---------------------------
def chunk_list(lst, size):
    """Split a list into evenly sized chunks."""
    for i in range(0, len(lst), size):
        yield lst[i:i + size]


def chunked_snapshot(chunk, acc_index, kind, spark, tag, fast_mode):
    """Process one chunk of accessions to generate hash DataFrame."""
    try:
        sdf = snapshot_hashes_for_accessions(
            chunk, acc_index=acc_index, kind=kind,
            spark=spark, fast_mode=fast_mode
        )
        if sdf is None or sdf.rdd.isEmpty():
            return None
        return sdf.withColumn("tag", F.lit(tag))
    except Exception as e:
        print(f"[snapshot-thread] Error in chunk: {e}")
        return None


# ---------------------------
# CLI entrypoint
# ---------------------------
@click.command()
@click.option("--database", required=True, help="Spark SQL database to write into.")
@click.option("--hash-table", default="assembly_hashes", show_default=True,
              help="Delta table name to write.")
@click.option("--tag", required=True, help="Snapshot tag (e.g., release version or date).")
@click.option("--kind", type=click.Choice(["annotation", "md5"]), default="annotation", show_default=True,
              help="Hashing strategy to use.")
@click.option("--index-path", default=None,
              help="Path to local assembly_summary_refseq.tsv. If not provided, auto-detect latest in bronze/refseq/indexes/")
@click.option("--chunk-size", default=1000, show_default=True,
              help="Number of accessions to process per chunk.")
@click.option("--max-workers", default=8, show_default=True,
              help="Number of parallel threads.")
@click.option("--fast-mode/--no-fast-mode", default=True, show_default=True,
              help="Use precomputed FTP hashes when available.")
@click.option("--data-dir", default=None, help="Optional path to store Delta table files.")


def main(database, hash_table, tag, kind, index_path, chunk_size, max_workers, fast_mode, data_dir):
    """Main entry for taking hash snapshot of RefSeq assemblies."""

    start = time.time()
    
    print(f"[snapshot] Target: {database}.{hash_table} | Tag: {tag}")
    print(f"[snapshot] Kind: {kind} | Chunk size: {chunk_size} | Max workers: {max_workers}")
    print(f"[snapshot] Fast mode: {fast_mode}")

    # --- Build Spark session ---
    spark = build_spark(database)

    # --- Load RefSeq TSV index ---
    if index_path:
        index_path = Path(index_path)
        print(f"[snapshot] Loading index from: {index_path}")
    else:
        bronze_dir = Path("bronze/refseq/indexes")
        all_files = sorted(bronze_dir.glob("assembly_summary_refseq.*.tsv"))
        if not all_files:
            raise FileNotFoundError("No local index TSV files found in 'bronze/refseq/indexes'")
        index_path = all_files[-1]
        print(f"[snapshot] Auto-selected latest index: {index_path}")

    # Read file content and parse
    txt = Path(index_path).read_text(encoding="utf-8")
    acc_index = parse_assembly_summary(txt)
    accessions = list(acc_index.keys())
    print(f"[snapshot] Parsed {len(accessions)} accessions from index.")

    # --- Split into chunks --- 
    chunks = list(chunk_list(accessions, chunk_size))
    print(f"[snapshot] Processing {len(chunks)} chunks...")

    results = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for result in executor.map(
            lambda c: chunked_snapshot(c, acc_index, kind, spark, tag, fast_mode),
            chunks
        ):
            if result:
                results.append(result)

    # --- Combine and write ---
    if results:
        final_df = results[0]
        for sdf in results[1:]:
            final_df = final_df.unionByName(sdf)

        write_delta_table(final_df, spark, database, hash_table, mode="append", data_dir=data_dir)
        print(f"[snapshot] Snapshot written to {database}.{hash_table} with tag '{tag}'.")
    else:
        print("[snapshot] No results to write. Exiting.")

    print(f"[snapshot] Completed in {round((time.time() - start) / 60, 2)} minutes.")


if __name__ == "__main__":
    main()



"""

PYTHONPATH=src/parsers python src/parsers/refseq_pipeline/cli/snapshot_hashes.py \
  --database refseq_api \
  --hash-table assembly_hashes \
  --tag 20251014 \
  --kind md5 \
  --index-path bronze/refseq/indexes/assembly_summary_refseq.20250930.tsv \
  --chunk-size 500 \
  --fast-mode \
  --max-workers 12 \
  --data-dir delta_data/refseq


"""


