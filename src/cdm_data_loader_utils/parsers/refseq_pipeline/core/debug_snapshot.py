from pyspark.sql import functions as F

from cdm_data_loader_utils.parsers.refseq_pipeline.core.hashes_snapshot import (
    snapshot_hashes_for_accessions,
    write_hash_snapshot,
)
from cdm_data_loader_utils.parsers.refseq_pipeline.core.refseq_io import load_refseq_assembly_index
from cdm_data_loader_utils.parsers.refseq_pipeline.core.spark_delta import build_spark

## python -m refseq_pipeline.core.debug_snapshot ##

# --- Setup ---
spark = build_spark("refseq_debug")
acc_index = load_refseq_assembly_index()

# Sample accessions to test
sample_accs = [
    "GCF_000001405.40",  # human
    "GCF_000001635.27",  # mouse
]
tag = "debug_test"

# --- Run snapshot hash fetch ---
sdf = snapshot_hashes_for_accessions(sample_accs, acc_index=acc_index, kind="auto", spark=spark, fast_mode=False)

# Add tag column
sdf = sdf.withColumn("tag", F.lit(tag))

# --- Show preview ---
print("=== Preview ===")
sdf.show(truncate=False)

# --- Write to Delta (optional) ---
write_hash_snapshot(
    spark=spark, sdf=sdf, database="refseq_debug", table="assembly_hashes_test", data_dir="output_delta"
)

# --- Query from Delta to verify ---
print("\n=== Verifying from Delta Table ===")
spark.sql("USE refseq_debug")
spark.sql(f"""
    SELECT accession, kind, tag, content_sha256, retrieved_at
    FROM assembly_hashes_test
    WHERE tag = '{tag}'
""").show(truncate=False)
