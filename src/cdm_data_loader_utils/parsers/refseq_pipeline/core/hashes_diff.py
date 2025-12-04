import logging
from typing import Optional, Dict, List
from pyspark.sql import SparkSession, DataFrame, functions as F

# Initialize logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")


def diff_hash_and_get_changed_taxids(
    spark: SparkSession,
    database: str,
    hash_table: str,
    acc_index: Dict[str, Dict[str, str]],
    tag_new: Optional[str] = None,
    tag_old: Optional[str] = None,
    debug: bool = False,
) -> List[str]:
    """
    Compare hash snapshots and return taxon IDs whose assemblies have changed.

    Args:
        spark (SparkSession): Active Spark session.
        database (str): Database name (Delta).
        hash_table (str): Table name storing hash snapshots.
        acc_index (dict): Mapping of accession → metadata (with taxid or species_taxid).
        tag_new (str, optional): Tag for newer snapshot.
        tag_old (str, optional): Tag for older snapshot.
        debug (bool): Whether to print changed accessions and taxids.

    Returns:
        List[str]: Unique taxon IDs affected by hash content changes.
    """
    full_table = f"{database}.{hash_table}"
    if not spark.catalog.tableExists(full_table):
        logger.warning(f"[diff] Table {full_table} not found.")
        return []

    hash_df: DataFrame = spark.table(full_table)
    if hash_df.rdd.isEmpty():
        logger.warning(f"[diff] Table {full_table} is empty.")
        return []

    # --- Tag or timestamp selection ---
    if tag_new and tag_old:
        cond_new = F.col("tag") == tag_new
        cond_old = F.col("tag") == tag_old
        logger.info(f"[diff] Using tag comparison: old='{tag_old}', new='{tag_new}'")
    else:
        ts = [r["retrieved_at"] for r in hash_df.select("retrieved_at").distinct().orderBy("retrieved_at").tail(2)]
        if len(ts) < 2:
            logger.warning(f"[diff] Only one snapshot available: {ts}")
            return []
        old_ts, new_ts = ts[0], ts[1]
        cond_new = F.col("retrieved_at") == new_ts
        cond_old = F.col("retrieved_at") == old_ts
        logger.info(f"[diff] Using timestamp comparison: old='{old_ts}', new='{new_ts}'")

    # --- Extract old/new hashes ---
    old_df = (
        hash_df.filter(cond_old)
        .select("accession", "kind", "content_sha256")
        .withColumnRenamed("content_sha256", "old_sha")
    )
    new_df = (
        hash_df.filter(cond_new)
        .select("accession", "kind", "content_sha256")
        .withColumnRenamed("content_sha256", "new_sha")
    )

    logger.info(f"[diff] Snapshot sizes: old={old_df.count()} records, new={new_df.count()} records")

    # --- Detect changed or missing hashes ---
    changed_df = (
        new_df.join(old_df, on=["accession", "kind"], how="fullouter")
        .where((F.col("old_sha").isNull()) | (F.col("new_sha").isNull()) | (F.col("old_sha") != F.col("new_sha")))
        .select("accession")
        .distinct()
    )

    changed_count = changed_df.count()
    logger.info(f"[diff] Changed accessions detected: {changed_count}")

    if changed_count == 0:
        logger.info("[diff] No hash differences detected — snapshots are identical.")
        return []

    # --- Collect changed accessions ---
    accessions = [r["accession"] for r in changed_df.collect() if r["accession"]]

    if debug:
        sample = [r["accession"] for r in changed_df.limit(10).collect()]
        logger.info(f"[debug] Example changed accessions: {sample}")

    # --- Map to taxids ---
    changed_taxids, seen = [], set()
    for acc in accessions:
        meta = acc_index.get(acc) or {}
        taxid = meta.get("species_taxid") or meta.get("taxid")
        if taxid and taxid not in seen:
            seen.add(taxid)
            changed_taxids.append(taxid)
            if debug:
                logger.info(f"[debug] accession={acc} → taxid={taxid}")
        elif debug and not taxid:
            logger.warning(f"[debug] accession={acc} → No taxid found in metadata.")

    logger.info(f"[diff] {len(accessions)} changed accessions → {len(changed_taxids)} unique taxids")
    return changed_taxids
