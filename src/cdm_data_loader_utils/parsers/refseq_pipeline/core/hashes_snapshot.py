from typing import Any
from datetime import datetime, timezone
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StructField, StringType

from refseq_pipeline.core.refseq_io import fetch_annotation_hash, fetch_md5_checksums, text_sha256
from refseq_pipeline.core.spark_delta import write_delta_table
from refseq_pipeline.core.config import DEFAULT_HASH_TABLE


# Delta schema for hash snapshots
HASH_SCHEMA = StructType(
    [
        StructField("accession", StringType(), True),
        StructField("ftp_path", StringType(), True),
        StructField("kind", StringType(), True),
        StructField("content_sha256", StringType(), True),
        StructField("retrieved_at", StringType(), True),
        StructField("tag", StringType(), True),
    ]
)


def snapshot_hashes_for_accessions(
    accessions: list[str],
    acc_index: dict[str, dict[str, str]],
    kind: str = "auto",
    timeout: int = 30,
    spark: SparkSession = None,
    fast_mode: bool = False,
):
    """
    Build a Spark DataFrame of hash snapshots for given accessions.

    Args:
        accessions: List of accession strings.
        acc_index: accession→metadata mapping (ftp_path, taxid, etc.)
        kind: "annotation", "md5", or "auto".
        timeout: Timeout for fetching FTP files.
        spark: SparkSession
        fast_mode: If True, only fetch md5 to reduce time.

    Returns:
        Spark DataFrame with HASH_SCHEMA.
    """
    if spark is None:
        raise ValueError("You must pass a SparkSession to snapshot_hashes_for_accessions().")

    ts = datetime.now(timezone.utc).isoformat(timespec="seconds")
    rows: list[dict[str, Any]] = []

    total, success, failed = 0, 0, 0

    for acc in accessions:
        total += 1
        meta = acc_index.get(acc) or {}
        ftp = (meta.get("ftp_path") or "").rstrip("/")
        if not ftp:
            failed += 1
            continue

        content = None
        actual_kind = None

        try:
            if fast_mode:
                # Fast mode: only use md5
                txt = fetch_md5_checksums(ftp, timeout=timeout)
                if txt:
                    content, actual_kind = txt, "md5"
            else:
                # Try annotation first, fallback to md5
                if kind in ("annotation", "auto"):
                    txt = fetch_annotation_hash(ftp, timeout=timeout)
                    if txt:
                        content, actual_kind = txt, "annotation"
                if content is None and kind in ("md5", "auto"):
                    txt = fetch_md5_checksums(ftp, timeout=timeout)
                    if txt:
                        content, actual_kind = txt, "md5"

        except Exception as e:
            print(f"[hash] Failed to fetch hash for {acc}: {e}")
            failed += 1
            continue

        if content is None:
            failed += 1
            continue

        rows.append(
            {
                "accession": acc,
                "ftp_path": ftp,
                "kind": actual_kind,
                "content_sha256": text_sha256(content),
                "retrieved_at": ts,
                "tag": None,
            }
        )
        success += 1

    print(f"[hash] Finished snapshot: total={total}, success={success}, failed={failed}")

    if not rows:
        return spark.createDataFrame([], schema=HASH_SCHEMA)

    return spark.createDataFrame(rows, schema=HASH_SCHEMA)


def write_hash_snapshot(
    spark: SparkSession,
    sdf,
    database: str,
    table: str = DEFAULT_HASH_TABLE,
    mode: str = "append",
    data_dir: str | None = None,
):
    """
    Write the hash snapshot Spark DataFrame into Delta table.
    Automatically handles empty dataframes and non-Delta locations.
    """
    if sdf is None or sdf.rdd.isEmpty():
        print("[hash] empty snapshot; skip write.")
        return

    if "tag" not in sdf.columns:
        from pyspark.sql import functions as F

        sdf = sdf.withColumn("tag", F.lit(None))

    try:
        write_delta_table(sdf, spark, database, table, mode=mode, data_dir=data_dir)
    except Exception as e:
        print(f"[hash] Failed to write hash snapshot to {database}.{table}: {e}")
        raise
