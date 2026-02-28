"""

PYTHONPATH=src python3 -m cdm_data_loader_utils.parsers.prefix_audit \
  --parquet "dataset/part-00000-0a0d0261-1fee-477d-90d8-1df048058fbf-c000.snappy.parquet" \
  --idmapping "dataset/ECOLI_83333_idmapping.dat" \
  --registry "dataset/registry.json"

"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Set
import gzip

from pyspark.sql import SparkSession
from pyspark.sql.functions import lower, col


# -----------------------------
# Load parquet db prefixes
# -----------------------------
def load_parquet_prefixes(
    spark: SparkSession,
    parquet_path: str,
) -> Set[str]:
    df = spark.read.parquet(parquet_path)

    prefixes = df.select(lower(col("db")).alias("db")).distinct().collect()
    return {row["db"] for row in prefixes if row["db"]}


# -----------------------------
# Load ID_type from idmapping.dat
# -----------------------------


def load_idmapping_prefixes(
    idmapping_path: Path,
) -> Set[str]:
    """
    Extract unique ID_type values (second column) from idmapping.dat.
    Supports both plain .dat files.
    """

    prefixes: Set[str] = set()

    # Support gz or plain text
    open_fn = gzip.open if idmapping_path.suffix == ".gz" else open

    with open_fn(idmapping_path, "rt") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                continue

            id_type = parts[1].strip().lower()
            if id_type:
                prefixes.add(id_type)

    return prefixes


def load_registry_prefixes(
    registry_path: Path,
) -> Set[str]:
    """
    Load canonical prefixes from Bioregistry registry.json.
    The registry file is expected to be a dict mapping
    canonical_prefix -> metadata.
    """

    if not registry_path.exists():
        raise FileNotFoundError(f"Registry prefix not found: {registry_path}")

    with registry_path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, dict):
        raise ValueError("Invalid registry format")

    return {prefix.strip().lower() for prefix in data if isinstance(prefix, str) and prefix.strip()}


# -----------------------------
# Compare universes
# -----------------------------
def compare_prefix_sets(
    parquet_set: Set[str],
    idmapping_set: Set[str],
) -> None:
    """
    Compare parquet prefix universe with idmapping universe.
    """

    print("\n=== Parquet vs IDMapping ===\n")

    total_parquet = len(parquet_set)
    total_idmapping = len(idmapping_set)

    in_both = parquet_set & idmapping_set
    only_in_parquet = parquet_set - idmapping_set
    only_in_idmapping = idmapping_set - parquet_set

    print(f"Total parquet prefixes: {total_parquet}")
    print(f"Total idmapping prefixes: {total_idmapping}\n")

    print(f"Shared prefixes: {len(in_both)} ({len(in_both) / total_parquet:.1%}")
    for p in sorted(in_both):
        print(f"{p}")

    print(f"\nPresent only in parquet: {len(only_in_parquet)}")
    for p in sorted(only_in_parquet):
        print(f"{p}")

    print(f"\nPresent only in idmapping: {len(only_in_idmapping)}")
    for p in sorted(only_in_idmapping):
        print(f"{p}")


def compare_with_registry(
    parquet_set: Set[str],
    registry_set: Set[str],
) -> None:
    print("\n=== Bioregistry ===\n")

    valid = parquet_set & registry_set
    missing = parquet_set - registry_set

    print(f"Valid namespaces in Bioregistry: {len(valid)}")
    for p in sorted(valid):
        print("  ", p)

    print(f"\nNot found in Bioregistry: {len(missing)}")
    for p in sorted(missing):
        print("  ", p)


# -----------------------------
# CLI entry
# -----------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Audit identifier prefix universe against UniProt idmapping.")
    parser.add_argument(
        "--parquet",
        required=True,
    )
    parser.add_argument(
        "--idmapping",
        required=True,
    )
    parser.add_argument(
        "--registry",
        required=True,
    )

    args = parser.parse_args()
    spark = SparkSession.builder.appName("PrefixAudit").getOrCreate()

    parquet_set = load_parquet_prefixes(spark, args.parquet)
    idmapping_set = load_idmapping_prefixes(Path(args.idmapping))
    registry_set = load_registry_prefixes(Path(args.registry))

    compare_prefix_sets(parquet_set, idmapping_set)
    compare_with_registry(parquet_set, registry_set)

    spark.stop()


if __name__ == "__main__":
    main()
