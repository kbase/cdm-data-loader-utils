# ==========================================================
# Compare two Delta snapshot tags and extract changed TaxIDs

# This script automatically registers the Delta table if not found,
# then compares two hash snapshots (20250930 vs 20251014)
# to detect TaxIDs whose assemblies have been added, removed, or updated.
# ==========================================================


import json
import os
from pathlib import Path

import click
from refseq_pipeline.core.hashes_diff import diff_hash_and_get_changed_taxids
from refseq_pipeline.core.refseq_io import load_refseq_assembly_index
from refseq_pipeline.core.spark_delta import build_spark


def run_diff_changed_taxids(database: str, hash_table: str, new_tag: str, old_tag: str) -> list[str]:
    """
    Compare hash snapshots between two tags and return the list of changed TaxIDs.
    """
    spark = build_spark(database)

    # ---------- Locate Delta path ----------
    # PROJECT_ROOT = Path(__file__).resolve().parents[2]
    # delta_path = PROJECT_ROOT / "delta_data" / "refseq" / database / hash_table

    project_root = Path("/global_share/alinawang/cdm-data-loader-utils")
    delta_path = project_root / "delta_data" / "refseq" / "refseq_api" / "assembly_hashes"

    if not delta_path.exists():
        raise RuntimeError(
            f"[diff] Delta table path does not exist:\n  {delta_path}\n"
            f"Make sure you have written snapshots before running diff."
        )

    # ---------- Register table cleanly (no exceptions) ----------
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {database}")
    spark.sql(f"""
        CREATE TABLE IF NOT EXISTS {database}.{hash_table}
        USING DELTA
        LOCATION '{delta_path}'
    """)

    print(f"[register] Ensured table is registered at {delta_path}")

    # ---------- Load accession → taxid index ----------
    acc_index = load_refseq_assembly_index()

    # ---------- Compute changed TaxIDs ----------
    taxids = diff_hash_and_get_changed_taxids(
        spark,
        database,
        hash_table,
        acc_index,
        tag_new=new_tag,
        tag_old=old_tag,
    )
    return taxids


@click.command()
@click.option("--database", required=True, help="Spark database where the Delta table is stored.")
@click.option(
    "--hash-table", default="assembly_hashes", show_default=True, help="Delta table name containing hash snapshots."
)
@click.option("--new-tag", required=True, help="Tag for the newer snapshot (e.g., 20251014).")
@click.option("--old-tag", required=True, help="Tag for the older snapshot (e.g., 20250930).")
@click.option(
    "--out-path",
    default=None,
    type=click.Path(),
    help="Path to save changed TaxIDs JSON. If not supplied, auto-writes into compare_snapshot_data/<old>_vs_<new>/",
)
def main(database, hash_table, new_tag, old_tag, out_path):
    """CLI entrypoint for comparing snapshot hashes and exporting changed TaxIDs."""
    try:
        taxids = run_diff_changed_taxids(database, hash_table, new_tag, old_tag)

        if not taxids:
            click.echo(f"[diff] No changed TaxIDs found between {old_tag} and {new_tag}.")
            return

        PROJECT_ROOT = Path(__file__).resolve().parents[2]

        if out_path is None:
            output_dir = PROJECT_ROOT / "compare_snapshot_data" / f"{old_tag}_vs_{new_tag}"
            out_path = output_dir / "changed_taxids.json"
        else:
            output_dir = Path(out_path).parent
            out_path = Path(out_path)

        os.makedirs(output_dir, exist_ok=True)

        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(taxids, f, indent=2)

        click.echo(f"[diff] {len(taxids)} changed TaxIDs written to {out_path}")

    except Exception as e:
        click.echo(f"[diff] ERROR: {e}", err=True)
        raise SystemExit(1)


if __name__ == "__main__":
    main()


"""

  PYTHONPATH=src/parsers \
python -m refseq_pipeline.cli.diff_changed_taxids \
  --database refseq_api \
  --hash-table assembly_hashes \
  --old-tag 20250930 \
  --new-tag 20251120 \
  --out-path compare_snapshot_data/20250930_vs_20251120/changed_taxids.json

  taxid can be changed if:
  - new assembly added under taxid
  - assembly removed under taxid
  - assembly under taxid updated (hash changed)

"""
