# ==========================================================
# Compare two Delta snapshot tags and extract changed TaxIDs

# This script automatically registers the Delta table if not found,
# then compares two hash snapshots (20250930 vs 20251014)
# to detect TaxIDs whose assemblies have been added, removed, or updated.
# ==========================================================

# import os
# import json
# import click
# from pathlib import Path
# from pyspark.sql.utils import AnalysisException

# from refseq_pipeline.core.spark_delta import build_spark
# from refseq_pipeline.core.refseq_io import load_refseq_assembly_index
# from refseq_pipeline.core.hashes_diff import diff_hash_and_get_changed_taxids


# def run_diff_changed_taxids(
#         database: str,
#         hash_table: str,
#         new_tag: str,
#         old_tag: str) -> list[str]:
#     """
#     Run diff between two Delta table snapshots and return changed TaxIDs.
#     Automatically registers the Delta table if it is not found.

#     """

#     # Initialize Spark
#     spark = build_spark(database)

#     # Ensure the Delta table is registered
#     try:
#         spark.sql(f"DESCRIBE TABLE {database}.{hash_table}")
#         print(f"[diff] Found existing table: {database}.{hash_table}")
#     except AnalysisException:
#         PROJECT_ROOT = Path(__file__).resolve().parents[2]
#         #delta_path = os.path.abspath(f"delta_data/refseq/{database}/{hash_table}")
#         delta_path = PROJECT_ROOT / "delta_data" / "refseq" / database / hash_table
#         print(f"[register] Table not found. Registering from: {delta_path}")

#         spark.sql(f"CREATE DATABASE IF NOT EXISTS {database}")

#         spark.sql(f"""
#             CREATE TABLE IF NOT EXISTS {database}.{hash_table}
#             USING DELTA
#             LOCATION '{delta_path}'
#         """)

#         print(f"[register] Table {database}.{hash_table} successfully registered.")

#     # Load the latest RefSeq assembly index (for accession-to-TaxID mapping)
#     acc_index = load_refseq_assembly_index()

#     # Compute changed TaxIDs between the two tags
#     taxids = diff_hash_and_get_changed_taxids(
#         spark,
#         database,
#         hash_table,
#         acc_index,
#         tag_new=new_tag,
#         tag_old=old_tag,
#     )
#     return taxids


# @click.command()
# @click.option(
#     "--database",
#     required=True,
#     help="Spark SQL database name where the Delta table is stored (e.g., 'refseq_api')."
# )
# @click.option(
#     "--hash-table",
#     default="assembly_hashes",
#     show_default=True,
#     help="Delta table name containing hash snapshots."
# )
# @click.option(
#     "--new-tag",
#     required=True,
#     help="Tag for the newer snapshot (e.g., 20251014)."
# )
# @click.option(
#     "--old-tag",
#     required=True,
#     help="Tag for the older snapshot (e.g., 20250930)."
# )

# @click.option(
#     "--out-path",
#     default=None,
#     type=click.Path(),
#     help="Path to save changed TaxIDs as a JSON file."
# )

# def main(database, hash_table, new_tag, old_tag, out_path):
#     """
#     Compare two Delta snapshots (by tag) and output a list of changed TaxIDs.
#     A TaxID is considered changed if:
#       - a new assembly was added under the TaxID
#       - an assembly was removed under the TaxID
#       - or an existing assembly’s hash changed
#     """
#     try:
#         taxids = run_diff_changed_taxids(database, hash_table, new_tag, old_tag)

#         if not taxids:
#             click.echo(f"[diff] No changed TaxIDs found between {old_tag} and {new_tag}.")
#             return 
        
        
#         PROJECT_ROOT = Path(__file__).resolve().parents[2]
#         if out_path is None:
#             output_dir = PROJECT_ROOT / "compare_snapshot_data" / f"{old_tag}_vs_{new_tag}"
#             out_path = output_dir / "changed_taxids.json"
#         else:
#             output_dir = Path(os.path.dirname(out_path)) if os.path.dirname(out_path) else PROJECT_ROOT / "compare_snapshot_data"
#             out_path = Path(out_path)

#         os.makedirs(output_dir, exist_ok=True)

#         # Save JSON
#         with open(out_path, "w", encoding="utf-8") as f:
#             json.dump(taxids, f, indent=2)
#         click.echo(f"[diff] {len(taxids)} changed TaxIDs written to {out_path}")

#     except Exception as e:
#         click.echo(f"[diff] ERROR: {e}", err=True)
#         raise SystemExit(1)


# if __name__ == "__main__":
#     main()

import os
import json
import click
from pathlib import Path

from refseq_pipeline.core.spark_delta import build_spark
from refseq_pipeline.core.refseq_io import load_refseq_assembly_index
from refseq_pipeline.core.hashes_diff import diff_hash_and_get_changed_taxids


def run_diff_changed_taxids(
        database: str,
        hash_table: str,
        new_tag: str,
        old_tag: str) -> list[str]:
    """
    Compare hash snapshots between two tags and return the list of changed TaxIDs.
    Avoids Spark's AnalysisException by ensuring the Delta table is registered
    *before* calling DESCRIBE TABLE.
    """

    spark = build_spark(database)

    # ---------- Locate Delta path ----------
    PROJECT_ROOT = Path(__file__).resolve().parents[2]
    delta_path = PROJECT_ROOT / "delta_data" / "refseq" / database / hash_table

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
@click.option(
    "--database", required=True,
    help="Spark database where the Delta table is stored (e.g., 'refseq_api')."
)
@click.option(
    "--hash-table", default="assembly_hashes", show_default=True,
    help="Delta table name containing hash snapshots."
)
@click.option(
    "--new-tag", required=True,
    help="Tag for the newer snapshot (e.g., 20251014)."
)
@click.option(
    "--old-tag", required=True,
    help="Tag for the older snapshot (e.g., 20250930)."
)
@click.option(
    "--out-path", default=None, type=click.Path(),
    help="Path to save changed TaxIDs JSON. If not supplied, auto-writes into compare_snapshot_data/<old>_vs_<new>/"
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

python -m refseq_pipeline.cli.diff_changed_taxids \
  --database refseq_api \
  --hash-table assembly_hashes \
  --new-tag 20251014 \
  --old-tag 20250930 \
  --out-path changed_taxids_20251014_vs_20250930.json

  taxid can be changed if:
  - new assembly added under taxid
  - assembly removed under taxid
  - assembly under taxid updated (hash changed)

"""

