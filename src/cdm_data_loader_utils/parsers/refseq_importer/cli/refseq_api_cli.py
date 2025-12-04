"""
Example usage:

PYTHONPATH=src/parsers \
python -m cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli \
  --taxid "224325, 2741724, 193567" \
  --database refseq_api \
  --mode overwrite \
  --debug \
  --unique-per-taxon \
  --data-dir /global_share/alinawang/cdm-data-loader-utils/output/taxon_data

"""

import os
import re
from typing import Optional

import click

from cdm_data_loader_utils.parsers.refseq_importer.core.cdm_builders import build_cdm_datasource
from cdm_data_loader_utils.parsers.refseq_importer.core.spark_delta import build_spark, write_delta
from cdm_data_loader_utils.parsers.refseq_importer.core.tables_finalize import finalize_tables, write_and_preview
from cdm_data_loader_utils.parsers.refseq_importer.core.taxon_processing import process_taxon

# ---------------- Helpers ----------------


def parse_taxid_args(taxid_arg: Optional[str], taxid_file: Optional[str]) -> list[str]:
    """
    Parse and collect valid numeric TaxIDs from command-line arguments and file.
    """

    # empty list to collect taxids，avoid the duplicate TaxIDs
    taxids: list[str] = []

    # Parse taxid argument: --taxid "224325, 2741724"
    if taxid_arg:
        id_list = taxid_arg.split(",")  # separate them into a list using commas
        for num in id_list:
            # Keep only digits
            id = re.sub(r"\D+", "", num.strip())  # Remove all non-numeric characters
            if id:
                taxids.append(id)

    # Parse --taxid-file
    if taxid_file:
        if not os.path.exists(taxid_file):
            raise click.BadParameter(f"Path '{taxid_file}' does not exist.", param_hint="--taxid-file")
        with open(taxid_file, "r", encoding="utf-8") as f:
            for line in f:
                id = re.sub(r"\D+", "", line.strip())
                if id:
                    taxids.append(id)

    # Deduplicate while preserving order
    seen = set()
    unique_taxids = []
    for id in taxids:
        if id not in seen:
            seen.add(id)
            unique_taxids.append(id)

    return unique_taxids


# ---------------- CLI ----------------
@click.command()
@click.option("--taxid", required=True, help="Comma-separated NCBI TaxIDs, e.g. '224325,2741724'.")
@click.option("--api-key", default=None, help="Optional NCBI API key (increases rate limits).")
@click.option("--database", default="refseq_api", show_default=True, help="Delta schema/database.")
@click.option(
    "--mode",
    default="overwrite",
    type=click.Choice(["overwrite", "append"]),
    show_default=True,
    help="Write mode for Delta tables.",
)
@click.option(
    "--debug/--no-debug", default=False, show_default=True, help="Print per-record parsed fields for troubleshooting."
)
@click.option(
    "--allow-genbank-date/--no-allow-genbank-date",
    default=False,
    show_default=True,
    help="Allow using GenBank submissionDate as fallback for RefSeq created date.",
)
@click.option(
    "--unique-per-taxon/--all-assemblies",
    default=False,
    show_default=True,
    help="Keep only one assembly per taxon (latest by release_date).",
)
@click.option(
    "--data-dir",
    required=True,
    help="Base directory for Delta tables, e.g. /Users/yuewang/Documents/RefSeq/entities_data",
)
def cli(taxid, api_key, database, mode, debug, allow_genbank_date, unique_per_taxon, data_dir):
    main(
        taxid=taxid,
        api_key=api_key,
        database=database,
        mode=mode,
        debug=debug,
        allow_genbank_date=allow_genbank_date,
        unique_per_taxon=unique_per_taxon,
        data_dir=data_dir,
    )


def main(
    taxid,
    api_key,
    database,
    mode,
    debug,
    allow_genbank_date: bool = False,
    unique_per_taxon: bool = False,
    data_dir: str | None = None,
):
    spark = build_spark(database)

    df_ds = build_cdm_datasource(spark)
    write_delta(spark, df_ds, database, "datasource", mode, data_dir=data_dir)

    entities, collections, names, identifiers = [], [], [], []
    seen = set()

    taxids = [t.strip() for t in taxid.split(",") if t.strip()]
    print(f"Using TaxIDs: {taxids}")

    # ---- Process each taxon ----
    for tx in taxids:
        print(f"Fetching taxon={tx}")

        e, c, n, i = process_taxon(
            spark,
            tx,
            api_key,
            debug,
            allow_genbank_date,
            unique_per_taxon,
            seen,
        )
        entities.extend(e)
        collections.extend(c)
        names.extend(n)
        identifiers.extend(i)

    # ---- FINALIZE with Spark ----
    df_entity, df_coll, df_name, df_ident = finalize_tables(
        spark,
        entities,
        collections,
        names,
        identifiers,
    )

    # ---- Write results ----
    write_and_preview(
        spark,
        database,
        mode,
        df_entity,
        df_coll,
        df_name,
        df_ident,
        data_dir=data_dir,
    )


if __name__ == "__main__":
    cli()
