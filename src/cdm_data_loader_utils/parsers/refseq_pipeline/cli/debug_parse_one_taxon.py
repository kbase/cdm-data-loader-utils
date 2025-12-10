import click
from refseq_pipeline.core.cdm_parse import parse_reports
from refseq_pipeline.core.datasets_api import fetch_reports_by_taxon
from refseq_pipeline.core.driver import process_and_write_reports
from refseq_pipeline.core.spark_delta import build_spark

"""

    PYTHONPATH=src/parsers python src/parsers/refseq_pipeline/cli/debug_parse_one_taxon.py \
  --taxid 752788 \
  --database refseq_api \
  --table assembly_stats \
  --mode overwrite \
  --data-dir /global_share/alinawang/cdm-data-loader-utils/delta_data \
  --preview-only

 """


@click.command()
@click.option("--taxid", required=True, help="NCBI Taxonomy ID to fetch.")
@click.option("--database", default="refseq_api", show_default=True)
@click.option("--table", default="assembly_stats", show_default=True)
@click.option("--mode", default="overwrite", type=click.Choice(["append", "overwrite"]))
@click.option("--data-dir", default=None, help="Optional external delta path.")
@click.option("--preview-only", is_flag=True, help="Only show parsed results, don't write.")
def main(taxid, database, table, mode, data_dir, preview_only):
    """
    Debug parser for a single TaxID.

    The function is designed to:
    - Fetch taxonomic genome reports from NCBI for a specific TaxID
    - Parse them into a structured Spark DataFrame using CDM schema
    - Optionally preview the parsed result in the console
    - Optionally write the parsed result to a Delta Lake table

    """
    # Initialize Spark Session
    spark = build_spark(database)
    print(f"[DEBUG] Fetching reports for TaxID: {taxid}")

    # Fetch NCBI genome reports for this taxon
    reports = list(fetch_reports_by_taxon(taxon=taxid, debug=True))
    if not reports:
        print("[DEBUG] No reports found.")
        return

    # Parse the reports into Spark DataFrame using CDM structure
    df = parse_reports(reports, return_spark=True, spark=spark)
    print(f"[DEBUG] Parsed {df.count()} CDM rows for TaxID={taxid}")

    # Preview the parsed DataFrame — convert to string and fill nulls
    df.select([df[c].cast("string").alias(c) for c in df.columns]).na.fill("None").show(truncate=False)

    # Write parsed data to Delta Lake
    if not preview_only:
        print(f"[debug] Writing to Delta table: {database}.{table}")
        process_and_write_reports(
            spark=spark,
            reports=reports,
            database=database,
            table=table,
            mode=mode,
            data_dir=data_dir,
        )


if __name__ == "__main__":
    main()
