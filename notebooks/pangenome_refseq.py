"""
Utility script to identify missing RefSeq assemblies relative to GTDB.

This script:
1. Reads a GTDB metastore table.
2. Removes GB_/RS_ prefixes from genome_id.
3. Downloads the latest RefSeq assembly summary.
4. Computes missing GCF assemblies.
5. Outputs two text files using Spark distributed write:
   - r214_assemblies
   - missing_refseq_ids
"""

from __future__ import annotations

import logging
import tempfile
import urllib.request
from pathlib import Path
from typing import TYPE_CHECKING

import click
from pyspark.sql.functions import regexp_replace

from berdl_notebook_utils.setup_spark_session import get_spark_session

if TYPE_CHECKING:
    from pyspark.sql import SparkSession

logger = logging.getLogger(__name__)

REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"


def download_refseq_summary(output_path: Path) -> Path:
    """
    Download RefSeq assembly summary file.
    """
    logger.info("Downloading RefSeq assembly summary from %s", REFSEQ_URL)
    urllib.request.urlretrieve(REFSEQ_URL, output_path)  # noqa: S310
    return output_path


def parse_refseq_gcf_ids(file_path: Path) -> list[str]:
    """
    Parse all GCF_ assembly accessions from the RefSeq summary file.
    """
    assembly_ids: list[str] = []

    with file_path.open(encoding="utf-8") as file:
        for line in file:
            if line.startswith("#"):
                continue

            accession = line.split("\t", 1)[0]
            if accession.startswith("GCF_"):
                assembly_ids.append(accession)

    return assembly_ids


@click.command()
@click.option(
    "--gtdb-table",
    required=True,
    help="Metastore table containing genome_id column",
)
@click.option(
    "--output-dir",
    required=True,
    help="Output directory (e.g. s3a://...) where text files will be written",
)
def main(gtdb_table: str, output_dir: str) -> None:
    """
    Run the missing RefSeq assembly detection pipeline.
    """
    logging.basicConfig(level=logging.INFO)

    spark: SparkSession = get_spark_session()

    # ------------------------------------------------------------------
    # 1. Read GTDB genome table
    # ------------------------------------------------------------------
    r214_df = spark.table(gtdb_table).select("genome_id").distinct()

    rm_prefix_df = (
        r214_df.withColumn(
            "assembly_id",
            regexp_replace("genome_id", r"^(GB_|RS_)", ""),
        )
        .select("assembly_id")
        .distinct()
    )

    logger.info("GTDB assemblies: %d", rm_prefix_df.count())

    # ------------------------------------------------------------------
    # 2. Download RefSeq summary securely
    # ------------------------------------------------------------------
    with tempfile.NamedTemporaryFile(delete=False) as tmp:
        summary_path = Path(tmp.name)

    download_refseq_summary(summary_path)

    # ------------------------------------------------------------------
    # 3. Parse RefSeq GCF IDs
    # ------------------------------------------------------------------
    refseq_ids = parse_refseq_gcf_ids(summary_path)

    refseq_df = spark.createDataFrame(
        [(x,) for x in refseq_ids],
        ["assembly_id"],
    )

    logger.info("RefSeq assemblies: %d", refseq_df.count())

    # ------------------------------------------------------------------
    # 4. Compute missing assemblies
    # ------------------------------------------------------------------
    missing_df = refseq_df.join(
        rm_prefix_df,
        on="assembly_id",
        how="left_anti",
    )

    logger.info("Missing RefSeq assemblies: %d", missing_df.count())

    # ------------------------------------------------------------------
    # 5. Distributed Spark text output
    # ------------------------------------------------------------------

    # Output 1: All GTDB assemblies
    rm_prefix_df.select("assembly_id").orderBy("assembly_id").coalesce(1).write.mode("overwrite").text(
        f"{output_dir}/r214_assemblies"
    )

    # Output 2: Missing RefSeq assemblies
    missing_df.select("assembly_id").orderBy("assembly_id").coalesce(1).write.mode("overwrite").text(
        f"{output_dir}/missing_refseq_ids"
    )

    logger.info("Output files successfully written to %s", output_dir)


if __name__ == "__main__":
    main()
