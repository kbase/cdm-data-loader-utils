import logging
import urllib.request
from pathlib import Path
import click
from typing import List, Optional
from pyspark.sql import SparkSession, DataFrame
from berdl_notebook_utils.setup_spark_session import get_spark_session


logger = logging.getLogger(__name__)

REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"


def download_refseq_summary(output_path: Path) -> Path:
    logger.info("Downloading RefSeq assembly summary from %s", REFSEQ_URL)
    urllib.request.urlretrieve(REFSEQ_URL, output_path)
    return output_path


def parse_refseq_gcf_ids(file_path: Path) -> List[str]:
    assembly_ids: List[str] = []

    with open(file_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue

            accession = line.split("\t", 1)[0]

            if accession.startswith("GCF_"):
                assembly_ids.append(accession)

    logger.info("Parsed %d GCF assemblies", len(assembly_ids))
    return assembly_ids


def build_refseq_df(spark: SparkSession, assembly_ids: List[str]) -> DataFrame:
    return spark.createDataFrame(
        [(x,) for x in assembly_ids],
        ["assembly_id"],
    )


def compute_missing_refseq(
    refseq_df: DataFrame,
    existing_df: DataFrame,
) -> DataFrame:
    return refseq_df.join(
        existing_df.select("assembly_id"),
        on="assembly_id",
        how="left_anti",
    )


def read_existing_df(
    spark: SparkSession,
    existing_table: Optional[str],
    existing_path: Optional[str],
) -> DataFrame:
    if existing_table:
        logger.info("Reading existing table from metastore: %s", existing_table)
        return spark.table(existing_table)

    if existing_path:
        logger.info("Reading existing Delta table from path: %s", existing_path)
        return spark.read.format("delta").load(existing_path)

    raise ValueError("Either --existing-table or --existing-path must be provided.")


def write_output(df: DataFrame, output_path: str) -> None:
    logger.info("Writing missing assemblies to %s", output_path)

    (df.coalesce(1).write.format("delta").mode("overwrite").save(output_path))


# -------------------------------------------------------------------------
# Pipeline Orchestration
# -------------------------------------------------------------------------
def run_pipeline(
    spark: SparkSession,
    existing_table: Optional[str],
    existing_path: Optional[str],
    output_path: str,
) -> None:
    summary_path = Path("assembly_summary_refseq.txt")

    download_refseq_summary(summary_path)

    assembly_ids = parse_refseq_gcf_ids(summary_path)

    refseq_df = build_refseq_df(spark, assembly_ids)

    existing_df = read_existing_df(
        spark,
        existing_table=existing_table,
        existing_path=existing_path,
    )

    missing_df = compute_missing_refseq(refseq_df, existing_df)

    missing_count = missing_df.count()
    logger.info("Missing RefSeq assemblies: %d", missing_count)

    write_output(missing_df, output_path)


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------
@click.command()
@click.option("--existing-table", help="Existing table in metastore containing RefSeq assemblies")
@click.option("--existing-path", help="Existing Delta path containing RefSeq assemblies")
@click.option("--output-path", required=True, help="Output path for missing RefSeq assemblies (Delta format)")
def main(
    existing_table: Optional[str],
    existing_path: Optional[str],
    output_path: str,
) -> None:
    spark = get_spark_session()

    run_pipeline(
        spark=spark,
        existing_table=existing_table,
        existing_path=existing_path,
        output_path=output_path,
    )


def cli():
    logging.basicConfig(level=logging.INFO)
    main(standalone_mode=False)


if __name__ == "__main__":
    cli()
