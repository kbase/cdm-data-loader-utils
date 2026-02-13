import logging
import urllib.request
from pathlib import Path

import click
from pyspark.sql import SparkSession
from pyspark.sql.functions import regexp_replace

from berdl_notebook_utils.setup_spark_session import get_spark_session


logger = logging.getLogger(__name__)

REFSEQ_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"


def download_refseq_summary(output_path: Path) -> Path:
    logger.info("Downloading RefSeq assembly summary from %s", REFSEQ_URL)
    urllib.request.urlretrieve(REFSEQ_URL, output_path)  # noqa: S310
    return output_path


def parse_refseq_gcf_ids(file_path: Path) -> list[str]:
    assembly_ids: list[str] = []

    with open(file_path, encoding="utf-8") as f:
        for line in f:
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
    help="Directory where output text files will be written",
)
def main(gtdb_table: str, output_dir: str) -> None:
    logging.basicConfig(level=logging.INFO)

    spark: SparkSession = get_spark_session()

    # Read the GTDB genome table:
    r214_df = spark.table(gtdb_table).select("genome_id").distinct()

    rm_prefix_df = (
        r214_df.withColumn(
            "assembly_id",
            regexp_replace("genome_id", r"^(GB_|RS_)", ""),
        )
        .select("assembly_id")
        .distinct()
    )

    logger.info("Total GTDB assemblies: %d", rm_prefix_df.count())

    # Download RefSeq summary in BERDL temp directory
    summary_path = Path("/tmp/assembly_summary_refseq.txt")
    download_refseq_summary(summary_path)

    # Parse RefSeq GCF IDs
    refseq_ids = parse_refseq_gcf_ids(summary_path)

    refseq_df = spark.createDataFrame(
        [(x,) for x in refseq_ids],
        ["assembly_id"],
    )

    # Compute missing values in GTDB
    missing_df = refseq_df.join(
        rm_prefix_df,
        on="assembly_id",
        how="left_anti",
    )

    logger.info("Missing RefSeq assemblies: %d", missing_df.count())

    # Prepare output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Output 1: Write GTDB existing assemblies
    gtdb_ids = rm_prefix_df.select("assembly_id").rdd.map(lambda x: x[0]).collect()

    with open(output_path / "r214_assemblies.txt", "w") as f:
        for assembly_id in gtdb_ids:
            f.write(f"{assembly_id}\n")

    # Output 2: Write missing RefSeq assemblies
    missing_ids = missing_df.select("assembly_id").rdd.map(lambda x: x[0]).collect()

    with open(output_path / "missing_refseq_ids.txt", "w") as f:
        for assembly_id in missing_ids:
            f.write(f"{assembly_id}\n")

    logger.info("Output files written to %s", output_dir)


if __name__ == "__main__":
    main()
