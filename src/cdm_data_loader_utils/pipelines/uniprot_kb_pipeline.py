"""DLT pipeline to import UniProt data."""

import datetime
from collections.abc import Generator
from pathlib import Path
from typing import Any

import click
import dlt
from dlt.extract.items import DataItemWithMeta
from pydantic import Field
from pydantic_settings import BaseSettings

from cdm_data_loader_utils.parsers.uniprot.uniprot_kb import ENTRY_XML_TAG, parse_uniprot_entry
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger
from cdm_data_loader_utils.utils.xml_utils import stream_xml_file

logger = get_cdm_logger()

APP_NAME = "uniprot_kb_importer"

VALID_DESTINATIONS = ["local_fs", "minio"]

DEFAULT_UNIPROT_DIR = "s3://cdm-lake/tenant-general-warehouse/kbase/datasets/uniprot/"
SOURCE_DIR = Path("/global_share") / "uniprot" / "derived" / "2025_03" / "uniprot_kb"


class Settings(BaseSettings):
    """Configuration for running the UniRef import pipeline."""

    input_dir: Path = Field()
    destination: str = Field()
    start_at: int = Field(0)
    timestamp: datetime.datetime = Field(datetime.datetime.now(tz=datetime.UTC))


@dlt.resource(name="parse_uniprot", parallelized=True)
def parse_uniprot(file_path: str | Path, current_timestamp: datetime.datetime) -> Generator[DataItemWithMeta, Any]:
    """Parse the information from a UniProt file.

    :param file_path: path to a uniprot XML file
    :type file_path: str | Path
    :return: parsed uniprot data sorted to the appropriate tables
    :rtype: Generator[DataItemWithMeta, Any]
    """
    for n, entry in enumerate(stream_xml_file(file_path, ENTRY_XML_TAG)):
        parsed_entry = parse_uniprot_entry(entry, current_timestamp, file_path)
        for table, rows in parsed_entry.items():
            yield dlt.mark.with_table_name(rows, table)
        if n % 1000 == 0:
            print(f"Processed {n} entries")


def run_pipeline(config: Settings) -> None:
    """Execute the pipeline.

    :param config: config for running the pipeline.
    :type config: Settings
    """
    for uniprot_file in sorted(config.input_dir.glob("*.xml.gz")):
        if config.start_at:
            # get the integer part of the file name
            f_int = uniprot_file.stem.replace("parts_", "")
            if int(f_int) < config.start_at:
                logger.info("Skipping %s", str(uniprot_file))
                continue
        logger.info("Reading from %s", str(uniprot_file))
        pipeline = dlt.pipeline(
            pipeline_name="uniprot_kb",
            destination=dlt.destination(config.destination, max_table_nesting=0),
            dataset_name="uniprot_kb",
        )
        load_info = pipeline.run(parse_uniprot(uniprot_file, config.timestamp), table_format="delta")
        logger.info("Work complete: %s ingested!", str(uniprot_file))
        logger.info(load_info)


@click.command()
@click.option("-s", "--start", type=str, default=0, help="UniProt file to start import at")
@click.option("-i", "--input_dir", default=str(SOURCE_DIR), help="Location of UniProtKB XML files to import")
@click.option(
    "-d",
    "--destination",
    type=click.Choice(VALID_DESTINATIONS),
    default="local_fs",
    help="Destination configuration to use for data output",
)
@click.option(
    "-o",
    "--output",
    type=str,
    default="",
    help="Location to save imported data to, if different from the default supplied by the destination config",
)
def cli(input_dir: str, destination: str, output: str | None, start: int) -> None:
    """CLI interface for the UniProt importer pipeline.

    :param input_dir: Location of the directory containing the UniProt XML files
    :type input_dir: str
    :param destination: destination configuration to use
    :type destination: str | None
    :param output: location in the object store to save files to
    :type output: str
    :param start: if provided, which file to start the import at
    :type start: int
    """
    # check whether there is a custom output location; if so, set it in the config
    if output:
        dlt.config[f"destination.{destination}.bucket_url"] = output

    runtime_config = Settings(
        input_dir=Path(input_dir),
        destination=destination,
        start_at=start,
        timestamp=datetime.datetime.now(tz=datetime.UTC),
    )
    run_pipeline(runtime_config)


if __name__ == "__main__":
    cli()
