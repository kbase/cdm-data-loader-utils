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

from cdm_data_loader_utils.parsers.uniprot.uniref import UNIREF_URL, UNIREF_VARIANTS, parse_uniref_entry
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger
from cdm_data_loader_utils.utils.xml_utils import stream_xml_file

logger = get_cdm_logger()

APP_NAME = "uniref_importer"

VALID_DESTINATIONS = ["local_fs", "minio"]

DEFAULT_UNIPROT_DIR = "s3://cdm-lake/tenant-general-warehouse/kbase/datasets/uniprot/"
UNIREF_DIR = Path("/global_share") / "uniprot" / "derived" / "2025_03" / "uniref"


class Settings(BaseSettings):
    """Configuration for running the UniRef import pipeline."""

    input_dir: Path = Field()
    destination: str = Field()
    uniref_variant: int = Field()
    start_at: int = Field(0)
    timestamp: datetime.datetime = Field(datetime.datetime.now(tz=datetime.UTC))


@dlt.resource(name="parse_uniref", parallelized=True)
def parse_uniref(
    file_path: str | Path, current_timestamp: datetime.datetime, uniref_value: int
) -> Generator[DataItemWithMeta, Any]:
    """Parse the information from a UniProt entry.

    :param entry: _description_
    :type entry: Element
    :return: _description_
    :rtype: _type_
    """
    for n, entry in enumerate(stream_xml_file(file_path, f"{{{UNIREF_URL}}}entry")):
        parsed_entry = parse_uniref_entry(entry, current_timestamp, f"UniRef {uniref_value}", file_path)
        for table, rows in parsed_entry.items():
            yield dlt.mark.with_table_name(rows, table)
        if n + 1 % 10000 == 0:
            print(f"Processed {n + 1} entries")


def run_pipeline(config: Settings) -> None:
    """Execute the pipeline.

    :param config: config for running the pipeline.
    :type config: Settings
    """
    for uniref_file in sorted(config.input_dir.glob("*.xml.gz")):
        if config.start_at:
            # get the integer part of the file name
            f_int = uniref_file.stem.replace("part_", "")
            if int(f_int) < config.start_at:
                logger.info("Skipping %s", str(uniref_file))
                continue
        logger.info("Reading from %s", str(uniref_file))
        pipeline = dlt.pipeline(
            pipeline_name=f"uniref_{config.uniref_variant}",
            destination=dlt.destination(config.destination, max_table_nesting=0),
            dataset_name="uniprot_kb",
        )
        load_info = pipeline.run(
            parse_uniref(uniref_file, config.timestamp, config.uniref_variant), table_format="delta"
        )
        logger.info("Work complete!")
        logger.info(load_info)


@click.command()
@click.option(
    "-n", "--uniref", required=True, type=click.Choice(UNIREF_VARIANTS), help="Which UniRef variant to import"
)
@click.option("-s", "--start", type=int, default=0, help="File to start import at")
@click.option("-i", "--input_dir", default=str(UNIREF_DIR), help="Location of UniRef XML files to import")
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
def cli(input_dir: str, destination: str, output: str | None, start: int, uniref: int) -> None:
    """CLI interface for the UniRef importer pipeline.

    :param input_dir: Location of the directory containing the UniRef XML files
    :type input_dir: str
    :param destination: destination configuration to use
    :type destination: str | None
    :param output: location in the object store to save files to
    :type output: str
    :param start: if provided, which file to start the import at
    :type start: int
    :param uniref: which UniRef dataset to import (50 / 90 / 100)
    :type uniref: int
    """
    # check whether there is a custom output location; if so, set it in the config
    if output:
        dlt.config[f"destination.{destination}.bucket_url"] = output

    runtime_config = Settings(
        input_dir=Path(input_dir),
        destination=destination,
        uniref_variant=uniref,
        start_at=start,
        timestamp=datetime.datetime.now(tz=datetime.UTC),
    )
    run_pipeline(runtime_config)


if __name__ == "__main__":
    cli()
