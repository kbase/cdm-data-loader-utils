"""DLT pipeline to import UniProt data."""

from datetime import UTC, datetime
from typing import Annotated, Final

from dlt.common.pipeline import LoadInfo
from dlt.extract import DltResource
from pydantic import Field
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import LogInterval
from cdm_data_loaders.core.settings import (
    CLI_SHORTCUTS,
    CtsSettings,
    default_settings_with_shortcuts,
)
from cdm_data_loaders.parsers.uniprot.uniprot_kb import ENTRY_XML_TAG, parse_uniprot_entry
from cdm_data_loaders.pipelines.core import (
    PipelineCli,
    run_pipeline,
)
from cdm_data_loaders.readers.xml import build_xml_file_resource

PIPELINE_NAME: Final[str] = "uniprot_kb"
UNIPROT_LOG_INTERVAL: Final[int] = 1000


class UniProtSettings(CtsSettings):
    """Configuration for running the UniProt KB import pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name="uniprot", cli_shortcuts=CLI_SHORTCUTS
    )

    log_interval: Annotated[
        LogInterval,
        Field(
            default=UNIPROT_LOG_INTERVAL,
        ),
    ]


def parse_uniprot(settings: UniProtSettings) -> DltResource:
    """Build the resource that parses the information from UniProt files.

    :param settings: config for running the pipeline.
    :type settings: UniProtSettings
    :return: resource yielding parsed UniProt entries
    :rtype: DltResource
    """
    # a single timestamp is used to mark every entity parsed in this run
    timestamp = datetime.now(UTC)
    return build_xml_file_resource(
        settings=settings,
        xml_tag=ENTRY_XML_TAG,
        parse_fn=lambda entry, file_path: parse_uniprot_entry(entry=entry, timestamp=timestamp, file_path=file_path),
        resource_name="parse_uniprot",
    )


def run_uniprot_pipeline(settings: UniProtSettings) -> LoadInfo | None:
    """Execute the UniProt KB pipeline."""
    return run_pipeline(
        settings=settings,
        resource=parse_uniprot(settings),
        pipeline_kwargs={
            "pipeline_name": PIPELINE_NAME,
            "dataset_name": PIPELINE_NAME,
        },
    )


cli = PipelineCli(UniProtSettings, run_uniprot_pipeline)


if __name__ == "__main__":
    cli()
