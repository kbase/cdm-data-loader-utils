"""DLT pipeline to import UniRef data."""

from datetime import UTC, datetime
from enum import StrEnum
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
from cdm_data_loaders.parsers.uniprot.uniref import (
    ENTRY_XML_TAG,
    parse_uniref_entry,
)
from cdm_data_loaders.pipelines.core import (
    PipelineCli,
    run_pipeline,
)
from cdm_data_loaders.readers.xml import build_xml_file_resource

PIPELINE_NAME: Final[str] = "uniref"
UNIREF_LOG_INTERVAL: Final[int] = 10000
VARIANT: Final[str] = "variant"
FIFTY: Final[str] = "50"
NINETY: Final[str] = "90"
HUNDRED: Final[str] = "100"


class UnirefVariantEnum(StrEnum):
    """Valid values for the Uniref variant."""

    FIFTY = FIFTY
    NINETY = NINETY
    HUNDRED = HUNDRED


UNIREF_VARIANTS: Final[list[str]] = [member.value for member in UnirefVariantEnum.__members__.values()]


class UnirefSettings(CtsSettings):
    """Configuration for running the UniRef import pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts={**CLI_SHORTCUTS, "variant": "v"},
    )

    variant: Annotated[
        UnirefVariantEnum,
        Field(
            description=f"Which UniRef variant to import. Choices: {UNIREF_VARIANTS}",
        ),
    ]
    log_interval: Annotated[
        LogInterval,
        Field(
            default=UNIREF_LOG_INTERVAL,
        ),
    ]


def parse_uniref(settings: UnirefSettings) -> DltResource:
    """Build the resource that parses the information from UniRef files.

    :param settings: config for running the pipeline.
    :type settings: UnirefSettings
    :return: resource yielding parsed UniRef entries
    :rtype: DltResource
    """
    # a single timestamp is used to mark every entity parsed in this run
    timestamp = datetime.now(UTC)
    return build_xml_file_resource(
        settings=settings,
        xml_tag=ENTRY_XML_TAG,
        parse_fn=lambda entry, file_path: parse_uniref_entry(
            entry=entry,
            timestamp=timestamp,
            file_path=file_path,
            uniref_variant=f"UniRef {settings.variant}",
        ),
        resource_name="parse_uniref",
    )


def run_uniref_pipeline(settings: UnirefSettings) -> LoadInfo | None:
    """Execute the Uniref pipeline.

    :param settings: config for running the pipeline.
    :type settings: UnirefSettings
    """
    return run_pipeline(
        settings=settings,
        resource=parse_uniref(settings),
        pipeline_kwargs={
            "pipeline_name": PIPELINE_NAME,
            "dataset_name": f"uniref_{settings.variant}",
        },
    )


cli = PipelineCli(UnirefSettings, run_uniref_pipeline)


if __name__ == "__main__":
    cli()
