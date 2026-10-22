"""Settings class for the XML to dictionary pipeline for the KBase CTS."""

from typing import Annotated, Final

from pydantic import Field
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import (
    NonEmptyStr,
    TableName,
)
from cdm_data_loaders.core.settings import (
    XML_CLI_SHORTCUTS,
    XmlIngestSettings,
    default_settings_with_shortcuts,
)

PIPELINE_NAME: Final[str] = "xmltodict"


class XmlToDictSettings(XmlIngestSettings):
    """Settings for the XML ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts=XML_CLI_SHORTCUTS,
    )

    table_name: TableName
    xml_tag: Annotated[
        NonEmptyStr,
        Field(
            description="XML tag to capture the contents of",
        ),
    ]
