"""Settings class for the XML to dictionary pipeline for the KBase CTS."""

from typing import Annotated, Final

from pydantic import Field
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import (
    DEFAULT_XML_FILE_GLOB,
    FILE_GLOB,
    FileGlob,
    NonEmptyStr,
    TableName,
)
from cdm_data_loaders.core.settings import (
    FileFormatIngestSettings,
    default_settings_with_shortcuts,
)

PIPELINE_NAME: Final[str] = "xmltodict"


class XmlToDictSettings(FileFormatIngestSettings):
    """Settings for the XML ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts={
            **FileFormatIngestSettings.model_config["cli_shortcuts"],
            FILE_GLOB: "g",
        },
    )

    file_glob: Annotated[
        FileGlob,
        Field(
            default=DEFAULT_XML_FILE_GLOB,
            description="Glob pattern for XML files inside each entity's input subdirectory.",
        ),
    ]
    table_name: TableName
    xml_tag: Annotated[
        NonEmptyStr,
        Field(
            description="XML tag to capture the contents of",
        ),
    ]
