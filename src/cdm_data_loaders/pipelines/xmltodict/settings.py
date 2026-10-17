"""Settings class for the XML to dictionary pipeline for the KBase CTS."""

from functools import cached_property
from logging import Logger, getLogger
from pathlib import Path
from typing import Annotated, Any, Final, Self

from pydantic import Field, PrivateAttr, computed_field, model_validator
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import (
    DEFAULT_XML_FILE_GLOB,
    FILE_GLOB,
    BufferSize,
    DatasetName,
    FileGlob,
    LoaderFileFormat,
    LogInterval,
    NonEmptyStr,
    TableName,
)
from cdm_data_loaders.core.settings import (
    CLI_SHORTCUTS,
    CtsSettings,
    default_settings_with_shortcuts,
)
from cdm_data_loaders.readers.xsd import find_list_and_single_child_paths, load_schema

PIPELINE_NAME: Final[str] = "xmltodict"

logger: Logger = getLogger(__name__)


class XmlToDictSettings(CtsSettings):
    """Settings for the XML ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts={
            **CLI_SHORTCUTS,
            FILE_GLOB: "g",
        },
    )

    _paths: Annotated[dict[str, Any], Field(description="parent-child paths in the schema")] = PrivateAttr(
        default_factory=dict
    )

    buffer_size: BufferSize
    dataset_name: DatasetName
    file_glob: Annotated[
        FileGlob,
        Field(
            default=DEFAULT_XML_FILE_GLOB,
            description="Glob pattern for XML files inside each entity's input subdirectory.",
        ),
    ]
    loader_file_format: LoaderFileFormat
    log_interval: LogInterval
    table_name: TableName
    xml_tag: Annotated[
        NonEmptyStr,
        Field(
            description="XML tag to capture the contents of",
        ),
    ]
    xsd_file: Annotated[
        NonEmptyStr | None,
        Field(
            None,
            description="XSD schema file for the dataset; this should be a relative path from the input_dir. Including the schema will ensure that the appropriate relationships (1-to-1 vs 1-to-many) are created between parent and child elements.",
        ),
    ]

    @model_validator(mode="after")
    def check_xsd_file(self) -> Self:
        """Check whether the schema file is present and correct; if so, init self._paths."""
        if self.xsd_file:
            try:
                schema_path = Path(self.input_dir) / self.xsd_file
                schema = load_schema(schema_path)
                list_paths, single_paths = find_list_and_single_child_paths(schema)
                self._paths = {"list": list_paths, "single": single_paths}
            except Exception as e:
                err_msg = f"Could not generate parent-child relationships for schema: {e!s}"
                logger.exception(err_msg)
                raise RuntimeError(err_msg) from e

        return self

    @computed_field
    @cached_property
    def xmltodict_args(self) -> dict[str, Any]:
        """Extra parser configuration arguments for running xmltodict over a dataset.

        :return: kwargs for configuring xmltodict.parse()
        :rtype: dict[str, Any]
        """
        if not self._paths:
            return {}

        single_paths = self._paths["single"]

        def force_list(path: list[tuple[str, str | None]], key: str, _: Any) -> bool:
            if not path:
                return False

            return (path[-1][0], key) not in single_paths

        return {"force_list": force_list}
