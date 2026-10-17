"""Settings for the XSV (CSV/TSV) ingestion pipeline for the KBase CTS."""

from functools import cached_property
from logging import Logger, getLogger
from pathlib import Path
from typing import Annotated, Any, Final, Self

from pydantic import Field, PrivateAttr, model_validator
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import (
    BufferSize,
    DatasetName,
    FileGlob,
    LoaderFileFormat,
    LogInterval,
    NonEmptyStr,
    TableName,
)
from cdm_data_loaders.core.settings import CLI_SHORTCUTS, CtsSettings, default_settings_with_shortcuts
from cdm_data_loaders.readers.jsonschema_xsv.xsv_reader import resolve_xsv_parsing_config
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import (
    ValidatedSchema,
    generate_first_pass_schema,
    validate_jsonschema,
)

PIPELINE_NAME: Final[str] = "xsv"

logger: Logger = getLogger(__name__)


SHORTCUTS: Final[dict[str, str]] = {
    **CLI_SHORTCUTS,
    "file_glob": "g",
    "loader_file_format": "f",
    "schema_file": "s",
    "table_name": "t",
}


class XsvIngestSettings(CtsSettings):
    """Settings for the XSV (CSV/TSV) ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts=SHORTCUTS,
    )

    _validated_schema: ValidatedSchema | None = PrivateAttr(default=None)

    buffer_size: BufferSize
    dataset_name: DatasetName
    file_glob: FileGlob
    loader_file_format: LoaderFileFormat
    log_interval: LogInterval
    schema_file: Annotated[
        NonEmptyStr,
        Field(
            description=(
                "JSON Schema file describing the XSV data's columns and types; relative path "
                "from input_dir. May include an x-xsv-config block describing the file's "
                "delimiter, comment character, quoting, and null value handling."
            ),
        ),
    ]
    table_name: TableName

    @model_validator(mode="after")
    def check_schema_file(self) -> Self:
        """Load and validate the schema file."""
        try:
            self._validated_schema = validate_jsonschema(Path(self.input_dir) / self.schema_file)
        except Exception as e:
            err_msg = f"Could not load JSON Schema for XSV data: {e!s}"
            logger.exception(err_msg)
            raise RuntimeError(err_msg) from e
        return self

    @property
    def validated_schema(self) -> ValidatedSchema:
        """The parsed, validated JSON Schema describing the XSV data."""
        if self._validated_schema is None:
            err_msg = "Schema has not been validated"
            raise RuntimeError(err_msg)
        return self._validated_schema

    @cached_property
    def first_pass_schema(self) -> ValidatedSchema:
        """A loose schema used to validate only that each row has the correct number of columns."""
        return generate_first_pass_schema(self.validated_schema)

    @cached_property
    def parsing_config(self) -> dict[str, Any]:
        """Cleaning/validation parameters for qsv, derived from the schema's x-xsv-config block."""
        return resolve_xsv_parsing_config(self.validated_schema)
