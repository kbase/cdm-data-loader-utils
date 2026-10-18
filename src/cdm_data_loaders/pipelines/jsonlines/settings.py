"""Setting for the JSONL validate-and-load pipeline for the KBase CTS."""

from typing import Annotated, Final

from pydantic import Field, field_validator
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import (
    DEFAULT_JSONL_FILE_GLOB,
    DEFAULTS,
    FILE_GLOB,
    LOADER_FILE_FORMAT,
    TABLE_NAME,
    FileGlob,
    LoaderFileFormat,
    NonEmptyStr,
    TableName,
)
from cdm_data_loaders.core.settings import (
    CLI_SHORTCUTS,
    FileFormatIngestSettings,
    default_settings_with_shortcuts,
)

EXTRACT_PIPELINE_NAME: Final[str] = "jsonlines"
PYDANTIC_PIPELINE_NAME: Final[str] = "jsonlines_pydantic_validator"
JSONSCHEMA_PIPELINE_NAME: Final[str] = "jsonlines_jsonschema_validator"

ENTITY_MODELS_MODULE: Final[str] = "entity_models_module"
SCHEMA_FILES_MODULE: Final[str] = "schema_files_module"
TABLE_NAMES: Final[str] = "table_names"

LOCAL_CLI_SHORTCUTS: Final[dict[str, str]] = {
    ENTITY_MODELS_MODULE: "m",
    FILE_GLOB: "g",
    LOADER_FILE_FORMAT: "f",
    SCHEMA_FILES_MODULE: "m",
    TABLE_NAME: "t",
    TABLE_NAMES: "t",
}


class JsonlIngestSettings(FileFormatIngestSettings):
    """Settings for the JSONL ingest pipeline."""

    file_glob: Annotated[
        FileGlob,
        Field(
            default=DEFAULT_JSONL_FILE_GLOB,
            description="Glob pattern for JSONL files",
        ),
    ]
    loader_file_format: LoaderFileFormat = Field(default=DEFAULTS[LOADER_FILE_FORMAT])


class JsonlExtractSettings(JsonlIngestSettings):
    """Settings for the JSONL extraction-only pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=EXTRACT_PIPELINE_NAME,
        cli_shortcuts={
            **CLI_SHORTCUTS,
            **{k: v for k, v in LOCAL_CLI_SHORTCUTS.items() if k in [FILE_GLOB, LOADER_FILE_FORMAT, TABLE_NAME]},
        },
    )

    table_name: TableName


class ValidatedJsonlIngestSettings(JsonlIngestSettings):
    """JSONL ingestion pipeline with validation."""

    table_names: list[NonEmptyStr] | None = Field(
        default=None,
        description="Table names to process. Defaults to every table in the configured registry.",
    )

    @field_validator(TABLE_NAMES, mode="before")
    @classmethod
    def split_table_names(cls, v: str | list[str] | None) -> list[str] | None:
        """Split a comma-separated string into a list. Pass lists and None through unchanged."""
        if v is None or isinstance(v, list):
            return v
        return [name.strip() for name in v.split(",") if name.strip()] or None


class JsonlPydanticIngestSettings(ValidatedJsonlIngestSettings):
    """Pipeline settings for JSONL pipeline that uses Pydantic for validation."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PYDANTIC_PIPELINE_NAME,
        cli_shortcuts={
            **CLI_SHORTCUTS,
            **{
                k: v
                for k, v in LOCAL_CLI_SHORTCUTS.items()
                if k in [ENTITY_MODELS_MODULE, FILE_GLOB, LOADER_FILE_FORMAT, TABLE_NAMES]
            },
        },
    )

    entity_models_module: Annotated[
        str,
        Field(
            description=(
                "Dotted import path to a module that defines `ENTITY_MODELS: dict[str, type[BaseModel]]`. "
                "Each key is a table name and a subdirectory of `input_dir`. "
                "Each value is the Pydantic model for that table."
            ),
        ),
    ]


class JsonlJsonschemaIngestSettings(ValidatedJsonlIngestSettings):
    """Pipeline settings for JSONL ingestion with JSON Schema validation."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=JSONSCHEMA_PIPELINE_NAME,
        cli_shortcuts={
            **CLI_SHORTCUTS,
            **{
                k: v
                for k, v in LOCAL_CLI_SHORTCUTS.items()
                if k
                in [
                    FILE_GLOB,
                    LOADER_FILE_FORMAT,
                    SCHEMA_FILES_MODULE,
                ]
            },
        },
    )

    schema_files_module: Annotated[
        str,
        Field(
            description=(
                "Dotted import path to a module that defines `ENTITY_SCHEMAS: dict[str, dict]`. "
                "Each key is a table name and a subdirectory of `input_dir`. "
                "Each value is a JSON Schema dictionary declaring its draft with `$schema`."
            ),
        ),
    ]
