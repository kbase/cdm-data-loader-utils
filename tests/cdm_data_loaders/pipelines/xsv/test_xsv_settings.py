"""Tests for XsvIngestSettings field validation, schema loading, and derived properties."""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any, Final

import pytest
from pydantic import ValidationError

from cdm_data_loaders.core.fields import DEFAULTS
from cdm_data_loaders.pipelines.xsv.settings import PIPELINE_NAME, XsvIngestSettings

VALID_SCHEMA_URI: Final[str] = "https://json-schema.org/draft/2020-12/schema"
COLUMNS: Final[list[str]] = ["number", "date", "float", "boolean", "string"]

VALID_SCHEMA: Final[dict[str, Any]] = {
    "$schema": VALID_SCHEMA_URI,
    "title": "test schema",
    "type": "object",
    "required": COLUMNS,
    "properties": {
        "number": {"type": "integer"},
        "date": {"type": "string"},
        "float": {"type": ["string", "null"]},
        "boolean": {"type": "boolean"},
        "string": {"type": "string"},
    },
}

MINIMAL_SETTINGS_KWARGS: dict[str, object] = {
    "buffer_size": 100,
    "log_interval": 1000,
    "dataset_name": "xsv_test_dataset",
    "loader_file_format": "parquet",
    "table_name": "my_table",
    "log_config_file": None,
    "dlt_dev_mode": False,
    "use_destination": "local_fs",
    "use_output_dir_for_pipeline_metadata": False,
}


@pytest.fixture
def settings_factory(tmp_path: Path) -> Callable[..., XsvIngestSettings]:
    """Return a factory that builds a valid XsvIngestSettings, with fields open to override."""

    def _factory(schema: dict[str, Any] | None = VALID_SCHEMA, **overrides: object) -> XsvIngestSettings:
        input_dir = tmp_path / "input"
        input_dir.mkdir(exist_ok=True)
        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        if schema is not None:
            (input_dir / "schema.json").write_text(json.dumps(schema))

        kwargs: dict[str, object] = {
            **MINIMAL_SETTINGS_KWARGS,
            "input_dir": str(input_dir),
            "output_dir": str(output_dir),
            "file_glob": DEFAULTS["file_glob"],
            "schema_file": str(input_dir / "schema.json"),
        }
        kwargs.update(overrides)
        return XsvIngestSettings(**kwargs)

    return _factory


"""defaults and basic construction"""


def test_xsv_settings_pass_defaults(settings_factory: Callable[..., XsvIngestSettings]) -> None:
    """file_glob, buffer_size, and log_interval default to the common CTS defaults."""
    settings = settings_factory()
    assert settings.file_glob == DEFAULTS["file_glob"]
    assert settings.buffer_size == DEFAULTS["buffer_size"]
    assert settings.log_interval == DEFAULTS["log_interval"]
    assert settings.model_config["cli_prog_name"] == PIPELINE_NAME


def test_xsv_settings_pass_loads_and_validates_schema(
    settings_factory: Callable[..., XsvIngestSettings],
) -> None:
    """validated_schema exposes the parsed schema's required columns."""
    settings = settings_factory()
    assert settings.validated_schema.required_cols == COLUMNS


def test_xsv_settings_fail_missing_schema_file(settings_factory: Callable[..., XsvIngestSettings]) -> None:
    """A schema_file that does not exist on disk raises ValidationError."""
    with pytest.raises(ValidationError, match="Path does not point to a file"):
        settings_factory(schema=None, schema_file="does_not_exist.json")


@pytest.mark.parametrize(
    "invalid_schema",
    [
        pytest.param({"required": COLUMNS}, id="missing-dollar-schema"),
        pytest.param({"$schema": VALID_SCHEMA_URI}, id="missing-required"),
        pytest.param({"$schema": VALID_SCHEMA_URI, "required": []}, id="empty-required"),
    ],
)
def test_xsv_settings_fail_invalid_schema(
    settings_factory: Callable[..., XsvIngestSettings], invalid_schema: dict[str, Any]
) -> None:
    """A structurally invalid schema (missing $schema or required) raises RuntimeError."""
    with pytest.raises(RuntimeError, match="Could not load JSON Schema"):
        settings_factory(schema=invalid_schema)


def test_xsv_settings_fail_empty_schema_file_name(settings_factory: Callable[..., XsvIngestSettings]) -> None:
    """An empty schema_file string is rejected by the NonEmptyStr constraint."""
    with pytest.raises(ValidationError):
        settings_factory(schema_file="")


@pytest.mark.parametrize("field", ["buffer_size", "log_interval"])
def test_xsv_settings_fail_non_positive_int_fields(
    settings_factory: Callable[..., XsvIngestSettings], field: str
) -> None:
    """buffer_size and log_interval must be positive integers."""
    with pytest.raises(ValidationError):
        settings_factory(**{field: 0})


"""derived properties: first_pass_schema and parsing_config"""


def test_xsv_settings_pass_first_pass_schema_is_loose(
    settings_factory: Callable[..., XsvIngestSettings],
) -> None:
    """first_pass_schema derives a loose schema typing every required column as string-or-null."""
    settings = settings_factory()
    first_pass = settings.first_pass_schema
    assert first_pass.required_cols == COLUMNS
    assert all(first_pass.jsonschema["properties"][col]["type"] == ["string", "null"] for col in COLUMNS)


def test_xsv_settings_pass_first_pass_schema_is_cached(
    settings_factory: Callable[..., XsvIngestSettings],
) -> None:
    """first_pass_schema is computed once and reused on subsequent access."""
    settings = settings_factory()
    assert settings.first_pass_schema is settings.first_pass_schema


def test_xsv_settings_pass_parsing_config_empty_without_xsv_config_block(
    settings_factory: Callable[..., XsvIngestSettings],
) -> None:
    """parsing_config is empty when the schema carries no x-xsv-config block."""
    settings = settings_factory()
    assert settings.parsing_config == {}


def test_xsv_settings_pass_parsing_config_resolved_from_xsv_config_block(
    settings_factory: Callable[..., XsvIngestSettings],
) -> None:
    """parsing_config resolves qsv kwargs from the schema's x-xsv-config block."""
    schema_with_config = {**VALID_SCHEMA, "x-xsv-config": {"x-delimiter": ",", "x-has-header": True}}
    settings = settings_factory(schema=schema_with_config)
    assert settings.parsing_config == {"delimiter": ",", "missing_header": False}


def test_xsv_settings_pass_parsing_config_is_cached(settings_factory: Callable[..., XsvIngestSettings]) -> None:
    """parsing_config is computed once and reused on subsequent access."""
    settings = settings_factory()
    assert settings.parsing_config is settings.parsing_config


def test_xsv_settings_fail_validated_schema_unset_raises() -> None:
    """Accessing validated_schema before it has been set raises RuntimeError."""
    settings = XsvIngestSettings.model_construct()
    with pytest.raises(RuntimeError, match="Schema has not been validated"):
        _ = settings.validated_schema
