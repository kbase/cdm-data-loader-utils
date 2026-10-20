"""Tests for XsvIngestSettings field validation, schema loading, and derived properties."""

from collections.abc import Callable
from typing import Any

import pytest
from pydantic import ValidationError

from cdm_data_loaders.core.fields import DEFAULTS
from cdm_data_loaders.pipelines.xsv.settings import PIPELINE_NAME, XsvIngestSettings
from tests.cdm_data_loaders.pipelines.xsv.conftest import COLUMNS, DEFAULT_XSV_SCHEMA, VALID_SCHEMA_URI

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
    schema_without_config = {k: v for k, v in DEFAULT_XSV_SCHEMA.items() if k != "x-xsv-config"}
    settings = settings_factory(schema=schema_without_config)
    assert settings.parsing_config == {}


def test_xsv_settings_pass_parsing_config_resolved_from_xsv_config_block(
    settings_factory: Callable[..., XsvIngestSettings],
) -> None:
    """parsing_config resolves qsv kwargs from the schema's x-xsv-config block."""
    settings = settings_factory()
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
