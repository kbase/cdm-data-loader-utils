"""Shared fixtures for the jsonschema_to_pyspark test suite."""

from typing import Any

import pytest
from jsonschema import Draft7Validator, Draft202012Validator

from cdm_data_loaders.converters.jsonschema_to_pyspark.converter import ConversionContext, JSONSchemaToPySpark


@pytest.fixture
def converter() -> JSONSchemaToPySpark:
    """A JSONSchemaToPySpark instance configured with all default settings."""
    return JSONSchemaToPySpark()


@pytest.fixture
def strict_converter() -> JSONSchemaToPySpark:
    """A JSONSchemaToPySpark instance that raises instead of falling back to StringType."""
    return JSONSchemaToPySpark(treat_unknown_as_string=False)


@pytest.fixture
def ctx(converter: JSONSchemaToPySpark) -> ConversionContext:
    """A default `converter` fixture (Draft 2020-12, no extra metadata keywords)."""
    return converter._build_context(Draft202012Validator)  # noqa: SLF001


@pytest.fixture
def strict_ctx(strict_converter: JSONSchemaToPySpark) -> ConversionContext:
    """Strict ConversionContext fixture: (Draft 2020-12, no extra keywords, treat_unknown_as_string=False.

    For use with the `strict_converter` fixture.
    """
    return strict_converter._build_context(Draft202012Validator)  # noqa: SLF001


@pytest.fixture
def draft7_ctx(converter: JSONSchemaToPySpark) -> ConversionContext:
    """A ConversionContext using the Draft-07 validator class."""
    return converter._build_context(validator_cls=Draft7Validator)  # noqa: SLF001


def base_object_schema(**overrides: Any) -> dict[str, Any]:  # noqa: ANN401
    """Build a minimal valid root object schema, applying keyword overrides."""
    schema: dict[str, Any] = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {},
    }
    schema.update(overrides)
    return schema
