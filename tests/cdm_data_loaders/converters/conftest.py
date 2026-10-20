"""Shared fixtures and helpers for the converters test suite."""

import re
from collections.abc import Iterator
from typing import Any, Final

import pytest
from jsonschema import Draft7Validator, Draft202012Validator
from pyiceberg.partitioning import PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.table import Table
from pyiceberg.table.metadata import new_table_metadata
from pyiceberg.table.sorting import SortOrder
from pyiceberg.types import LongType, NestedField, StringType

from cdm_data_loaders.converters.dlt_to_jsonschema import DltToJSONSchema
from cdm_data_loaders.converters.jsonschema_to_dlt import JSONSchemaToDlt
from cdm_data_loaders.converters.jsonschema_to_pyarrow import JSONSchemaToPyArrow
from cdm_data_loaders.converters.jsonschema_to_pyspark import ConversionContext, JSONSchemaToPySpark

JSON_SCHEMA_KEYWORDS: Final = frozenset(
    {
        "$schema",
        "$id",
        "$defs",
        "title",
        "description",
        "type",
        "format",
        "pattern",
        "items",
        "properties",
        "required",
        "additionalProperties",
        "anyOf",
        "minimum",
        "maximum",
        "contentEncoding",
    }
)


@pytest.fixture
def dlt_to_jsonschema_converter() -> DltToJSONSchema:
    """A DltToJSONSchema instance configured with all default settings."""
    return DltToJSONSchema()


@pytest.fixture
def jsonschema_to_dlt_converter() -> JSONSchemaToDlt:
    """A JSONSchemaToDlt instance configured with all default settings."""
    return JSONSchemaToDlt(schema_name="test_source")


@pytest.fixture
def nested_json_converter() -> JSONSchemaToDlt:
    """A JSONSchemaToDlt instance that keeps nested structures as json columns."""
    return JSONSchemaToDlt(schema_name="test_source", skip_nested_types=True)


@pytest.fixture
def jsonschema_to_pyspark_converter() -> JSONSchemaToPySpark:
    """A JSONSchemaToPySpark instance configured with all default settings."""
    return JSONSchemaToPySpark()


@pytest.fixture
def jsonschema_to_pyarrow_converter() -> JSONSchemaToPyArrow:
    """A JSONSchemaToPyArrow instance configured with all default settings."""
    return JSONSchemaToPyArrow()


@pytest.fixture
def strict_converter() -> JSONSchemaToPySpark:
    """A JSONSchemaToPySpark instance that raises instead of falling back to StringType."""
    return JSONSchemaToPySpark(treat_unknown_as_string=False)


@pytest.fixture
def ctx(jsonschema_to_pyspark_converter: JSONSchemaToPySpark) -> ConversionContext:
    """A default `jsonschema_to_pyspark_converter` fixture (Draft 2020-12, no extra metadata keywords)."""
    return jsonschema_to_pyspark_converter._build_context(Draft202012Validator)  # noqa: SLF001


@pytest.fixture
def strict_ctx(strict_converter: JSONSchemaToPySpark) -> ConversionContext:
    """Strict ConversionContext fixture: (Draft 2020-12, no extra keywords, treat_unknown_as_string=False.

    For use with the `strict_converter` fixture.
    """
    return strict_converter._build_context(Draft202012Validator)  # noqa: SLF001


@pytest.fixture
def draft7_ctx(jsonschema_to_pyspark_converter: JSONSchemaToPySpark) -> ConversionContext:
    """A ConversionContext using the Draft-07 validator class."""
    return jsonschema_to_pyspark_converter._build_context(validator_cls=Draft7Validator)  # noqa: SLF001


def base_stored_schema(tables: dict[str, Any]) -> dict[str, Any]:
    """Build a minimal stored-schema dict wrapping the given tables."""
    return {"name": "test_source", "tables": tables}


def base_table(**overrides: Any) -> dict[str, Any]:  # noqa: ANN401
    """Build a minimal TTableSchema with a single non-nullable text column, applying overrides."""
    table: dict[str, Any] = {
        "columns": {"col": {"name": "col", "data_type": "text", "nullable": False}},
    }
    table.update(overrides)
    return table


def base_object_schema(**overrides: Any) -> dict[str, Any]:  # noqa: ANN401
    """Build a minimal valid root object schema, applying keyword overrides."""
    schema: dict[str, Any] = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {},
    }
    schema.update(overrides)
    return schema


def make_table(
    identifier: tuple[str, ...],
    schema: Schema,
    partition_spec: PartitionSpec | None = None,
    properties: dict[str, str] | None = None,
) -> Table:
    """Build a real pyiceberg Table backed by freshly generated metadata."""
    metadata = new_table_metadata(
        schema,
        partition_spec if partition_spec is not None else PartitionSpec(),
        SortOrder(),
        f"file:///tmp/{'/'.join(identifier)}",
        properties or {},
    )
    return Table(identifier, metadata, f"{metadata.location}/metadata/00000.metadata.json", None, None)


@pytest.fixture
def simple_schema() -> Schema:
    """A schema with one required and one optional field."""
    return Schema(
        NestedField(1, "id", LongType(), required=True, doc="primary key"),
        NestedField(2, "name", StringType(), doc="display name"),
    )


def iter_schema_keywords(node: dict) -> Iterator[str]:
    """Yield every schema keyword key in a schema document, descending into subschemas."""
    for key, value in node.items():
        yield key
        if key in {"properties", "$defs"} and isinstance(value, dict):
            for subschema in value.values():
                yield from iter_schema_keywords(subschema)
        elif key == "items" and isinstance(value, dict):
            yield from iter_schema_keywords(value)
        elif key == "anyOf" and isinstance(value, list):
            for subschema in value:
                if isinstance(subschema, dict):
                    yield from iter_schema_keywords(subschema)


def iter_property_names(node: dict) -> Iterator[str]:
    """Yield every property name in a schema document, descending into subschemas."""
    for key, value in node.items():
        if key == "properties" and isinstance(value, dict):
            for name, subschema in value.items():
                yield name
                yield from iter_property_names(subschema)
        elif key == "items" and isinstance(value, dict):
            yield from iter_property_names(value)
        elif key == "anyOf" and isinstance(value, list):
            for subschema in value:
                if isinstance(subschema, dict):
                    yield from iter_property_names(subschema)


def iter_extension_keys(node: dict) -> Iterator[str]:
    """Yield every x- prefixed key and every key inside extension blocks."""
    for key, value in node.items():
        if key.startswith("x-"):
            yield key
            if isinstance(value, dict):
                yield from iter_extension_keys(value)
        elif isinstance(value, dict):
            yield from iter_extension_keys(value)
        elif isinstance(value, list):
            for item in value:
                if isinstance(item, dict):
                    yield from iter_extension_keys(item)


def to_snake_case(key: str) -> str:
    """Return the key unchanged if already snake case."""
    return re.sub(r"(?<!^)(?=[A-Z])", "_", key).lower()
