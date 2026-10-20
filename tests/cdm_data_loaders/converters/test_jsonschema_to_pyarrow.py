"""Unit tests for jsonschema_to_pyarrow."""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pyarrow as pa
import pytest
import yaml
from pydantic import ValidationError

from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionSpec
from cdm_data_loaders.converters.jsonschema_to_pyarrow import (
    InvalidJSONSchemaError,
    JSONSchemaToPyArrow,
    JSONSchemaToPyArrowError,
)
from cdm_data_loaders.utils.jsonschema.dereferencer import dereference_schema
from tests.cdm_data_loaders.converters.conftest import base_object_schema


def test_jsonschema_to_pyarrow_fail_instance_is_frozen(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """Assigning to a JSONSchemaToPyArrow attribute raises ValidationError."""
    with pytest.raises(ValidationError, match="Instance is frozen"):
        jsonschema_to_pyarrow_converter.treat_unknown_as_string = False


def test_convert_fail_missing_schema_keyword(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """convert() rejects a schema with no top-level '$schema' keyword."""
    with pytest.raises(InvalidJSONSchemaError, match="missing a '\\$schema'"):
        jsonschema_to_pyarrow_converter.convert({"type": "object", "properties": {}})


@pytest.mark.parametrize(
    "root_type",
    ["string", "array", "integer", "number", "boolean", "null", ["object", "null"]],
    ids=["string", "array", "integer", "number", "boolean", "null", "nullable-object"],
)
def test_convert_fail_non_object_root_type(
    jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow, root_type: str | list[str]
) -> None:
    """convert() rejects any root schema whose declared 'type' isn't literally 'object'."""
    with pytest.raises(JSONSchemaToPyArrowError, match="must be of type 'object'"):
        jsonschema_to_pyarrow_converter.convert(base_object_schema(type=root_type))


@pytest.mark.parametrize(
    ("override", "match"),
    [
        pytest.param({"$ref": "#/$defs/Foo"}, "unresolved \\$ref", id="ref"),
        pytest.param({"allOf": [{"type": "string"}]}, "unmerged 'allOf'", id="all-of"),
    ],
)
def test_convert_fail_unresolved_references(
    jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow, override: dict[str, Any], match: str
) -> None:
    """convert() refuses a root schema with un-dereferenced $ref or allOf."""
    schema = base_object_schema(**override)
    del schema["type"]
    with pytest.raises(JSONSchemaToPyArrowError, match=match):
        jsonschema_to_pyarrow_converter.convert(schema)


def test_convert_fail_root_resolves_to_map(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """convert() rejects a root schema that has no 'properties' and resolves to a map."""
    schema = base_object_schema(patternProperties={"^x-": {"type": "string"}})
    del schema["properties"]
    with pytest.raises(JSONSchemaToPyArrowError, match="did not resolve to a struct"):
        jsonschema_to_pyarrow_converter.convert(schema)


def test_convert_fail_strict_unknown_type_is_wrapped() -> None:
    """Emitter errors surface as the facade's own error type."""
    schema = base_object_schema(properties={"value": {"type": ["string", "integer"]}})
    with pytest.raises(JSONSchemaToPyArrowError, match="multi-type union"):
        JSONSchemaToPyArrow(treat_unknown_as_string=False).convert(schema)


def test_convert_pass_simple_object(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """convert() produces a schema with correct nullability and no metadata for a flat object."""
    schema = base_object_schema(properties={"name": {"type": "string"}, "age": {"type": "integer"}}, required=["name"])
    assert jsonschema_to_pyarrow_converter.convert(schema) == pa.schema(
        [pa.field("name", pa.string(), nullable=False), pa.field("age", pa.int64())]
    )


def test_convert_pass_root_without_type_infers_object(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """convert() accepts a root that omits 'type' but declares 'properties'."""
    schema = {"$schema": "https://json-schema.org/draft/2020-12/schema", "properties": {"a": {"type": "string"}}}
    assert jsonschema_to_pyarrow_converter.convert(schema) == pa.schema([pa.field("a", pa.string())])


def test_convert_pass_nested_schema() -> None:
    """Nested objects, arrays, dynamic maps, formats and decimals convert together."""
    schema = base_object_schema(
        properties={
            "id": {"type": "integer", "minimum": 0, "maximum": 10, "title": "Identifier"},
            "created": {"type": "string", "format": "date-time"},
            "price": {"type": "number", "multipleOf": 0.01},
            "tags": {"type": "array", "items": {"type": "string"}},
            "labels": {"type": "object", "additionalProperties": {"type": "string"}},
            "owner": {
                "type": "object",
                "properties": {"name": {"type": "string"}, "active": {"type": "boolean"}},
                "required": ["name"],
            },
        },
        required=["id"],
    )
    schema_out = JSONSchemaToPyArrow().convert(schema)
    assert schema_out == pa.schema(
        [
            pa.field("id", pa.int32(), nullable=False),
            pa.field("created", pa.timestamp("us", tz="UTC")),
            pa.field("price", pa.decimal128(38, 2)),
            pa.field("tags", pa.list_(pa.string())),
            pa.field("labels", pa.map_(pa.string(), pa.string())),
            pa.field(
                "owner",
                pa.struct([pa.field("name", pa.string(), nullable=False), pa.field("active", pa.bool_())]),
            ),
        ]
    )
    assert schema_out.field("id").metadata == {b"jsonschema": b'{"title": "Identifier"}'}


def test_convert_pass_dereferenced_schema() -> None:
    """A dereferenced $ref schema converts using the referenced definition."""
    schema = base_object_schema(
        properties={"owner": {"$ref": "#/$defs/Person"}},
        **{"$defs": {"Person": {"type": "object", "properties": {"name": {"type": "string"}}}}},
    )
    result = JSONSchemaToPyArrow().convert(dereference_schema(schema))
    assert result == pa.schema([pa.field("owner", pa.struct([pa.field("name", pa.string())]))])


def test_convert_to_table_pass_empty_table(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """convert_to_table() returns a zero-row table with the converted schema."""
    schema = base_object_schema(properties={"id": {"type": "integer"}, "tags": {"type": "array"}}, required=["id"])
    table = jsonschema_to_pyarrow_converter.convert_to_table(schema)
    assert table.num_rows == 0
    assert table.schema == jsonschema_to_pyarrow_converter.convert(schema)
    assert table.schema.field("id").nullable is False
    assert table.append_column("extra", pa.array([], pa.string())).column_names == ["id", "tags", "extra"]


def test_convert_to_table_pass_accepts_rows(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """The emitted schema can build a populated table with the matching nested values."""
    schema = base_object_schema(
        properties={"id": {"type": "integer"}, "tags": {"type": "array", "items": {"type": "string"}}},
    )
    arrow_schema = jsonschema_to_pyarrow_converter.convert(schema)
    table = pa.Table.from_pylist([{"id": 1, "tags": ["a", "b"]}, {"id": 2, "tags": None}], schema=arrow_schema)
    assert table.to_pylist() == [{"id": 1, "tags": ["a", "b"]}, {"id": 2, "tags": None}]


def test_convert_from_string_pass_valid_json(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """convert_from_string() parses JSON text and converts it identically to convert()."""
    schema = base_object_schema(properties={"id": {"type": "integer"}})
    result = jsonschema_to_pyarrow_converter.convert_from_string(json.dumps(schema))
    assert result == jsonschema_to_pyarrow_converter.convert(schema)


def test_convert_from_string_fail_invalid_json(jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow) -> None:
    """convert_from_string() propagates a JSONDecodeError for malformed JSON text."""
    with pytest.raises(json.JSONDecodeError):
        jsonschema_to_pyarrow_converter.convert_from_string("{not valid json")


@pytest.mark.parametrize(
    ("suffix", "dumper"),
    [(".json", json.dumps), (".yaml", yaml.safe_dump)],
    ids=["json", "yaml"],
)
def test_convert_from_file_pass_supported_extensions(
    jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow,
    tmp_path: Path,
    suffix: str,
    dumper: Callable[[dict[str, Any]], str],
) -> None:
    """convert_from_file() loads both .json and YAML files."""
    schema = base_object_schema(properties={"id": {"type": "integer"}})
    path = tmp_path / f"schema{suffix}"
    path.write_text(dumper(schema))
    assert jsonschema_to_pyarrow_converter.convert_from_file(str(path)) == jsonschema_to_pyarrow_converter.convert(
        schema
    )


def test_convert_from_file_fail_missing_file(
    jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow, tmp_path: Path
) -> None:
    """convert_from_file() raises FileNotFoundError for a nonexistent path."""
    with pytest.raises(FileNotFoundError):
        jsonschema_to_pyarrow_converter.convert_from_file(str(tmp_path / "does_not_exist.json"))


def test_jsonschema_to_pyarrow_pass_options_reach_emitter() -> None:
    """Format overrides, metadata keywords, union metadata and nesting limits configure the emitter."""
    converter = JSONSchemaToPyArrow(
        format_map={"date-time": pa.string()},
        extra_metadata_keywords=["pattern"],
        emit_unions_as_json=True,
        max_nesting=1,
    )
    schema = base_object_schema(
        properties={
            "when": {"type": "string", "format": "date-time", "pattern": "^2"},
            "either": {"type": ["string", "integer"]},
            "nested": {"type": "object", "properties": {"x": {"type": "integer"}}},
        }
    )
    result = converter.convert(schema)
    assert converter.extra_metadata_keywords == frozenset({"pattern"})
    assert result.field("when").type == pa.string()
    assert result.field("when").metadata == {b"jsonschema": b'{"pattern": "^2"}'}
    assert result.field("either").metadata == {b"original_union": b'["string", "integer"]'}
    assert result.field("nested").type == pa.string()


def test_jsonschema_to_pyarrow_pass_extension_registry() -> None:
    """A custom extension registry validates x- namespaces selected as metadata."""
    registry = DEFAULT_EXTENSIONS.register(
        ExtensionSpec("x-vendor", {"type": "object", "properties": {"level": {"type": "string"}}})
    )
    converter = JSONSchemaToPyArrow(extension_registry=registry, extra_metadata_keywords=frozenset({"x-vendor"}))
    schema = base_object_schema(properties={"name": {"type": "string", "x-vendor": {"level": "public"}}})
    metadata = converter.convert(schema).field("name").metadata
    assert json.loads(metadata[b"jsonschema"]) == {"x-vendor": {"level": "public"}}


def test_jsonschema_to_pyarrow_fail_unregistered_extension(
    jsonschema_to_pyarrow_converter: JSONSchemaToPyArrow,
) -> None:
    """An unregistered x- namespace in the schema is rejected through the facade error."""
    schema = base_object_schema(properties={"name": {"type": "string", "x-vendor": {"level": "public"}}})
    with pytest.raises(JSONSchemaToPyArrowError):
        jsonschema_to_pyarrow_converter.convert(schema)
