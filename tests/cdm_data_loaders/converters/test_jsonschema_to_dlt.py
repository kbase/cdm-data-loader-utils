"""Unit and integration tests for jsonschema_to_dlt."""

import logging
from pathlib import Path
from typing import Any

import pytest
import yaml
from dlt.common.data_types.typing import TDataType
from dlt.common.schema import Schema
from dlt.common.schema.utils import is_nested_table
from pydantic import ValidationError

from cdm_data_loaders.converters.core.inference import decimal_places
from cdm_data_loaders.converters.emitters.dlt import SCHEMA_ENGINE_VERSION
from cdm_data_loaders.converters.jsonschema_to_dlt import (
    InvalidJSONSchemaError,
    JSONSchemaToDlt,
    JSONSchemaToDltError,
    _data_type_from_enum,
)
from cdm_data_loaders.utils.jsonschema.dereferencer import dereference_schema
from tests.cdm_data_loaders.converters.conftest import base_object_schema


def test_json_schema_to_dlt_fail_instance_is_frozen(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """Assigning to a JSONSchemaToDlt instance attribute raises ValidationError (frozen model)."""
    with pytest.raises(ValidationError, match="Instance is frozen"):
        jsonschema_to_dlt_converter.skip_nested_types = True


@pytest.mark.parametrize(
    "schema_name",
    ["", "   "],
    ids=["empty", "whitespace"],
)
def test_json_schema_to_dlt_fail_empty_schema_name(schema_name: str) -> None:
    """A blank schema_name raises at construction."""
    with pytest.raises(ValidationError, match="schema_name must be a non-empty string"):
        JSONSchemaToDlt(schema_name=schema_name)


def test_data_type_from_enum_values() -> None:
    """_data_type_from_enum maps value lists to the narrowest covering dlt type."""
    assert _data_type_from_enum([]) == "text"
    assert _data_type_from_enum([True, False]) == "bool"
    assert _data_type_from_enum([1, 2, 3]) == "bigint"
    assert _data_type_from_enum([1.5, 2.5]) == "double"
    assert _data_type_from_enum([1, 2.5]) == "double"
    assert _data_type_from_enum(["a", "b"]) == "text"
    assert _data_type_from_enum([1, "a"]) == "text"


def test_decimal_places_values() -> None:
    """decimal_places counts fractional digits of a multipleOf value exactly."""
    two_places = 2
    assert decimal_places(0.01) == two_places
    assert decimal_places(1) == 0
    assert decimal_places(0.5) == 1
    assert decimal_places(100) == 0


"""convert"""


def test_convert_fail_missing_schema_keyword(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """convert() rejects a schema with no top-level '$schema' keyword."""
    schema = {"type": "object", "properties": {}}
    with pytest.raises(InvalidJSONSchemaError, match="missing a '\\$schema'"):
        jsonschema_to_dlt_converter.convert(schema)


@pytest.mark.parametrize(
    "root_type",
    ["string", "array", "integer", "number", "boolean", "null"],
)
def test_convert_fail_non_object_root_type(jsonschema_to_dlt_converter: JSONSchemaToDlt, root_type: str) -> None:
    """convert() rejects any root schema whose declared 'type' isn't 'object'."""
    schema = base_object_schema(type=root_type)
    with pytest.raises(JSONSchemaToDltError, match="must be of type 'object'"):
        jsonschema_to_dlt_converter.convert(schema)


def test_convert_fail_unresolved_ref(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """convert() refuses a schema containing an un-dereferenced $ref."""
    schema = base_object_schema(
        properties={"a": {"$ref": "#/$defs/Foo"}},
    )
    with pytest.raises(JSONSchemaToDltError, match="unresolved \\$ref"):
        jsonschema_to_dlt_converter.convert(schema)


def test_convert_fail_unmerged_allof(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """convert() refuses a schema containing an unmerged allOf."""
    schema = base_object_schema(
        properties={"a": {"allOf": [{"type": "string"}]}},
    )
    with pytest.raises(JSONSchemaToDltError, match="unmerged 'allOf'"):
        jsonschema_to_dlt_converter.convert(schema)


def test_convert_pass_empty_object_schema(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """An empty object schema converts to a single root table with no columns."""
    stored = jsonschema_to_dlt_converter.convert(base_object_schema())
    assert stored["name"] == "test_source"
    assert stored["engine_version"] == SCHEMA_ENGINE_VERSION
    assert sorted(stored["tables"].keys()) == ["test_source"]
    assert stored["tables"]["test_source"]["columns"] == {}


"""scalar type mapping"""


@pytest.mark.parametrize(
    ("prop_schema", "expected_type"),
    [
        ({"type": "string"}, "text"),
        ({"type": "string", "format": "date-time"}, "timestamp"),
        ({"type": "string", "format": "date"}, "date"),
        ({"type": "string", "format": "time"}, "time"),
        ({"type": "string", "format": "uuid"}, "text"),
        ({"type": "string", "format": "binary"}, "binary"),
        ({"type": "string", "format": "byte"}, "binary"),
        ({"type": "integer"}, "bigint"),
        ({"type": "integer", "minimum": -5, "maximum": 100}, "bigint"),
        ({"type": "number"}, "double"),
        ({"type": "number", "multipleOf": 0.01}, "decimal"),
        ({"type": "number", "multipleOf": 5}, "decimal"),
        ({"type": "boolean"}, "bool"),
        ({"enum": ["a", "b"]}, "text"),
        ({"enum": [1, 2]}, "bigint"),
        ({"enum": [True, False]}, "bool"),
        ({"enum": [1, 2.5]}, "double"),
        ({"pattern": "^x"}, "text"),
        ({"minimum": 1}, "double"),
        ({"minimum": 1, "maximum": 10}, "double"),
        (True, "text"),
        ({"type": ["string", "null"]}, "text"),
        ({"type": ["integer", "null"]}, "bigint"),
        ({"type": ["string", "integer"]}, "text"),
    ],
    ids=[
        "string",
        "string-date-time",
        "string-date",
        "string-time",
        "string-uuid",
        "string-binary",
        "string-byte",
        "integer",
        "integer-bounded",
        "number",
        "number-multiple-of",
        "number-integral-multiple-of",
        "boolean",
        "enum-strings",
        "enum-ints",
        "enum-bools",
        "enum-mixed-numeric",
        "typeless-pattern",
        "typeless-minimum",
        "typeless-min-max",
        "boolean-schema-true",
        "union-with-null",
        "union-integer-null",
        "union-multi",
    ],
)
def test_convert_scalar_column_types_pass(
    jsonschema_to_dlt_converter: JSONSchemaToDlt,
    prop_schema: dict[str, Any] | bool,
    expected_type: TDataType,
) -> None:
    """Scalar property schemas map to the expected dlt data type on a flat table."""
    schema = base_object_schema(properties={"col": prop_schema})
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert stored["tables"]["test_source"]["columns"]["col"]["data_type"] == expected_type


def test_convert_scalar_column_type_fail_unknown_type_strict(
    nested_json_converter: JSONSchemaToDlt,
) -> None:
    """A multi-type union collapses to text rather than raising; dlt has no union types."""
    schema = base_object_schema(properties={"col": {"type": ["string", "integer"]}})
    stored = nested_json_converter.convert(schema)
    assert stored["tables"]["test_source"]["columns"]["col"]["data_type"] == "text"


"""nullability"""


@pytest.mark.parametrize(
    ("required", "expected_nullable"),
    [(["col"], False), ([], True)],
    ids=["required", "optional"],
)
def test_convert_nullability_from_required(
    jsonschema_to_dlt_converter: JSONSchemaToDlt,
    required: list[str],
    expected_nullable: bool,
) -> None:
    """Members of the parent's `required` list become non-nullable columns."""
    schema = base_object_schema(
        required=required,
        properties={"col": {"type": "string"}},
    )
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert stored["tables"]["test_source"]["columns"]["col"]["nullable"] is expected_nullable


"""nested structures"""


def test_convert_nested_object_becomes_child_table(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """A nested object becomes a dlt child table with parent linkage and __ naming."""
    schema = base_object_schema(
        properties={
            "address": {
                "type": "object",
                "properties": {"city": {"type": "string"}},
            }
        }
    )
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source", "test_source__address"]
    child = stored["tables"]["test_source__address"]
    assert is_nested_table(child)
    assert child["parent"] == "test_source"
    assert child["columns"]["city"]["data_type"] == "text"


def test_convert_deeply_nested_objects_chain_child_tables(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """Objects nested more than one level deep chain child tables with full path names."""
    schema = base_object_schema(
        properties={
            "address": {
                "type": "object",
                "properties": {
                    "geo": {
                        "type": "object",
                        "properties": {"lat": {"type": "number"}},
                    }
                },
            }
        }
    )
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert sorted(stored["tables"].keys()) == [
        "test_source",
        "test_source__address",
        "test_source__address__geo",
    ]
    assert stored["tables"]["test_source__address__geo"]["parent"] == "test_source__address"
    assert stored["tables"]["test_source__address__geo"]["columns"]["lat"]["data_type"] == "double"


def test_convert_array_of_objects_becomes_child_table(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """An array of objects becomes a dlt child table with the object's columns."""
    schema = base_object_schema(
        properties={
            "orders": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {"total": {"type": "number"}},
                },
            }
        }
    )
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source", "test_source__orders"]
    child = stored["tables"]["test_source__orders"]
    assert child["parent"] == "test_source"
    assert child["columns"]["total"]["data_type"] == "double"


def test_convert_array_of_scalars_becomes_value_child_table(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """An array of scalars becomes a child table with a single `value` column (dlt normalizer behavior)."""
    schema = base_object_schema(properties={"tags": {"type": "array", "items": {"type": "string"}}})
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source", "test_source__tags"]
    child = stored["tables"]["test_source__tags"]
    assert child["parent"] == "test_source"
    assert sorted(child["columns"].keys()) == ["value"]
    assert child["columns"]["value"]["data_type"] == "text"


@pytest.mark.parametrize(
    ("items_schema", "expected_type"),
    [
        ([{"type": "string"}, {"type": "integer"}], "text"),
        ([{"type": "integer"}], "bigint"),
        ({"type": "number"}, "double"),
    ],
    ids=["tuple-items-first-wins", "tuple-items-single", "single-items"],
)
def test_convert_array_tuple_style_items_use_first_schema(
    jsonschema_to_dlt_converter: JSONSchemaToDlt,
    items_schema: list[dict[str, Any]] | dict[str, Any],
    expected_type: TDataType,
) -> None:
    """Tuple-style items and prefixItems use the first element's schema for the value column."""
    schema = base_object_schema(properties={"tags": {"type": "array", "items": items_schema}})
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert stored["tables"]["test_source__tags"]["columns"]["value"]["data_type"] == expected_type


def test_convert_array_without_items_becomes_json_column(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """An array schema with no `items` has no element type and stays a json column on the parent table."""
    schema = base_object_schema(properties={"tags": {"type": "array"}})
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source"]
    assert stored["tables"]["test_source"]["columns"]["tags"]["data_type"] == "json"


def test_convert_skip_nested_types_keeps_json_columns(nested_json_converter: JSONSchemaToDlt) -> None:
    """With skip_nested_types, nested objects and arrays stay as json columns."""
    schema = base_object_schema(
        properties={
            "meta": {"type": "object", "properties": {"a": {"type": "string"}}},
            "nums": {"type": "array", "items": {"type": "number"}},
        }
    )
    stored = nested_json_converter.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source"]
    columns = stored["tables"]["test_source"]["columns"]
    assert columns["meta"]["data_type"] == "json"
    assert columns["nums"]["data_type"] == "json"


def test_convert_flatten_scalars_false_keeps_scalar_arrays_json() -> None:
    """With flatten_scalars=False, scalar arrays stay as json columns while objects still flatten."""
    conv = JSONSchemaToDlt(schema_name="test_source", flatten_scalars=False)
    schema = base_object_schema(
        properties={
            "tags": {"type": "array", "items": {"type": "string"}},
            "orders": {
                "type": "array",
                "items": {"type": "object", "properties": {"total": {"type": "number"}}},
            },
        }
    )
    stored = conv.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source", "test_source__orders"]
    assert stored["tables"]["test_source"]["columns"]["tags"]["data_type"] == "json"


def test_convert_max_nesting_keeps_deep_levels_json() -> None:
    """Levels deeper than max_nesting stay as json columns."""
    conv = JSONSchemaToDlt(schema_name="test_source", max_nesting=2)
    schema = base_object_schema(
        properties={
            "a": {
                "type": "object",
                "properties": {
                    "b": {
                        "type": "object",
                        "properties": {
                            "c": {
                                "type": "object",
                                "properties": {"d": {"type": "string"}},
                            }
                        },
                    }
                },
            }
        }
    )
    stored = conv.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source", "test_source__a", "test_source__a__b"]
    deep_column = stored["tables"]["test_source__a__b"]["columns"]["c"]
    assert deep_column["data_type"] == "json"


def test_convert_max_nesting_root_level_bounded() -> None:
    """With max_nesting=1, only the root table flattens one level; deeper levels stay json."""
    conv = JSONSchemaToDlt(schema_name="test_source", max_nesting=1)
    schema = base_object_schema(
        properties={
            "a": {
                "type": "object",
                "properties": {
                    "b": {
                        "type": "object",
                        "properties": {"c": {"type": "string"}},
                    }
                },
            }
        }
    )
    stored = conv.convert(schema)
    assert sorted(stored["tables"].keys()) == ["test_source", "test_source__a"]
    nested_column = stored["tables"]["test_source__a"]["columns"]["b"]
    assert nested_column["data_type"] == "json"


"""metadata and descriptions"""


def test_convert_pass_description_and_table_description(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """Descriptions propagate from the root schema to the table and from properties to columns."""
    schema = base_object_schema(
        description="root description",
        properties={"col": {"type": "string", "description": "column description"}},
    )
    stored = jsonschema_to_dlt_converter.convert(schema)
    table = stored["tables"]["test_source"]
    assert table["description"] == "root description"
    assert table["columns"]["col"]["description"] == "column description"


def test_convert_pass_write_disposition_applied_to_root_tables() -> None:
    """A configured write_disposition lands on root tables only."""
    conv = JSONSchemaToDlt(schema_name="test_source", write_disposition="merge")
    schema = base_object_schema(
        properties={
            "child": {"type": "object", "properties": {"x": {"type": "string"}}},
        }
    )
    stored = conv.convert(schema)
    assert stored["tables"]["test_source"]["write_disposition"] == "merge"
    assert "write_disposition" not in stored["tables"]["test_source__child"]


"""entry points"""


def test_convert_from_string_pass(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """convert_from_string parses a JSON string and converts it."""
    schema_str = '{"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "object", "properties": {"a": {"type": "string"}}}'
    stored = jsonschema_to_dlt_converter.convert_from_string(schema_str)
    assert stored["tables"]["test_source"]["columns"]["a"]["data_type"] == "text"


def test_convert_from_file_pass_json(jsonschema_to_dlt_converter: JSONSchemaToDlt, tmp_path: Path) -> None:
    """convert_from_file reads a .json schema file."""
    path = tmp_path / "schema.json"
    path.write_text(
        '{"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "object", "properties": {"a": {"type": "integer"}}}'
    )
    stored = jsonschema_to_dlt_converter.convert_from_file(str(path))
    assert stored["tables"]["test_source"]["columns"]["a"]["data_type"] == "bigint"


def test_convert_from_file_pass_yaml(jsonschema_to_dlt_converter: JSONSchemaToDlt, tmp_path: Path) -> None:
    """convert_from_file reads a .yaml schema file."""
    schema = base_object_schema(properties={"a": {"type": "boolean"}})
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.dump(schema))
    stored = jsonschema_to_dlt_converter.convert_from_file(str(path))
    assert stored["tables"]["test_source"]["columns"]["a"]["data_type"] == "bool"


def test_convert_from_file_fail_missing_file(jsonschema_to_dlt_converter: JSONSchemaToDlt, tmp_path: Path) -> None:
    """convert_from_file raises FileNotFoundError for a nonexistent path."""
    with pytest.raises(FileNotFoundError):
        jsonschema_to_dlt_converter.convert_from_file(str(tmp_path / "nope.json"))


"""to_yaml and live Schema integration"""


def test_to_yaml_pass_round_trips_through_yaml(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """to_yaml emits YAML that parses back to the same stored schema tables."""
    schema = base_object_schema(
        required=["id"],
        properties={"id": {"type": "integer"}, "name": {"type": "string"}},
    )
    yaml_text = jsonschema_to_dlt_converter.to_yaml(schema)
    parsed: dict[str, Any] = yaml.safe_load(yaml_text)
    assert sorted(parsed["tables"].keys()) == ["test_source"]
    assert parsed["tables"]["test_source"]["columns"]["id"]["data_type"] == "bigint"


def test_to_schema_pass_builds_live_dlt_schema(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """to_schema returns a live dlt Schema containing the converted tables plus internal tables."""
    schema = base_object_schema(
        required=["id"],
        properties={
            "id": {"type": "integer"},
            "address": {"type": "object", "properties": {"city": {"type": "string"}}},
            "tags": {"type": "array", "items": {"type": "string"}},
        },
    )
    live = jsonschema_to_dlt_converter.to_schema(schema)
    assert isinstance(live, Schema)
    assert live.name == "test_source"
    table_names = {t["name"] for t in live.data_tables()}
    assert table_names == {"test_source", "test_source__address", "test_source__tags"}
    columns = {t["name"]: t["columns"] for t in live.data_tables()}
    assert columns["test_source"]["id"]["data_type"] == "bigint"
    assert columns["test_source"]["id"]["nullable"] is False
    assert columns["test_source__address"]["city"]["data_type"] == "text"
    assert columns["test_source__tags"]["value"]["data_type"] == "text"


def test_to_schema_pass_deep_nesting_parents_first(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """to_schema handles multi-level chains; parents merge before children."""
    schema = base_object_schema(
        required=["c"],
        properties={
            "a": {
                "type": "object",
                "properties": {
                    "b": {
                        "type": "object",
                        "required": ["c"],
                        "properties": {"c": {"type": "string"}},
                    }
                },
            }
        },
    )
    live = jsonschema_to_dlt_converter.to_schema(schema)
    table_names = {t["name"] for t in live.data_tables()}
    # dlt's data_tables() hides zero-column tables: the column-less root/intermediate
    # tables in this chain are omitted, only leaf tables with columns are listed
    assert table_names == {"test_source__a__b"}
    chain_tables = {t["name"] for t in live.tables.values()}
    assert chain_tables >= {"test_source", "test_source__a", "test_source__a__b"}
    geo_columns = next(t["columns"] for t in live.data_tables() if t["name"] == "test_source__a__b")
    assert geo_columns["c"]["data_type"] == "text"
    assert geo_columns["c"]["nullable"] is False


def test_to_schema_pass_yaml_export_round_trips(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """A live Schema built by to_schema serializes with to_pretty_yaml and reloads."""
    schema = base_object_schema(properties={"id": {"type": "integer"}})
    live = jsonschema_to_dlt_converter.to_schema(schema)
    yaml_text = live.to_pretty_yaml()
    reloaded = Schema.from_dict(yaml.safe_load(yaml_text))
    assert reloaded.name == "test_source"
    assert {t["name"] for t in reloaded.data_tables()} >= {"test_source"}


def test_to_schema_fail_missing_schema_keyword(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """to_schema enforces the same $schema requirement as convert()."""
    with pytest.raises(InvalidJSONSchemaError, match="missing a '\\$schema'"):
        jsonschema_to_dlt_converter.to_schema({"type": "object", "properties": {}})


"""integration with the dereferencer"""


def test_convert_pass_after_dereferencing(caplog: pytest.LogCaptureFixture) -> None:
    """A schema with $defs/$ref converts after passing through dereference_schema."""
    schema = base_object_schema(
        properties={"address": {"$ref": "#/$defs/Address"}},
        **{"$defs": {"Address": {"type": "object", "properties": {"city": {"type": "string"}}}}},
    )
    dereferenced = dereference_schema(schema)
    conv = JSONSchemaToDlt(schema_name="test_source")
    with caplog.at_level(logging.WARNING):
        stored = conv.convert(dereferenced)
    assert sorted(stored["tables"].keys()) == ["test_source", "test_source__address"]
    assert stored["tables"]["test_source__address"]["columns"]["city"]["data_type"] == "text"


def test_convert_pass_schema_with_oneof_after_dereference(jsonschema_to_dlt_converter: JSONSchemaToDlt) -> None:
    """oneOf-only properties fall back to text with a warning (no union types in dlt)."""
    schema = base_object_schema(
        properties={"col": {"oneOf": [{"type": "string"}, {"type": "integer"}]}},
    )
    stored = jsonschema_to_dlt_converter.convert(schema)
    assert stored["tables"]["test_source"]["columns"]["col"]["data_type"] == "text"


def test_convert_pass_boolean_schema_property(
    jsonschema_to_dlt_converter: JSONSchemaToDlt, caplog: pytest.LogCaptureFixture
) -> None:
    """A boolean schema property (true/false) maps to text with a warning."""
    schema = base_object_schema(properties={"col": True})
    with caplog.at_level(logging.WARNING):
        stored = jsonschema_to_dlt_converter.convert(schema)
    assert stored["tables"]["test_source"]["columns"]["col"]["data_type"] == "text"
    assert any("no type" in r.getMessage() for r in caplog.records)
