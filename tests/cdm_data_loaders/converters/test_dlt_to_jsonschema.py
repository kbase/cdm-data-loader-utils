"""Unit and integration tests for dlt_to_jsonschema."""

import json
from copy import deepcopy
from pathlib import Path
from typing import Any

import pytest
import yaml
from jsonschema import Draft202012Validator

from cdm_data_loaders.converters.dlt_to_jsonschema import DltToJSONSchema, DltToJSONSchemaError
from cdm_data_loaders.converters.emitters.json_schema import DLT_TYPES as TYPE_MAP
from cdm_data_loaders.converters.emitters.json_schema import JSON_SCHEMA_DIALECT
from cdm_data_loaders.converters.jsonschema_to_dlt import JSONSchemaToDlt
from tests.cdm_data_loaders.converters.conftest import (
    base_stored_schema,
    base_table,
)


def test_dlt_to_jsonschema_fail_instance_is_frozen(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Assigning to a DltToJSONSchema instance attribute raises ValidationError (frozen model)."""
    with pytest.raises(Exception, match="Instance is frozen"):
        dlt_to_jsonschema_converter.include_dlt_columns = True


"""convert"""


def test_convert_fail_empty_tables(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """An empty tables dict raises."""
    with pytest.raises(DltToJSONSchemaError, match="contains no tables"):
        dlt_to_jsonschema_converter.convert({"name": "s", "tables": {}})


def test_convert_fail_no_root_tables(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A schema where every table has a parent raises."""
    schema = base_stored_schema(
        {
            "child": {"parent": "parent", "columns": {}},
            "parent": {"parent": "grandparent", "columns": {}},
            "grandparent": {"parent": "child", "columns": {}},
        }
    )
    with pytest.raises(DltToJSONSchemaError, match="No root tables found"):
        dlt_to_jsonschema_converter.convert(schema)


def test_convert_fail_unrecognized_shape(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A dict that is neither a stored schema nor a tables mapping raises."""
    with pytest.raises(DltToJSONSchemaError, match="Input must be a dlt stored schema"):
        dlt_to_jsonschema_converter.convert({"foo": {"bar": 1}})


def test_convert_fail_dangling_parent(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A nested table referencing a missing parent raises."""
    schema = base_stored_schema({"t": base_table(), "t__child": {"parent": "missing", "columns": {}}})
    with pytest.raises(DltToJSONSchemaError, match="references parent 'missing'"):
        dlt_to_jsonschema_converter.convert(schema)


@pytest.mark.parametrize(
    ("tables", "message"),
    [
        ({"root": {"columns": {}}, "loop": {"parent": "loop", "columns": {}}}, "Cycle"),
        (
            {"root": {"columns": {}}, "first": {"parent": "second"}, "second": {"parent": "first"}},
            "Cycle",
        ),
        ({"root": None}, "must be a mapping"),
        ({"root": {"parent": []}}, "parent must be"),
        ({"root": {"columns": None}}, "columns must be a mapping"),
        ({"root": {"columns": {"broken": None}}}, "must be a mapping"),
    ],
    ids=["disconnected-self-cycle", "disconnected-cycle", "null-table", "list-parent", "null-columns", "null-column"],
)
def test_convert_fail_invalid_table_graph(
    dlt_to_jsonschema_converter: DltToJSONSchema, tables: dict[str, Any], message: str
) -> None:
    """Reject disconnected cycles and malformed tables with the converter's error type."""
    with pytest.raises(DltToJSONSchemaError, match=message):
        dlt_to_jsonschema_converter.convert(base_stored_schema(tables))


def test_convert_fail_unknown_data_type(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A column with an unknown data_type raises."""
    schema = base_stored_schema({"t": {"columns": {"col": {"name": "col", "data_type": "hypercube"}}}})
    with pytest.raises(DltToJSONSchemaError, match="unknown dlt data_type 'hypercube'"):
        dlt_to_jsonschema_converter.convert(schema)


def test_convert_pass_draft_2020_12_declared(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Every document declares the draft 2020-12 dialect and a urn $id."""
    schema = base_stored_schema({"t": base_table()})
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["$schema"] == JSON_SCHEMA_DIALECT
    assert docs["t"]["$id"] == "urn:dlt:t"


@pytest.mark.parametrize(
    "schema",
    [
        base_stored_schema({"t": base_table()}),
        {"t": base_table()},
    ],
    ids=["stored-schema", "plain-tables-dict"],
)
def test_convert_accepts_both_input_shapes(
    dlt_to_jsonschema_converter: DltToJSONSchema, schema: dict[str, Any]
) -> None:
    """Both stored-schema dicts and plain tables mappings convert."""
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["properties"]["col"]["type"] == "string"


def test_convert_validates_against_metaschema(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Generated documents are valid draft 2020-12 documents."""
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "a": {"name": "a", "data_type": "bigint", "nullable": False},
                    "b": {"name": "b", "data_type": "text"},
                }
            }
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    validator = Draft202012Validator
    assert validator.check_schema(docs["t"]) is None


"""column type mapping"""


@pytest.mark.parametrize(
    ("data_type", "expected"),
    [
        ("text", {"type": "string"}),
        ("bigint", {"type": "integer"}),
        ("double", {"type": "number"}),
        ("bool", {"type": "boolean"}),
        ("date", {"type": "string", "format": "date"}),
        ("time", {"type": "string", "format": "time"}),
        (
            "timestamp",
            {"type": "string", "format": "date-time", "x-dlt": {"data_type": "timestamp"}},
        ),
        (
            "decimal",
            {"type": "string", "pattern": r"^-?\d+(\.\d+)?$", "x-dlt": {"data_type": "decimal"}},
        ),
        ("binary", {"type": "string", "contentEncoding": "base64"}),
        ("json", {}),
        ("wei", {"type": "integer", "x-dlt": {"data_type": "wei"}}),
    ],
    ids=[
        "text",
        "bigint",
        "double",
        "bool",
        "date",
        "time",
        "timestamp",
        "decimal",
        "binary",
        "json",
        "wei",
    ],
)
def test_convert_column_types_pass(
    dlt_to_jsonschema_converter: DltToJSONSchema,
    data_type: str,
    expected: dict[str, Any],
) -> None:
    """Each dlt data type maps to the expected draft 2020-12 property schema (non-nullable)."""
    schema = base_stored_schema({"t": {"columns": {"col": {"name": "col", "data_type": data_type, "nullable": False}}}})
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["properties"]["col"] == expected


@pytest.mark.parametrize("data_type", ["decimal", "timestamp", "wei"], ids=["decimal", "timestamp", "wei"])
def test_convert_pass_type_map_metadata_isolation(dlt_to_jsonschema_converter: DltToJSONSchema, data_type: str) -> None:
    """Keep templates, sibling columns, repeated calls and returned metadata independent."""
    template = deepcopy(TYPE_MAP[data_type])
    hints = {"precision": 18, "scale": 4, "primary_key": True}
    plain_column = {"name": "plain", "data_type": data_type, "nullable": False}
    plain_schema = {"table": {"columns": {"plain": plain_column}}}
    schema = {
        "table": {
            "columns": {
                "hinted": {**plain_column, "name": "hinted", **hints},
                "plain": plain_column,
            }
        }
    }
    original_schema = deepcopy(schema)
    try:
        properties = dlt_to_jsonschema_converter.convert(schema)["table"]["properties"]
        expected_hinted = {**template, "x-dlt": {"data_type": data_type, **hints}}
        assert properties == {"hinted": expected_hinted, "plain": template}
        assert TYPE_MAP[data_type] == template
        assert schema == original_schema
        properties["hinted"]["x-dlt"]["scale"] = 2
        assert properties["plain"] == template
        assert dlt_to_jsonschema_converter.convert(plain_schema)["table"]["properties"] == {"plain": template}
        assert DltToJSONSchema().convert(plain_schema)["table"]["properties"] == {"plain": template}
        assert TYPE_MAP[data_type] == template
    finally:
        TYPE_MAP[data_type].clear()
        TYPE_MAP[data_type].update(template)


def test_convert_incomplete_column_accepts_anything(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A column with no data_type (incomplete) converts to an unconstrained schema."""
    schema = base_stored_schema({"t": {"columns": {"col": {"name": "col"}}}})
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["properties"]["col"] == {}


def test_convert_nullable_uses_anyof_null(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Nullable columns wrap the type in anyOf with null (repo convention)."""
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "req": {"name": "req", "data_type": "text", "nullable": False},
                    "opt": {"name": "opt", "data_type": "text", "nullable": True},
                }
            }
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    props = docs["t"]["properties"]
    assert props["req"] == {"type": "string"}
    assert props["opt"] == {"anyOf": [{"type": "string"}, {"type": "null"}]}
    assert docs["t"]["required"] == ["req"]


def test_convert_description_propagates(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Column and table descriptions propagate to the JSON Schema document."""
    schema = base_stored_schema(
        {
            "t": {
                "description": "table description",
                "columns": {
                    "col": {
                        "name": "col",
                        "data_type": "text",
                        "nullable": False,
                        "description": "column description",
                    }
                },
            }
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["description"] == "table description"
    assert docs["t"]["title"] == "t"
    assert docs["t"]["properties"]["col"]["description"] == "column description"


def test_convert_precision_scale_in_x_dlt(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Precision and scale hints land in the x-dlt metadata block."""
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "col": {"name": "col", "data_type": "decimal", "nullable": False, "precision": 38, "scale": 9}
                }
            }
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["properties"]["col"]["x-dlt"] == {"data_type": "decimal", "precision": 38, "scale": 9}


def test_convert_unknown_hints_preserved(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Non-standard column hints survive in the x-dlt block."""
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "col": {
                        "name": "col",
                        "data_type": "text",
                        "nullable": False,
                        "primary_key": True,
                        "unique": True,
                    },
                }
            }
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["properties"]["col"]["x-dlt"] == {"primary_key": True, "unique": True}


def test_convert_unknown_hints_dropped_when_disabled() -> None:
    """With preserve_unknown_hints=False, non-standard hints are dropped."""
    conv = DltToJSONSchema(preserve_unknown_hints=False)
    schema = base_stored_schema({"t": {"columns": {"col": {"name": "col", "data_type": "text", "primary_key": True}}}})
    docs = conv.convert(schema)
    assert "x-dlt" not in docs["t"]["properties"]["col"]


"""dlt internal columns"""


def test_convert_drops_dlt_internal_columns(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Dlt's `_dlt_*` bookkeeping columns are dropped by default."""
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "col": {"name": "col", "data_type": "text"},
                    "_dlt_id": {"name": "_dlt_id", "data_type": "text", "nullable": False},
                    "_dlt_load_id": {"name": "_dlt_load_id", "data_type": "text", "nullable": False},
                }
            }
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert set(docs["t"]["properties"].keys()) == {"col"}
    # the remaining col is nullable, so no required list is emitted
    assert "required" not in docs["t"]


def test_convert_keeps_dlt_columns_when_enabled() -> None:
    """With include_dlt_columns=True, _dlt_* columns are kept."""
    conv = DltToJSONSchema(include_dlt_columns=True)
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "col": {"name": "col", "data_type": "text"},
                    "_dlt_id": {"name": "_dlt_id", "data_type": "text", "nullable": False},
                }
            }
        }
    )
    docs = conv.convert(schema)
    assert "_dlt_id" in docs["t"]["properties"]
    assert docs["t"]["required"] == ["_dlt_id"]


def test_convert_drops_variant_columns(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Variant columns (col__v_text) are dropped by default."""
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "col": {"name": "col", "data_type": "text"},
                    "col__v_text": {"name": "col__v_text", "data_type": "text"},
                }
            }
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert set(docs["t"]["properties"].keys()) == {"col"}


def test_convert_keeps_variant_columns_when_enabled() -> None:
    """With include_variant_columns=True, variant columns are kept."""
    conv = DltToJSONSchema(include_variant_columns=True)
    schema = base_stored_schema(
        {
            "t": {
                "columns": {
                    "col": {"name": "col", "data_type": "text"},
                    "col__v_text": {"name": "col__v_text", "data_type": "text"},
                }
            }
        }
    )
    docs = conv.convert(schema)
    assert "col__v_text" in docs["t"]["properties"]


"""nested child tables"""


def test_convert_child_table_becomes_object_property(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A nested child table folds back into its parent as an object property."""
    schema = base_stored_schema(
        {
            "t": {"columns": {"col": {"name": "col", "data_type": "text", "nullable": False}}},
            "t__address": {
                "parent": "t",
                "columns": {"city": {"name": "city", "data_type": "text", "nullable": False}},
            },
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    props = docs["t"]["properties"]
    assert props["address"]["type"] == "object"
    assert props["address"]["properties"]["city"]["type"] == "string"
    assert set(docs.keys()) == {"t"}


def test_convert_deep_child_chain_folds_recursively(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """Multi-level child chains fold recursively into the root document."""
    schema = base_stored_schema(
        {
            "t": {"columns": {}},
            "t__a": {"parent": "t", "columns": {"x": {"name": "x", "data_type": "text", "nullable": False}}},
            "t__a__b": {
                "parent": "t__a",
                "columns": {"y": {"name": "y", "data_type": "bigint", "nullable": False}},
            },
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    props = docs["t"]["properties"]
    assert props["a"]["properties"]["b"]["properties"]["y"]["type"] == "integer"


def test_convert_scalar_value_table_becomes_array(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A child table with a single `value` column becomes an array of that value type."""
    schema = base_stored_schema(
        {
            "t": {"columns": {}},
            "t__tags": {
                "parent": "t",
                "columns": {"value": {"name": "value", "data_type": "bigint", "nullable": False}},
            },
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    props = docs["t"]["properties"]
    assert props["tags"]["type"] == "array"
    assert props["tags"]["items"] == {"type": "integer"}


def test_convert_scalar_value_table_nullable_items(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A nullable scalar-value column yields anyOf null items."""
    schema = base_stored_schema(
        {
            "t": {"columns": {}},
            "t__tags": {"parent": "t", "columns": {"value": {"name": "value", "data_type": "text"}}},
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert docs["t"]["properties"]["tags"]["items"] == {"anyOf": [{"type": "string"}, {"type": "null"}]}


def test_convert_child_table_mode_array() -> None:
    """With child_table_mode='array', object children render as arrays of objects."""
    conv = DltToJSONSchema(child_table_mode="array")
    schema = base_stored_schema(
        {
            "t": {"columns": {}},
            "t__orders": {
                "parent": "t",
                "columns": {"total": {"name": "total", "data_type": "double", "nullable": False}},
            },
        }
    )
    docs = conv.convert(schema)
    props = docs["t"]["properties"]
    assert props["orders"]["type"] == "array"
    assert props["orders"]["items"]["properties"]["total"]["type"] == "number"


def test_convert_child_required_propagates_to_parent(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A child with a non-nullable column appears in the parent's required list."""
    schema = base_stored_schema(
        {
            "t": {"columns": {}},
            "t__orders": {
                "parent": "t",
                "columns": {"total": {"name": "total", "data_type": "double", "nullable": False}},
            },
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    assert "orders" in docs["t"]["required"]


def test_convert_child_table_dlt_columns_skipped_in_required(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A child with only _dlt_* columns contributes nothing to the parent's required list."""
    schema = base_stored_schema(
        {
            "t": {"columns": {}},
            "t__child": {
                "parent": "t",
                "columns": {
                    "_dlt_id": {"name": "_dlt_id", "data_type": "text", "nullable": False},
                },
            },
        }
    )
    docs = dlt_to_jsonschema_converter.convert(schema)
    # _dlt_* columns are skipped, so nothing contributes to required
    assert "required" not in docs["t"]


"""entry points"""


def test_convert_from_string_pass_json(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """convert_from_string parses JSON text."""
    schema = base_stored_schema({"t": base_table()})
    docs = dlt_to_jsonschema_converter.convert_from_string(json.dumps(schema))
    assert docs["t"]["properties"]["col"]["type"] == "string"


def test_convert_from_string_pass_yaml(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """convert_from_string falls back to YAML when JSON parsing fails."""
    schema = base_stored_schema({"t": base_table()})
    docs = dlt_to_jsonschema_converter.convert_from_string(yaml.dump(schema))
    assert docs["t"]["properties"]["col"]["type"] == "string"


def test_convert_from_file_pass_json(dlt_to_jsonschema_converter: DltToJSONSchema, tmp_path: Path) -> None:
    """convert_from_file reads a .json schema file."""
    path = tmp_path / "schema.json"
    path.write_text(json.dumps(base_stored_schema({"t": base_table()})))
    docs = dlt_to_jsonschema_converter.convert_from_file(str(path))
    assert docs["t"]["$schema"] == JSON_SCHEMA_DIALECT


def test_convert_from_file_pass_yaml(dlt_to_jsonschema_converter: DltToJSONSchema, tmp_path: Path) -> None:
    """convert_from_file reads a .yaml schema file."""
    path = tmp_path / "schema.yaml"
    path.write_text(yaml.dump(base_stored_schema({"t": base_table()})))
    docs = dlt_to_jsonschema_converter.convert_from_file(str(path))
    assert docs["t"]["$schema"] == JSON_SCHEMA_DIALECT


def test_convert_from_file_fail_missing_file(dlt_to_jsonschema_converter: DltToJSONSchema, tmp_path: Path) -> None:
    """convert_from_file raises FileNotFoundError for a nonexistent path."""
    with pytest.raises(FileNotFoundError):
        dlt_to_jsonschema_converter.convert_from_file(str(tmp_path / "nope.json"))


"""round trip with jsonschema_to_dlt"""


def test_round_trip_with_forward_converter(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A JSON Schema -> dlt -> JSON Schema round trip preserves types and structure."""
    original = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "required": ["id", "name"],
        "properties": {
            "id": {"type": "integer"},
            "name": {"type": "string"},
            "joined": {"type": "string", "format": "date-time"},
            "active": {"type": "boolean"},
            "tags": {"type": "array", "items": {"type": "string"}},
        },
    }
    stored = JSONSchemaToDlt(schema_name="rt").convert(original)
    docs = dlt_to_jsonschema_converter.convert(stored)
    doc = docs["rt"]

    props = doc["properties"]
    assert props["id"] == {"type": "integer"}
    assert props["name"] == {"type": "string"}
    assert props["joined"]["anyOf"][0]["format"] == "date-time"
    assert props["active"]["anyOf"][0]["type"] == "boolean"
    assert props["tags"]["type"] == "array"
    assert props["tags"]["items"]["anyOf"][0] == {"type": "string"}
    assert sorted(doc["required"]) == ["id", "name"]


def test_round_trip_nested_object(dlt_to_jsonschema_converter: DltToJSONSchema) -> None:
    """A nested object survives the round trip as a folded property."""
    original = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "address": {
                "type": "object",
                "properties": {"city": {"type": "string"}},
                "required": ["city"],
            },
        },
    }
    stored = JSONSchemaToDlt(schema_name="rt").convert(original)
    docs = dlt_to_jsonschema_converter.convert(stored)
    props = docs["rt"]["properties"]
    assert props["address"]["type"] == "object"
    assert props["address"]["properties"]["city"]["type"] == "string"
