"""Typed dlt emission and structured compatibility with the legacy converter."""

import json
from dataclasses import replace
from decimal import Decimal
from pathlib import Path
from typing import Any, Literal

import pytest
import yaml
from dlt.common.schema import Schema
from pydantic import ValidationError
from pyiceberg.schema import Schema as IcebergSchema
from pyiceberg.types import BinaryType, DecimalType, ListType, NestedField, TimestampType, TimestamptzType

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.ir import Field, NodeHints, SchemaDocument, TypedNode
from cdm_data_loaders.converters.emitters.dlt import DltEmitter
from cdm_data_loaders.converters.readers.dlt import DltReader
from cdm_data_loaders.converters.readers.iceberg import IcebergReader
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader
from tests.cdm_data_loaders.converters.conftest import base_object_schema


@pytest.mark.parametrize(
    "options",
    [
        {},
        {"skip_nested_types": True},
        {"flatten_scalars": False},
        {"max_nesting": 1},
        {"write_disposition": "merge"},
    ],
    ids=["defaults", "skip-nested", "keep-scalar-arrays", "depth-one", "merge-root"],
)
def test_emit_pass_legacy_structure_and_dispatch(options: dict[str, Any]) -> None:
    """Dispatch, outer constraints, array nullability and depth match legacy output."""
    source = base_object_schema(
        description="root",
        required=["flag", "tags"],
        properties={
            "flag": False,
            "tags": {"type": "array", "items": {"type": "integer", "description": "ignored"}},
            "missing": {"type": "array"},
            "enum_first": {"enum": [1, 2], "oneOf": [{"type": "boolean"}]},
            "implicit_first": {"minimum": 1, "oneOf": [{"type": "boolean"}]},
            "outer_number": {"type": "string", "multipleOf": 0.01, "oneOf": [{"type": "number"}]},
            "outer_format": {"format": "date", "oneOf": [{"type": "string", "format": "binary"}]},
            "outer_object": {
                "type": "string",
                "properties": {"outer": {"type": "integer"}},
                "oneOf": [{"type": "object", "properties": {"inner": {"type": "boolean"}}}],
            },
            "deep": {
                "type": "object",
                "properties": {
                    "object": {"type": "object", "properties": {"leaf": {"type": "integer"}}},
                    "array": {
                        "type": "array",
                        "items": {"type": "object", "properties": {"leaf": {"type": "integer"}}},
                    },
                    "scalars": {"type": "array", "items": {"type": "boolean"}},
                },
            },
        },
    )
    document = JsonSchemaReader().read(source)
    columns = {
        "flag": {"name": "flag", "data_type": "text", "nullable": False},
        "missing": {"name": "missing", "data_type": "json", "nullable": True},
        "enum_first": {"name": "enum_first", "data_type": "bigint", "nullable": True},
        "implicit_first": {"name": "implicit_first", "data_type": "double", "nullable": True},
        "outer_number": {"name": "outer_number", "data_type": "decimal", "nullable": True},
        "outer_format": {"name": "outer_format", "data_type": "date", "nullable": True},
    }
    tables = {
        "test_source": {
            "name": "test_source",
            "description": "root",
            "columns": columns,
            "write_disposition": options.get("write_disposition", "append"),
            "resource": "test_source",
        }
    }
    if options.get("skip_nested_types"):
        columns.update(
            {
                "tags": {"name": "tags", "data_type": "json", "nullable": False},
                "outer_object": {"name": "outer_object", "data_type": "json", "nullable": True},
                "deep": {"name": "deep", "data_type": "json", "nullable": True},
            }
        )
    else:
        child_columns = {
            "test_source__outer_object": {"outer": {"name": "outer", "data_type": "bigint", "nullable": True}},
            "test_source__deep": {},
            "test_source__deep__array": {"leaf": {"name": "leaf", "data_type": "bigint", "nullable": True}},
        }
        if options.get("max_nesting") == 1:
            child_columns["test_source__deep"]["object"] = {"name": "object", "data_type": "json", "nullable": True}
        else:
            child_columns["test_source__deep__object"] = {
                "leaf": {"name": "leaf", "data_type": "bigint", "nullable": True}
            }
        if options.get("flatten_scalars") is False:
            columns["tags"] = {"name": "tags", "data_type": "json", "nullable": False}
            child_columns["test_source__deep"]["scalars"] = {"name": "scalars", "data_type": "json", "nullable": True}
        else:
            child_columns["test_source__tags"] = {"value": {"name": "value", "data_type": "bigint", "nullable": True}}
            child_columns["test_source__deep__scalars"] = {
                "value": {"name": "value", "data_type": "bool", "nullable": True}
            }
        tables.update(
            {
                name: {"name": name, "parent": name.rsplit("__", 1)[0], "columns": values}
                for name, values in child_columns.items()
            }
        )
    assert DltEmitter(schema_name="test_source", **options).emit(document) == {
        "name": "test_source",
        "version": 1,
        "previous_hashes": [],
        "engine_version": 11,
        "tables": tables,
    }


@pytest.mark.parametrize("case_id", ["json-dlt-nested", "json-dlt-union-dynamic"], ids=["nested", "union-dynamic"])
def test_emit_pass_captured_goldens(case_id: str) -> None:
    """Captured dlt outputs match exactly as structured stored schemas."""
    cases = json.loads(Path("tests/data/converters/ir/legacy_outputs.json").read_text(encoding="utf-8"))
    case = next(case for case in cases if case["id"] == case_id)
    assert DltEmitter(**case["options"]).emit(JsonSchemaReader().read(case["input"])) == case["expected"]


@pytest.mark.parametrize(
    ("prop", "expected_type", "child"),
    [
        (True, "text", False),
        (False, "text", False),
        ({}, "text", False),
        ({"type": "null"}, "text", False),
        ({"type": "boolean"}, "bool", False),
        ({"type": "integer", "minimum": -5, "maximum": 100}, "bigint", False),
        ({"type": "number"}, "double", False),
        ({"type": "number", "multipleOf": 5}, "decimal", False),
        ({"type": "number", "multipleOf": "invalid"}, "double", False),
        ({"type": "string"}, "text", False),
        ({"type": "string", "format": "date-time"}, "timestamp", False),
        ({"type": "string", "format": "datetime"}, "timestamp", False),
        ({"type": "string", "format": "timestamp"}, "timestamp", False),
        ({"type": "string", "format": "date"}, "date", False),
        ({"type": "string", "format": "time"}, "time", False),
        ({"type": "string", "format": "byte"}, "binary", False),
        ({"type": "string", "format": "binary"}, "binary", False),
        ({"type": "string", "format": "base64"}, "binary", False),
        ({"type": "string", "format": "uuid"}, "text", False),
        ({"type": "unrecognized"}, "text", False),
        ({"type": ["null"]}, "text", False),
        ({"type": []}, "text", False),
        ({"type": ["null", "integer"]}, "bigint", False),
        ({"type": ["string", "integer"], "oneOf": [{"type": "boolean"}]}, "text", False),
        ({"enum": []}, "text", False),
        ({"enum": [True, False]}, "bool", False),
        ({"enum": [1, 2]}, "bigint", False),
        ({"enum": [1, 2.5]}, "double", False),
        ({"enum": [1, "mixed"]}, "text", False),
        ({"enum": [None]}, "text", False),
        ({"pattern": "^x"}, "text", False),
        ({"minimum": 1}, "double", False),
        ({"anyOf": [{"type": "number", "multipleOf": 0.1}]}, "double", False),
        ({"type": "integer", "anyOf": [{"type": "string"}]}, "text", False),
        ({"oneOf": [{"type": "boolean"}], "anyOf": [{"type": "integer"}]}, "bool", False),
        ({"oneOf": [], "anyOf": [{"type": "boolean"}]}, "bool", False),
        ({"anyOf": [False, {"type": "integer"}]}, "text", False),
        ({"anyOf": [True, {"type": "integer"}]}, "text", False),
        ({"type": "object", "oneOf": [{"type": "integer"}]}, "bigint", False),
        ({"type": "array", "oneOf": [{"type": "integer"}]}, "bigint", False),
        ({"type": "array", "items": {}}, "json", False),
        ({"type": "array", "items": False}, "json", False),
        ({"type": "array", "items": []}, "json", False),
        ({"type": "array", "items": [{"type": "integer"}, {"type": "string"}]}, "bigint", True),
        ({"type": "array", "items": {"type": "integer"}, "prefixItems": []}, "json", False),
        ({"type": "array", "prefixItems": [{"type": "boolean"}, {"type": "string"}]}, "bool", True),
        ({"type": "array", "items": {"type": "array", "items": {"type": "integer"}}}, "json", True),
        ({"type": "object", "additionalProperties": {"type": "integer"}}, None, True),
    ],
    ids=[
        "true-schema",
        "false-schema",
        "empty",
        "null",
        "boolean",
        "bounded-integer",
        "number",
        "integral-multiple",
        "invalid-multiple",
        "string",
        "date-time",
        "datetime",
        "timestamp",
        "date",
        "time",
        "byte",
        "binary",
        "base64",
        "uuid",
        "unknown",
        "null-union",
        "empty-union",
        "nullable-integer",
        "multi-union-before-combiner",
        "empty-enum",
        "boolean-enum",
        "integer-enum",
        "numeric-enum",
        "mixed-enum",
        "null-enum",
        "implicit-string",
        "implicit-number",
        "branch-constraints-ignored",
        "declared-before-anyof",
        "oneof-before-anyof",
        "empty-oneof",
        "false-branch",
        "true-branch",
        "object-overridden",
        "array-overridden",
        "empty-items",
        "false-items",
        "empty-tuple",
        "tuple-first",
        "empty-prefix",
        "prefix-first",
        "nested-array",
        "dynamic-object",
    ],
)
@pytest.mark.parametrize("required", [False, True], ids=["optional", "required"])
def test_emit_pass_legacy_property_cases(
    prop: dict[str, Any] | bool, expected_type: str | None, child: bool, required: bool
) -> None:
    """Scalar, tuple, enum and combiner cases preserve legacy property semantics."""
    source = base_object_schema(properties={"col": prop}, required=["col"] if required else [])
    columns = {} if child else {"col": {"name": "col", "data_type": expected_type, "nullable": not required}}
    tables = {"records": {"name": "records", "columns": columns, "write_disposition": "append", "resource": "records"}}
    if child:
        tables["records__col"] = {
            "name": "records__col",
            "parent": "records",
            "columns": {"value": {"name": "value", "data_type": expected_type, "nullable": True}}
            if expected_type
            else {},
        }
    assert DltEmitter(schema_name="records").emit(JsonSchemaReader().read(source)) == {
        "name": "records",
        "version": 1,
        "engine_version": 11,
        "previous_hashes": [],
        "tables": tables,
    }


@pytest.mark.parametrize("mode", ["object", "array"], ids=["object-children", "array-children"])
def test_emit_pass_native_forest_round_trip(mode: Literal["object", "array"]) -> None:
    """Each native root retains scoped hints and explicit non-prefix child names."""
    source = {
        "name": "source",
        "version": 4,
        "previous_hashes": ["previous"],
        "engine_version": 11,
        "settings": {"schema_contract": "evolve"},
        "tables": {
            "root": {
                "name": "root",
                "description": "root rows",
                "write_disposition": "replace",
                "resource": "original_resource",
                "columns": {
                    "id": {"data_type": "bigint", "nullable": False, "precision": 32, "primary_key": True},
                    "amount": {"data_type": "decimal", "precision": 12, "scale": 3},
                    "when": {"data_type": "timestamp", "precision": 3, "timezone": False, "description": None},
                    "wei": {"data_type": "wei", "nullable": None},
                    "binary": {"data_type": "binary", "precision": 16},
                    "incomplete": {},
                    "explicit_null": {"data_type": None, "precision": None, "scale": None, "timezone": None},
                    "json": {"data_type": "json"},
                },
            },
            "rows": {"parent": "root", "description": None, "columns": {"label": {"data_type": "text"}}},
            "root__tags": {
                "parent": "root",
                "description": "tags",
                "columns": {"value": {"data_type": "text", "nullable": False, "description": "tag"}},
            },
            "other": {"columns": {}},
        },
    }
    documents = DltReader(child_table_mode=mode).read(source)
    emitter = DltEmitter(schema_name="source")
    assert emitter.emit(documents["root"]) == {
        **source,
        "tables": {key: value for key, value in source["tables"].items() if key != "other"},
    }
    assert emitter.emit(documents["other"]) == {**source, "tables": {"other": source["tables"]["other"]}}
    overridden = DltEmitter(schema_name="renamed", write_disposition="merge").emit(documents["root"])
    assert overridden["name"] == "renamed"
    assert overridden["tables"]["root"]["write_disposition"] == "merge"
    assert set(overridden["tables"]) == {"root", "rows", "root__tags"}


def test_emit_pass_native_edits_use_nodes_not_raw_tables() -> None:
    """Removed fields and changed types cannot be resurrected from provenance."""
    document = DltReader().read({"root": {"columns": {"old": {"data_type": "text"}, "drop": {"data_type": "bool"}}}})[
        "root"
    ]
    original = document.root.properties[0]
    changed = replace(
        original, name="new", node=replace(original.node, type="integer", hints=NodeHints(), nullable=False)
    )
    document = replace(document, root=replace(document.root, properties=(changed,)))
    assert DltEmitter(schema_name="source").emit(document)["tables"] == {
        "root": {"columns": {"new": {"data_type": "bigint", "nullable": False}}},
    }
    embedded = SchemaDocument(root=TypedNode("object", properties=(Field("nested", document.root),)))
    assert DltEmitter(schema_name="fresh").emit(embedded)["tables"] == {
        "fresh": {"name": "fresh", "columns": {}, "resource": "fresh", "write_disposition": "append"},
        "fresh__nested": {
            "name": "fresh__nested",
            "parent": "fresh",
            "columns": {"new": {"name": "new", "data_type": "bigint", "nullable": False}},
        },
    }


def test_emit_pass_iceberg_logical_types() -> None:
    """Native Iceberg decimals, timestamps, binary and element nullability survive."""
    source = IcebergSchema(
        NestedField(1, "amount", DecimalType(12, 3), required=True, doc="amount"),
        NestedField(2, "when", TimestampType(), required=False),
        NestedField(3, "zoned", TimestamptzType(), required=True),
        NestedField(4, "binary", BinaryType(), required=True),
        NestedField(5, "amounts", ListType(6, DecimalType(8, 2), element_required=True), required=True),
    )
    stored = DltEmitter(schema_name="iceberg").emit(IcebergReader().read(source))
    assert stored["tables"]["iceberg"]["columns"] == {
        "amount": {
            "name": "amount",
            "data_type": "decimal",
            "nullable": False,
            "precision": 12,
            "scale": 3,
            "description": "amount",
        },
        "when": {
            "name": "when",
            "data_type": "timestamp",
            "nullable": True,
            "precision": 6,
            "timezone": False,
            "description": None,
        },
        "zoned": {
            "name": "zoned",
            "data_type": "timestamp",
            "nullable": False,
            "precision": 6,
            "timezone": True,
            "description": None,
        },
        "binary": {"name": "binary", "data_type": "binary", "nullable": False, "description": None},
    }
    assert stored["tables"]["iceberg__amounts"]["columns"] == {
        "value": {"name": "value", "data_type": "decimal", "nullable": False, "precision": 8, "scale": 2},
    }


def test_emit_pass_decimal_constraint_has_no_precision() -> None:
    """Exact JSON multipleOf selects decimal without manufacturing precision."""
    document = JsonSchemaReader().read(
        base_object_schema(properties={"amount": {"type": "number", "multipleOf": Decimal("0.01")}})
    )
    assert DltEmitter(schema_name="records").emit(document)["tables"]["records"]["columns"] == {
        "amount": {"name": "amount", "nullable": True, "data_type": "decimal"},
    }


def test_to_schema_pass_parents_first_and_yaml() -> None:
    """The actual dlt Schema accepts a deep graph and YAML preserves stored output."""
    source = base_object_schema(
        properties={
            "first": {
                "type": "object",
                "properties": {"second": {"type": "object", "properties": {"id": {"type": "integer"}}}},
            }
        }
    )
    document = JsonSchemaReader().read(source)
    emitter = DltEmitter(schema_name="records")
    stored = emitter.emit(document)
    assert list(stored["tables"]) == ["records", "records__first", "records__first__second"]
    live = emitter.to_schema(document)
    assert isinstance(live, Schema)
    assert {name: live.tables[name] for name in stored["tables"]} == stored["tables"]
    assert yaml.safe_load(emitter.to_yaml(document)) == stored
    assert Schema.from_dict(yaml.safe_load(live.to_pretty_yaml())).name == "records"


@pytest.mark.parametrize(
    "options",
    [
        {"schema_name": ""},
        {"schema_name": "  "},
        {"schema_name": "valid", "max_nesting": 0},
        {"schema_name": "valid", "max_nesting": -1},
    ],
    ids=["empty-name", "blank-name", "zero-depth", "negative-depth"],
)
def test_dlt_emitter_fail_invalid_configuration(options: dict[str, Any]) -> None:
    """Invalid names and nonpositive nesting limits fail at construction."""
    with pytest.raises(ValidationError):
        DltEmitter(**options)


def test_dlt_emitter_fail_frozen_configuration() -> None:
    """Configuration cannot be changed after construction."""
    emitter = DltEmitter(schema_name="records")
    with pytest.raises(ValidationError, match="frozen"):
        emitter.max_nesting = 2


def test_emit_fail_unsupported_nodes() -> None:
    """Invalid roots and unknown native column types use the shared error base."""
    emitter = DltEmitter(schema_name="records")
    with pytest.raises(ConversionError, match="root"):
        emitter.emit(SchemaDocument(root=TypedNode("string")))
    with pytest.raises(ConversionError, match="No dlt mapping"):
        emitter.emit(SchemaDocument(root=TypedNode("object", properties=(Field("unknown", TypedNode("unknown")),))))


@pytest.mark.parametrize("include", [False, True], ids=["filtered", "included"])
def test_to_schema_pass_native_incomplete_metadata(include: bool) -> None:
    """Live schema construction supplies names without undoing reader filtering."""
    source = {
        "root": {
            "columns": {
                "id": {"data_type": "bigint", "nullable": False},
                "_dlt_id": {"data_type": "text", "row_key": True},
                "id__v_text": {"data_type": "text", "variant": True},
                "empty": {"nullable": None},
            },
        },
        "root__child": {"parent": "root", "columns": {"value": {"data_type": "text"}}},
    }
    document = DltReader(include_dlt_columns=include, include_variant_columns=include).read(source)["root"]
    emitter = DltEmitter(schema_name="source")
    columns = {name: column for name, column in source["root"]["columns"].items() if include or name in {"id", "empty"}}
    assert emitter.emit(document)["tables"] == {
        "root": {"columns": columns},
        "root__child": source["root__child"],
    }
    live = emitter.to_schema(document)
    assert live.tables["root"]["columns"] == {name: {**column, "name": name} for name, column in columns.items()}
    assert live.tables["root__child"]["columns"] == {"value": {"name": "value", "data_type": "text"}}


def test_emit_pass_native_empty_table_and_new_child() -> None:
    """Missing column mappings stay absent unless actual nodes introduce columns."""
    document = DltReader().read({"root": {"parent": None}})["root"]
    emitter = DltEmitter(schema_name="source")
    assert emitter.emit(document)["tables"] == {"root": {"parent": None}}
    new_child = Field("new", TypedNode("object", properties=(Field("value", TypedNode("string")),)))
    changed = replace(document, root=replace(document.root, properties=(new_child,)))
    assert emitter.emit(changed)["tables"] == {
        "root": {"parent": None},
        "root__new": {
            "name": "root__new",
            "parent": "root",
            "columns": {"value": {"name": "value", "data_type": "text", "nullable": True}},
        },
    }
