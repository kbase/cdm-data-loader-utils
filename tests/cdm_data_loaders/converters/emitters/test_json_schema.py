"""JSON Schema emission compatibility and structural preservation tests."""

import json
from copy import deepcopy
from dataclasses import replace
from decimal import Decimal
from pathlib import Path
from typing import Any, Final, Literal, cast

import pytest
from pyiceberg.schema import Schema
from pyiceberg.types import (
    BinaryType,
    BooleanType,
    DateType,
    DecimalType,
    DoubleType,
    FixedType,
    FloatType,
    IcebergType,
    IntegerType,
    ListType,
    LongType,
    MapType,
    NestedField,
    StringType,
    StructType,
    TimestampNanoType,
    TimestampType,
    TimestamptzNanoType,
    TimestamptzType,
    TimeType,
    UnknownType,
    UUIDType,
)

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionError, Extensions, ExtensionSpec
from cdm_data_loaders.converters.core.ir import Field, NodeHints, Provenance, SchemaDocument, TypedNode
from cdm_data_loaders.converters.dlt_to_jsonschema import DltToJSONSchema
from cdm_data_loaders.converters.emitters.json_schema import JSON_SCHEMA_DIALECT, JsonSchemaEmitter, decimal_pattern
from cdm_data_loaders.converters.readers.dlt import DltReader
from cdm_data_loaders.converters.readers.iceberg import IcebergReader
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader
from tests.cdm_data_loaders.converters.conftest import make_table

GOLDEN_PATH: Final = Path(__file__).parents[3] / "data" / "converters" / "ir" / "legacy_outputs.json"
CASES: Final[list[dict[str, Any]]] = [
    case
    for case in json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))
    if case["route"] in {"dlt-json", "iceberg-json"}
]

DLT_EXPECTED: Final[dict[str | None, dict[str, Any]]] = {
    None: {},
    "json": {},
    "text": {"type": "string"},
    "bigint": {"type": "integer"},
    "double": {"type": "number"},
    "bool": {"type": "boolean"},
    "decimal": {"type": "string", "pattern": r"^-?\d+(\.\d+)?$", "x-dlt": {"data_type": "decimal"}},
    "timestamp": {"type": "string", "format": "date-time", "x-dlt": {"data_type": "timestamp"}},
    "date": {"type": "string", "format": "date"},
    "time": {"type": "string", "format": "time"},
    "binary": {"type": "string", "contentEncoding": "base64"},
    "wei": {"type": "integer", "x-dlt": {"data_type": "wei"}},
}


@pytest.mark.parametrize(
    "data_type",
    [None, "text", "bigint", "double", "decimal", "bool", "timestamp", "date", "time", "binary", "json", "wei"],
    ids=[
        "incomplete",
        "text",
        "bigint",
        "double",
        "decimal",
        "bool",
        "timestamp",
        "date",
        "time",
        "binary",
        "json",
        "wei",
    ],
)
@pytest.mark.parametrize("nullable", [True, False, None], ids=["nullable", "required", "unresolved"])
@pytest.mark.parametrize("preserve", [True, False], ids=["keep-hints", "drop-hints"])
def test_emit_pass_dlt_scalar_compatibility(data_type: str | None, nullable: bool | None, preserve: bool) -> None:
    """Logical templates, nullability and hint filtering match the source converter."""
    column = {"name": "value", "data_type": data_type, "nullable": nullable}
    source: dict[str, Any] = {
        "records": {
            "columns": {
                "plain": column,
                "hinted": {**column, "precision": None, "scale": 2, "description": None, "primary_key": True},
            }
        }
    }
    document = DltReader().read(source)["records"]
    plain = deepcopy(DLT_EXPECTED[data_type])
    hinted = deepcopy(plain)
    hinted.setdefault("x-dlt", {}).update({"precision": None, "scale": 2})
    hinted["description"] = None
    if preserve:
        hinted["x-dlt"]["primary_key"] = True
    expected = {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:dlt:records",
        "title": "records",
        "type": "object",
        "additionalProperties": False,
        "properties": {
            "plain": {"anyOf": [plain, {"type": "null"}]} if nullable and plain else plain,
            "hinted": {"anyOf": [hinted, {"type": "null"}]} if nullable else hinted,
        },
    }
    if nullable is False:
        expected["required"] = ["plain", "hinted"]
    assert JsonSchemaEmitter(preserve_unknown_hints=preserve).emit(document) == expected


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_emit_pass_fixed_corpus(case: dict[str, Any]) -> None:
    """Reader and emitter composition reproduces the mechanically captured dictionaries."""
    source = case["input"]
    if case["route"] == "dlt-json":
        options = {key: value for key, value in case["options"].items() if key != "preserve_unknown_hints"}
        emitter = JsonSchemaEmitter(preserve_unknown_hints=case["options"].get("preserve_unknown_hints", True))
        actual = {name: emitter.emit(document) for name, document in DltReader(**options).read(source).items()}
    else:
        identifier = tuple(source["identifier"])
        table = make_table(identifier, Schema.model_validate(source["schema"]), properties=source["properties"])
        actual = JsonSchemaEmitter().emit(IcebergReader().read_table(table, identifier))
    assert actual == case["expected"]


@pytest.mark.parametrize(
    ("field_type", "scalar"),
    [
        (BooleanType(), {"type": "boolean"}),
        (IntegerType(), {"type": "integer", "minimum": -2147483648, "maximum": 2147483647}),
        (LongType(), {"type": "integer", "minimum": -9223372036854775808, "maximum": 9223372036854775807}),
        (FloatType(), {"type": "number"}),
        (DoubleType(), {"type": "number"}),
        (StringType(), {"type": "string"}),
        (DateType(), {"type": "string", "format": "date"}),
        (TimeType(), {"type": "string", "format": "time"}),
        (UUIDType(), {"type": "string", "format": "uuid"}),
        (
            TimestampType(),
            {
                "type": "string",
                "description": "ISO 8601 timestamp with no timezone attached.",
                "x-iceberg": {"logical_type": "timestamp-without-tz"},
            },
        ),
        (TimestamptzType(), {"type": "string", "format": "date-time"}),
        (
            DecimalType(10, 2),
            {
                "type": "string",
                "pattern": r"^-?\d{1,8}(\.\d{1,2})?$",
                "x-iceberg": {"logical_type": "decimal", "precision": 10, "scale": 2},
            },
        ),
        (
            DecimalType(2, 2),
            {
                "type": "string",
                "pattern": r"^-?0(\.\d{1,2})?$",
                "x-iceberg": {"logical_type": "decimal", "precision": 2, "scale": 2},
            },
        ),
        (
            DecimalType(10, 0),
            {
                "type": "string",
                "pattern": r"^-?\d{1,10}$",
                "x-iceberg": {"logical_type": "decimal", "precision": 10, "scale": 0},
            },
        ),
        (
            FixedType(16),
            {"type": "string", "contentEncoding": "base64", "x-iceberg": {"logical_type": "fixed", "length": 16}},
        ),
    ],
    ids=[
        "boolean",
        "integer",
        "long",
        "float",
        "double",
        "string",
        "date",
        "time",
        "uuid",
        "naive",
        "aware",
        "decimal",
        "fraction-only",
        "integer-decimal",
        "fixed",
    ],
)
@pytest.mark.parametrize("required", [True, False], ids=["required", "optional"])
@pytest.mark.parametrize("description", [None, "", "field documentation"], ids=["no-doc", "empty-doc", "doc"])
def test_emit_pass_iceberg_scalar_compatibility(
    field_type: IcebergType, scalar: dict[str, Any], required: bool, description: str | None
) -> None:
    """Scalar types and field descriptions retain exact metadata and null-wrapper placement."""
    schema = Schema(NestedField(1, "value", field_type, required=required, doc=description))
    value = deepcopy(scalar) if required else {"anyOf": [deepcopy(scalar), {"type": "null"}]}
    value.setdefault("x-iceberg", {}).update({"field_id": 1, "required": required})
    if description:
        value["description"] = description
    expected = {
        "type": "object",
        "properties": {"value": value},
        "additionalProperties": False,
        "$schema": JSON_SCHEMA_DIALECT,
        "x-iceberg": {"schema_id": schema.schema_id, "identifier_field_ids": []},
    }
    if required:
        expected["required"] = ["value"]
    assert JsonSchemaEmitter().emit(IcebergReader().read(schema)) == expected


@pytest.mark.parametrize("required", [True, False], ids=["required", "optional"])
@pytest.mark.parametrize("elements_required", [True, False], ids=["required-elements", "optional-elements"])
@pytest.mark.parametrize("description", [None, "", "nested data"], ids=["no-doc", "empty-doc", "doc"])
def test_emit_pass_iceberg_nested_compatibility(
    required: bool,
    elements_required: bool,
    description: str | None,
) -> None:
    """Nested map values, list elements and struct fields preserve independent wrappers."""
    schema = Schema(
        NestedField(
            1,
            "list",
            ListType(2, TimestampType(), element_required=elements_required),
            required=required,
            doc=description,
        ),
        NestedField(
            3,
            "map",
            MapType(
                4,
                IntegerType(),
                5,
                StructType(NestedField(6, "amount", DecimalType(9, 3), required=required, doc=description)),
                value_required=elements_required,
            ),
            required=required,
            doc=description,
        ),
    )
    identifier = ("namespace", "records")
    table = make_table(identifier, schema, properties={"comment": description or ""})
    timestamp = {
        "type": "string",
        "description": "ISO 8601 timestamp with no timezone attached.",
        "x-iceberg": {"logical_type": "timestamp-without-tz"},
    }
    amount = {
        "type": "string",
        "pattern": r"^-?\d{1,6}(\.\d{1,3})?$",
        "x-iceberg": {"logical_type": "decimal", "precision": 9, "scale": 3},
    }
    amount = amount if required else {"anyOf": [amount, {"type": "null"}]}
    amount.setdefault("x-iceberg", {}).update({"field_id": 6, "required": required})
    if description:
        amount["description"] = description
    value = {"type": "object", "properties": {"amount": amount}, "additionalProperties": False}
    if required:
        value["required"] = ["amount"]
    properties = {
        "list": {
            "type": "array",
            "items": timestamp if elements_required else {"anyOf": [timestamp, {"type": "null"}]},
        },
        "map": {
            "type": "array",
            "description": "Iceberg map rendered as key/value pairs (JSON object keys must be strings).",
            "items": {
                "type": "object",
                "properties": {
                    "key": {"type": "integer", "minimum": -2147483648, "maximum": 2147483647},
                    "value": value if elements_required else {"anyOf": [value, {"type": "null"}]},
                },
                "required": ["key", "value"],
                "additionalProperties": False,
            },
        },
    }
    for name, field_id in (("list", 1), ("map", 2)):
        if not required:
            properties[name] = {"anyOf": [properties[name], {"type": "null"}]}
        properties[name]["x-iceberg"] = {"field_id": field_id, "required": required}
        if description:
            properties[name]["description"] = description
    expected = {
        "type": "object",
        "properties": properties,
        "additionalProperties": False,
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:iceberg:namespace.records",
        "title": "records",
        "x-iceberg": {
            "identifier": ["namespace", "records"],
            "schema_id": table.schema().schema_id,
            "format_version": table.metadata.format_version,
            "current_snapshot_id": None,
            "location": "file:///tmp/namespace/records",
            "properties": {"comment": description or ""},
            "partition_spec": [],
            "identifier_field_ids": [],
        },
    }
    if required:
        expected["required"] = ["list", "map"]
    if description:
        expected["description"] = description
    assert JsonSchemaEmitter().emit(IcebergReader().read_table(table, identifier)) == expected


def test_emit_into_pass_adds_table_under_defs() -> None:
    """A table document is added to a namespace document's $defs, keyed by its name."""
    schema = Schema(NestedField(1, "id", LongType(), required=True))
    identifier = ("namespace", "records")
    table = make_table(identifier, schema)
    document = IcebergReader().read_table(table, identifier)
    namespace_doc: dict[str, Any] = {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:iceberg:namespace",
        "title": "namespace",
    }

    result = JsonSchemaEmitter().emit_into(document, namespace_doc)

    expected_entry = JsonSchemaEmitter().emit(document)
    del expected_entry["$schema"]

    assert result is namespace_doc
    assert set(namespace_doc["$defs"]) == {"records"}
    assert "$schema" not in namespace_doc["$defs"]["records"]
    assert namespace_doc["$defs"]["records"] == expected_entry


def test_emit_into_pass_preserves_existing_entries() -> None:
    """Adding a second table keeps the first table's entry under $defs."""
    schema = Schema(NestedField(1, "id", LongType(), required=True))
    first_table = make_table(("namespace", "first"), schema)
    second_table = make_table(("namespace", "second"), schema)
    namespace_doc: dict[str, Any] = {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:iceberg:namespace",
        "title": "namespace",
    }
    emitter = JsonSchemaEmitter()

    emitter.emit_into(IcebergReader().read_table(first_table, ("namespace", "first")), namespace_doc)
    emitter.emit_into(IcebergReader().read_table(second_table, ("namespace", "second")), namespace_doc)

    assert set(namespace_doc["$defs"]) == {"first", "second"}


def test_emit_into_fail_unnamed_document() -> None:
    """A document without a name cannot be added to a namespace document."""
    schema = Schema(NestedField(1, "id", LongType(), required=True))
    document = IcebergReader().read(schema)
    with pytest.raises(ConversionError, match="named document"):
        JsonSchemaEmitter().emit_into(document, {})


@pytest.mark.parametrize(
    "schema",
    [
        {},
        {"type": ["object"], "properties": {}, "required": [], "patternProperties": {}},
        {
            "type": ["null", "object"],
            "required": ["missing", "missing"],
            "properties": {"anything": True, "nothing": False},
        },
        {"properties": {"values": {"items": [True, False, {"type": "null"}], "additionalItems": False}}},
        {
            "properties": {
                "values": {
                    "prefixItems": [{"type": "boolean"}],
                    "items": True,
                    "contains": {"const": True},
                    "unevaluatedItems": False,
                }
            }
        },
        {
            "additionalProperties": {"type": "integer"},
            "patternProperties": {"^a": {"type": "string"}},
            "propertyNames": {"pattern": "^[a-z]+$"},
            "unevaluatedProperties": False,
        },
        {
            "anyOf": [{"properties": {}}, False],
            "oneOf": [{"required": ["a"]}, {"required": ["b"]}],
            "not": {"required": ["c"]},
            "if": True,
            "then": {},
            "else": False,
        },
        {
            "$defs": {"data": {"type": "string"}},
            "definitions": {},
            "dependentSchemas": {"a": False},
            "dependencies": {"a": ["b"], "b": {"required": ["c"]}},
            "contentSchema": {"type": "null"},
        },
        {
            "enum": [{"$ref": "literal", "allOf": [1]}],
            "const": {"$dynamicRef": "literal"},
            "default": {"$recursiveRef": "literal"},
            "examples": [{"properties": {"$ref": "data"}}],
        },
        {
            "properties": {
                "value": {
                    "type": "number",
                    "multipleOf": Decimal("0.00000000000000000001"),
                    "minimum": Decimal("-100.01"),
                    "exclusiveMaximum": Decimal("100.01"),
                }
            }
        },
    ],
    ids=[
        "absent-structure",
        "explicit-empty",
        "boolean-schemas",
        "tuple-items",
        "prefix-items",
        "dynamic-object",
        "combiners",
        "schema-maps",
        "literal-references",
        "decimal-constraints",
    ],
)
@pytest.mark.parametrize(
    "dialect", [JSON_SCHEMA_DIALECT, "http://json-schema.org/draft-07/schema#"], ids=["2020-12", "draft7"]
)
def test_emit_pass_json_structural_roundtrip(schema: dict[str, Any], dialect: str) -> None:
    """Retained schema structure, dialect and literal values survive without reinterpretation."""
    source = {"$schema": dialect, "$id": "urn:json:records", **schema}
    original = deepcopy(source)
    document = JsonSchemaReader().read(source)
    emitter = JsonSchemaEmitter()
    actual = emitter.emit(document)
    assert actual == original
    assert list(actual) == list(original)
    actual["examples"] = [{"changed": True}]
    assert emitter.emit(document) == original
    assert source == original


@pytest.mark.parametrize(
    ("node", "expected"),
    [
        (TypedNode("any"), {}),
        (TypedNode("never"), {"not": {}}),
        (TypedNode("null"), {"type": "null"}),
        (TypedNode("boolean", nullable=True), {"anyOf": [{"type": "boolean"}, {"type": "null"}]}),
        (
            TypedNode("object", properties=(Field("flag", TypedNode("boolean"), required=True),)),
            {"type": "object", "properties": {"flag": {"type": "boolean"}}, "required": ["flag"]},
        ),
        (TypedNode("array", items=TypedNode("any")), {"type": "array", "items": {}}),
        (TypedNode("any", annotations={"description": "unrestricted"}), {"description": "unrestricted"}),
        (
            TypedNode("map", key_type=TypedNode("integer"), value_type=TypedNode("any")),
            {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {"key": {"type": "integer"}, "value": {}},
                    "required": ["key", "value"],
                    "additionalProperties": False,
                },
            },
        ),
    ],
    ids=["any", "never", "null", "nullable-boolean", "object", "array", "annotated-any", "map"],
)
def test_emit_pass_source_neutral_nodes(node: TypedNode, expected: dict[str, Any]) -> None:
    """Direct typed nodes emit structural schemas without a fabricated source dialect."""
    assert JsonSchemaEmitter().emit(SchemaDocument(node)) == expected


@pytest.mark.parametrize("mode", ["object", "array"], ids=["objects", "arrays"])
@pytest.mark.parametrize("include_internal", [True, False], ids=["include-internal", "exclude-internal"])
@pytest.mark.parametrize("include_variants", [True, False], ids=["include-variants", "exclude-variants"])
def test_emit_pass_dlt_nested_compatibility(
    mode: Literal["object", "array"],
    include_internal: bool,
    include_variants: bool,
) -> None:
    """Child policies, descriptions, filtering and duplicate required order match legacy output."""
    source: dict[str, Any] = {
        "root": {
            "description": "root",
            "columns": {
                "child": {"data_type": "text", "nullable": False},
                "_dlt_id": {"data_type": "text", "nullable": False},
                "extra__v_text": {"data_type": "text", "nullable": False},
            },
        },
        "root__child": {
            "parent": "root",
            "description": None,
            "columns": {
                "name": {"data_type": "text", "nullable": False, "description": "name"},
            },
        },
        "root__child__values": {
            "parent": "root__child",
            "description": "values",
            "columns": {
                "value": {"data_type": "json"},
                "_dlt_id": {"data_type": "text", "nullable": False},
            },
        },
        "root__descendant": {"parent": "root", "columns": {}},
        "root__descendant__required": {
            "parent": "root__descendant",
            "columns": {
                "value": {"data_type": "bool", "nullable": False},
            },
        },
        "empty": {"description": "", "columns": {}},
    }
    original = deepcopy(source)
    documents = DltReader(
        child_table_mode=mode,
        include_dlt_columns=include_internal,
        include_variant_columns=include_variants,
    ).read(source)
    values = {"type": "array", "items": {}, "description": "values"}
    if include_internal:
        values_object = {
            "type": "object",
            "properties": {"value": {}, "_dlt_id": {"type": "string"}},
            "required": ["_dlt_id"],
            "additionalProperties": False,
        }
        values = {"type": "array", "items": values_object} if mode == "array" else values_object
        values["description"] = "values"
    child = {
        "type": "object",
        "properties": {"name": {"type": "string", "description": "name"}, "values": values},
        "required": ["name", "values"] if include_internal else ["name"],
        "additionalProperties": False,
    }
    child = {"type": "array", "items": child} if mode == "array" else child
    child["description"] = None
    descendant = {
        "type": "object",
        "properties": {"required": {"type": "array", "items": {"type": "boolean"}}},
        "required": ["required"],
        "additionalProperties": False,
    }
    properties = {"child": child}
    required = ["child"]
    for name, include in (("_dlt_id", include_internal), ("extra__v_text", include_variants)):
        if include:
            properties[name] = {"type": "string"}
            required.append(name)
    properties["descendant"] = {"type": "array", "items": descendant} if mode == "array" else descendant
    required.append("child")
    expected = {
        "root": {
            "type": "object",
            "properties": properties,
            "required": required,
            "additionalProperties": False,
            "description": "root",
            "$schema": JSON_SCHEMA_DIALECT,
            "$id": "urn:dlt:root",
            "title": "root",
        },
        "empty": {
            "type": "object",
            "properties": {},
            "additionalProperties": False,
            "description": "",
            "$schema": JSON_SCHEMA_DIALECT,
            "$id": "urn:dlt:empty",
            "title": "empty",
        },
    }
    actual = {name: JsonSchemaEmitter().emit(document) for name, document in documents.items()}
    assert actual == expected
    assert actual["root"]["required"].count("child") == len(("column", "child-table"))
    assert source == original


@pytest.mark.parametrize(
    "hints",
    [
        {},
        {"precision": None, "scale": None},
        {"precision": 0},
        {"scale": 0},
        {"precision": 38, "scale": 9},
        {"timezone": False, "primary_key": None, "unique": False},
        {"description": ""},
        {"description": None},
    ],
    ids=[
        "absent",
        "null-precision",
        "zero-precision",
        "zero-scale",
        "precision-scale",
        "boolean-hints",
        "empty-doc",
        "null-doc",
    ],
)
@pytest.mark.parametrize(
    "data_type", [None, "json", "decimal", "timestamp", "wei"], ids=["incomplete", "any", "decimal", "timestamp", "wei"]
)
@pytest.mark.parametrize("preserve", [True, False], ids=["keep-hints", "drop-hints"])
def test_emit_pass_dlt_hint_combinations(hints: dict[str, Any], data_type: str | None, preserve: bool) -> None:
    """Absent and null hints remain distinct without leaking source column keys."""
    source = {"records": {"columns": {"value": {"data_type": data_type, **hints}}}}
    document = DltReader().read(source)["records"]
    value = deepcopy(DLT_EXPECTED[data_type])
    precision = {
        "absent": {},
        "precision": {"precision": hints.get("precision")},
        "scale": {"precision": hints.get("precision"), "scale": hints.get("scale")},
    }
    precision_key = (
        "scale" if hints.get("scale") is not None else "precision" if hints.get("precision") is not None else "absent"
    )
    if precision[precision_key]:
        value.setdefault("x-dlt", {}).update(precision[precision_key])
    if "description" in hints:
        value["description"] = hints["description"]
    if preserve and "timezone" in hints:
        value.setdefault("x-dlt", {}).update({"timezone": False, "unique": False})
    expected = {
        "type": "object",
        "properties": {"value": {"anyOf": [value, {"type": "null"}]} if value else {}},
        "additionalProperties": False,
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:dlt:records",
        "title": "records",
    }
    assert JsonSchemaEmitter(preserve_unknown_hints=preserve).emit(document) == expected


@pytest.mark.parametrize("required", [True, False], ids=["required", "optional"])
@pytest.mark.parametrize(
    "default",
    [None, Decimal(0), Decimal("1.2300000000000000001")],
    ids=["null", "zero", "exact-decimal"],
)
def test_emit_pass_iceberg_decimal_defaults(required: bool, default: Decimal | None) -> None:
    """Non-null Decimal defaults remain Decimal objects at the field metadata scope."""
    schema = Schema(
        NestedField(1, "amount", DecimalType(38, 19), required=required, initial_default=default, write_default=default)
    )
    document = IcebergReader().read(schema)
    result = JsonSchemaEmitter().emit(document)
    scalar = {
        "type": "string",
        "pattern": r"^-?\d{1,19}(\.\d{1,19})?$",
        "x-iceberg": {"logical_type": "decimal", "precision": 38, "scale": 19},
    }
    value = scalar if required else {"anyOf": [scalar, {"type": "null"}]}
    metadata = {"field_id": 1, "required": required}
    if default is not None:
        metadata.update({"initial_default": default, "write_default": default})
    value.setdefault("x-iceberg", {}).update(metadata)
    assert result == {
        "type": "object",
        "properties": {"amount": value},
        "additionalProperties": False,
        **({"required": ["amount"]} if required else {}),
        "$schema": JSON_SCHEMA_DIALECT,
        "x-iceberg": {"schema_id": schema.schema_id, "identifier_field_ids": []},
    }
    if default is not None:
        assert isinstance(result["properties"]["amount"]["x-iceberg"]["initial_default"], Decimal)


@pytest.mark.parametrize("preserve", [True, False], ids=["keep-hints", "drop-hints"])
def test_emit_pass_registered_dlt_extensions(preserve: bool) -> None:
    """Registered vendor hints use their exact names and independently owned payloads."""
    registry = DEFAULT_EXTENSIONS.extend_payload(
        "x-dlt", {"x-vendor-hint": {"type": "array", "items": {"type": "integer"}}}
    )
    source = {"records": {"columns": {"value": {"data_type": "timestamp", "x-vendor-hint": [1, 2]}}}}
    document = DltReader(extension_registry=registry).read(source)["records"]
    emitter = JsonSchemaEmitter(preserve_unknown_hints=preserve)
    hints = {"data_type": "timestamp"}
    if preserve:
        hints["x-vendor-hint"] = [1, 2]
    expected = {
        "type": "object",
        "properties": {
            "value": {
                "anyOf": [
                    {"type": "string", "format": "date-time", "x-dlt": hints},
                    {"type": "null"},
                ]
            }
        },
        "additionalProperties": False,
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:dlt:records",
        "title": "records",
    }
    assert DltToJSONSchema(preserve_unknown_hints=preserve, extension_registry=registry).convert(source) == {
        "records": expected,
    }
    result = emitter.emit(document)
    assert result == expected
    result["properties"]["value"]["anyOf"][0]["x-dlt"]["x-vendor-hint"] = [3]
    assert emitter.emit(document) == expected
    assert source["records"]["columns"]["value"]["x-vendor-hint"] == [1, 2]


def test_emit_pass_json_registered_extensions_and_field_scopes() -> None:
    """Explicit extension contracts preserve nested literal data and field annotations."""
    registry = DEFAULT_EXTENSIONS.register(ExtensionSpec("x-custom", {"type": ["array", "null"]}))
    source = {"$schema": JSON_SCHEMA_DIALECT, "x-custom": [{"$ref": "literal"}], "properties": {"value": False}}
    document = JsonSchemaReader(extension_registry=registry).read(source)
    assert JsonSchemaEmitter().emit(document) == source
    field = replace(
        document.root.properties[0],
        annotations={"description": "denied"},
        extensions=Extensions({"x-custom": None}, registry),
    )
    document = replace(document, root=replace(document.root, properties=(field,)))
    expected = {**source, "properties": {"value": {"not": {}, "description": "denied", "x-custom": None}}}
    emitter = JsonSchemaEmitter()
    result = emitter.emit(document)
    assert result == expected
    result["x-custom"][0]["$ref"] = "changed"
    assert emitter.emit(document) == expected
    assert source["x-custom"] == [{"$ref": "literal"}]


@pytest.mark.parametrize(
    "payload",
    [{"x-unregistered": True}, {"x-dlt": {"unknown_hint": True}}, {"x-pii": "wrong"}],
    ids=["unknown-namespace", "unknown-payload", "invalid-payload"],
)
def test_emit_fail_invalid_extension_contract(payload: dict[str, Any]) -> None:
    """Malformed extension payloads fail at the IR boundary without silent registration."""
    with pytest.raises(ExtensionError):
        JsonSchemaEmitter().emit(SchemaDocument(TypedNode("any", extensions=payload)))


def test_emit_fail_unknown_dlt_type() -> None:
    """Unknown source types raise the shared conversion error for facade translation."""
    document = DltReader().read({"records": {"columns": {"value": {"data_type": "hypercube"}}}})["records"]
    with pytest.raises(ConversionError, match="Column 'value' has unknown dlt data_type 'hypercube'"):
        JsonSchemaEmitter().emit(document)


@pytest.mark.parametrize(
    "field_type",
    [BinaryType(), TimestampNanoType(), TimestamptzNanoType(), UnknownType()],
    ids=["binary", "timestamp-nano", "timestamptz-nano", "unknown"],
)
def test_emit_fail_unsupported_iceberg_types(field_type: IcebergType) -> None:
    """Binary emission and unsupported nanosecond types retain explicit unsupported errors."""
    with pytest.raises(ConversionError, match=r"No (JSON Schema|IR) mapping for Iceberg type"):
        JsonSchemaEmitter().emit(IcebergReader().read(Schema(NestedField(1, "value", field_type))))


@pytest.mark.parametrize(
    "node",
    [
        TypedNode("unknown"),
        TypedNode("map"),
        TypedNode("map", provenance=Provenance("iceberg")),
        TypedNode("number", hints=NodeHints(logical_type="decimal"), provenance=Provenance("iceberg")),
        TypedNode("string", hints=NodeHints(logical_type="fixed"), provenance=Provenance("iceberg")),
        TypedNode("unknown", hints=NodeHints(logical_type="timestamp_ns"), provenance=Provenance("iceberg")),
    ],
    ids=[
        "unknown",
        "incomplete-map",
        "incomplete-iceberg-map",
        "incomplete-decimal",
        "incomplete-fixed",
        "iceberg-nano",
    ],
)
def test_emit_fail_incomplete_type_facts(node: TypedNode) -> None:
    """Missing required logical facts and unknown types are not silently guessed."""
    with pytest.raises(ConversionError):
        JsonSchemaEmitter().emit(SchemaDocument(node))


@pytest.mark.parametrize("policy", [None, 0, "true"], ids=["none", "integer", "string"])
def test_json_schema_emitter_fail_nonboolean_policy(policy: object) -> None:
    """The hint policy accepts actual booleans only."""
    with pytest.raises(TypeError, match="must be a boolean"):
        JsonSchemaEmitter(preserve_unknown_hints=cast("bool", policy))


@pytest.mark.parametrize("precision", [1, 2, 10, 38], ids=["one-digit", "two-digits", "ten-digits", "max-precision"])
def test_decimal_pattern_pass_fixed_patterns(precision: int) -> None:
    """The public decimal helper produces consistent boundary patterns."""
    assert decimal_pattern(precision, 0) == rf"^-?\d{{1,{precision}}}$"
    assert decimal_pattern(precision, precision) == rf"^-?0(\.\d{{1,{precision}}})?$"


def test_emit_pass_typed_tree_over_provenance() -> None:
    """Rendering consumes edited typed facts even when source metadata describes different values."""
    source = {"records": {"columns": {"value": {"data_type": "text", "nullable": False}}}}
    document = DltReader().read(source)["records"]
    original_field = document.root.properties[0]
    edited = replace(
        original_field, node=replace(original_field.node, type="boolean", hints=NodeHints(logical_type="bool"))
    )
    document = replace(document, root=replace(document.root, properties=(edited,)))
    assert JsonSchemaEmitter().emit(document)["properties"] == {"value": {"type": "boolean"}}
    assert source["records"]["columns"]["value"]["data_type"] == "text"


@pytest.mark.parametrize("logical", ["decimal", "timestamp", "wei"], ids=["decimal", "timestamp", "wei"])
def test_emit_pass_sibling_and_call_isolation(logical: str) -> None:
    """Nested output mutations cannot affect siblings, typed inputs or subsequent emissions."""
    column = {"data_type": logical, "nullable": False}
    source = {"records": {"columns": {"hinted": {**column, "scale": 2, "primary_key": True}, "plain": column}}}
    original = deepcopy(source)
    document = DltReader().read(source)["records"]
    emitter = JsonSchemaEmitter()
    plain = deepcopy(DLT_EXPECTED[logical])
    hinted = deepcopy(plain)
    hinted["x-dlt"].update({"precision": None, "scale": 2, "primary_key": True})
    expected = {
        "type": "object",
        "properties": {"hinted": hinted, "plain": plain},
        "required": ["hinted", "plain"],
        "additionalProperties": False,
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:dlt:records",
        "title": "records",
    }
    result = emitter.emit(document)
    result["properties"]["hinted"]["x-dlt"]["data_type"] = "changed"
    assert result["properties"]["plain"] == expected["properties"]["plain"]
    assert emitter.emit(document) == expected
    assert JsonSchemaEmitter().emit(document) == expected
    assert source == original


def test_emit_pass_document_and_type_metadata() -> None:
    """Direct documents retain envelope overrides, type metadata and inferred field requirements."""
    emitter = JsonSchemaEmitter()
    assert emitter.emit(SchemaDocument(TypedNode("any"), name="records")) == {"title": "records"}
    root = TypedNode(
        "object", properties=(Field("flag", TypedNode("boolean", hints=NodeHints(logical_type="bool")), required=True),)
    )
    assert emitter.emit(SchemaDocument(root, name="records", provenance=Provenance("dlt"))) == {
        "type": "object",
        "properties": {"flag": {"type": "boolean"}},
        "required": ["flag"],
        "additionalProperties": False,
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:dlt:records",
        "title": "records",
    }
    iceberg = TypedNode(
        "integer",
        hints=NodeHints(logical_type="int", bit_width=32),
        nullable=True,
        extensions={"x-iceberg": {"logical_type": "int"}},
        provenance=Provenance("iceberg"),
    )
    document = SchemaDocument(iceberg, dialect="urn:dialect", identifier="urn:identifier", name="records")
    assert emitter.emit(document) == {
        "anyOf": [
            {"type": "integer", "minimum": -2147483648, "maximum": 2147483647, "x-iceberg": {"logical_type": "int"}},
            {"type": "null"},
        ],
        "$schema": "urn:dialect",
        "$id": "urn:identifier",
        "title": "records",
    }


@pytest.mark.parametrize("boolean", [True, False], ids=["true", "false"])
def test_emit_pass_boolean_node_metadata_and_edits(boolean: bool) -> None:
    """Boolean nodes become object schemas when annotations or structural facts are added."""
    document = JsonSchemaReader().read({"$schema": JSON_SCHEMA_DIALECT, "properties": {"value": boolean}})
    node = document.root.properties[0].node
    expected: dict[str, Any] = {} if boolean else {"not": {}}
    emitter = JsonSchemaEmitter()
    assert emitter.emit(SchemaDocument(node)) == expected
    assert emitter.emit(SchemaDocument(replace(node, annotations={"description": "value"}))) == {
        "description": "value",
        **expected,
    }
    if boolean:
        node = replace(node, properties=(Field("flag", TypedNode("boolean")),))
        assert emitter.emit(SchemaDocument(node)) == {"properties": {"flag": {"type": "boolean"}}}
