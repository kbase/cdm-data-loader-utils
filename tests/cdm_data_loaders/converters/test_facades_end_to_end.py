"""Cross-component conversion contracts without services or Spark sessions."""

import json
from collections import UserDict
from dataclasses import FrozenInstanceError
from decimal import Decimal
from pathlib import Path
from typing import Any, Final

import pytest
from frozendict import frozendict
from jsonschema import Draft202012Validator
from pydantic import ValidationError
from pyiceberg.schema import Schema
from pyiceberg.types import DecimalType, NestedField
from pyspark.sql.types import StringType

from cdm_data_loaders.converters import jsonschema_to_pyspark as spark_facade
from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionError, ExtensionSpec
from cdm_data_loaders.converters.core.ir import Field, TypedNode
from cdm_data_loaders.converters.dlt_to_jsonschema import DltToJSONSchema, DltToJSONSchemaError
from cdm_data_loaders.converters.emitters import pyspark as spark_emitter
from cdm_data_loaders.converters.emitters.dlt import DltEmitter
from cdm_data_loaders.converters.emitters.json_schema import JSON_SCHEMA_DIALECT, JsonSchemaEmitter
from cdm_data_loaders.converters.emitters.pyspark import PySparkEmitter
from cdm_data_loaders.converters.jsonschema_to_dlt import InvalidJSONSchemaError as InvalidDltSchemaError
from cdm_data_loaders.converters.jsonschema_to_dlt import JSONSchemaToDlt, JSONSchemaToDltError
from cdm_data_loaders.converters.jsonschema_to_pyspark import (
    InvalidJSONSchemaError as InvalidSparkSchemaError,
)
from cdm_data_loaders.converters.jsonschema_to_pyspark import JSONSchemaToPySpark, JSONSchemaToPySparkError
from cdm_data_loaders.converters.pyiceberg_to_jsonschema import (
    convert_field,
    convert_type,
    table_to_json_schema,
)
from cdm_data_loaders.converters.readers.dlt import DltReader
from cdm_data_loaders.converters.readers.iceberg import IcebergReader
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader
from tests.cdm_data_loaders.converters.conftest import make_table

GOLDEN_PATH: Final = Path(__file__).parents[2] / "data/converters/ir/legacy_outputs.json"
CASES: Final[list[dict[str, Any]]] = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_facades_end_to_end_pass_captured_outputs(case: dict[str, Any]) -> None:
    """All four public routes and their direct compositions match captured independent goldens."""
    source, options, expected = case["input"], case["options"], case["expected"]
    match case["route"]:
        case "json-dlt":
            assert JSONSchemaToDlt(**options).convert(source) == expected
            assert DltEmitter(**options).emit(JsonSchemaReader().read(source)) == expected
        case "json-spark":
            assert JSONSchemaToPySpark(**options).convert(source).jsonValue() == expected
            assert PySparkEmitter(**options).emit(JsonSchemaReader().read(source)).jsonValue() == expected
        case "dlt-json":
            assert DltToJSONSchema(**options).convert(source) == expected
            reader_options = {key: value for key, value in options.items() if key != "preserve_unknown_hints"}
            emitter = JsonSchemaEmitter(preserve_unknown_hints=options.get("preserve_unknown_hints", True))
            assert {
                name: emitter.emit(doc) for name, doc in DltReader(**reader_options).read(source).items()
            } == expected
        case "iceberg-json":
            identifier = tuple(source["identifier"])
            table = make_table(identifier, Schema.model_validate(source["schema"]), properties=source["properties"])
            assert table_to_json_schema(table, identifier) == expected
            assert JsonSchemaEmitter().emit(IcebergReader().read_table(table, identifier)) == expected
        case _:
            pytest.fail(f"Unhandled golden route: {case['route']}")


@pytest.mark.parametrize("converter_type", [JSONSchemaToPySpark, PySparkEmitter], ids=["facade", "direct-emitter"])
@pytest.mark.parametrize(
    ("options", "date_type", "datetime_type"),
    [
        ({"format_map": frozendict()}, "string", "string"),
        ({"format_map": {"date-time": StringType()}}, "date", "string"),
        ({}, "date", "timestamp"),
    ],
    ids=["empty-frozen-map", "override-keeps-known-defaults", "omitted-map-uses-builtins"],
)
def test_facades_end_to_end_pass_spark_format_map(
    converter_type: type[JSONSchemaToPySpark] | type[PySparkEmitter],
    options: dict[str, object],
    date_type: str,
    datetime_type: str,
) -> None:
    """Facade and direct emission preserve explicit empty maps, overrides and default formats."""
    source = {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "object",
        "properties": {
            "date": {"type": "string", "format": "date"},
            "datetime": {"type": "string", "format": "date-time"},
        },
    }
    converter = converter_type.model_validate(options)
    result = (
        converter.convert(source)
        if isinstance(converter, JSONSchemaToPySpark)
        else converter.emit(JsonSchemaReader().read(source))
    )
    assert result.jsonValue() == {
        "type": "struct",
        "fields": [
            {"name": "date", "type": date_type, "nullable": True, "metadata": {}},
            {"name": "datetime", "type": datetime_type, "nullable": True, "metadata": {}},
        ],
    }


@pytest.mark.parametrize("converter_type", [JSONSchemaToPySpark, PySparkEmitter], ids=["facade", "direct-emitter"])
@pytest.mark.parametrize("value", [{}, UserDict()], ids=["empty-dict", "empty-user-dict"])
def test_facades_end_to_end_fail_spark_empty_mutable_format_map(
    converter_type: type[JSONSchemaToPySpark] | type[PySparkEmitter], value: object
) -> None:
    """Empty mutable maps fail the frozendict field requirement for both Spark entry points."""
    with pytest.raises(ValidationError) as error:
        converter_type.model_validate({"format_map": value})
    assert [(item["loc"], item["type"]) for item in error.value.errors()] == [(("format_map",), "is_instance_of")]


def test_facades_end_to_end_pass_supported_json_dlt_json() -> None:
    """The supported required scalar subset retains structure through JSON, dlt and JSON."""
    source = {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:dlt:records",
        "title": "records",
        "type": "object",
        "properties": {"id": {"type": "integer"}, "name": {"type": "string"}},
        "required": ["id", "name"],
        "additionalProperties": False,
    }
    stored = JSONSchemaToDlt(schema_name="records").convert(source)
    assert stored == {
        "name": "records",
        "version": 1,
        "engine_version": 11,
        "previous_hashes": [],
        "tables": {
            "records": {
                "name": "records",
                "columns": {
                    "id": {"name": "id", "data_type": "bigint", "nullable": False},
                    "name": {"name": "name", "data_type": "text", "nullable": False},
                },
                "write_disposition": "append",
                "resource": "records",
            }
        },
    }
    assert DltToJSONSchema().convert(stored) == {"records": source}


@pytest.mark.parametrize(
    ("hints", "spark_type", "json_hints"),
    [
        ({}, "decimal(38,0)", {"data_type": "decimal"}),
        ({"precision": None, "scale": None}, "decimal(38,0)", {"data_type": "decimal"}),
        ({"precision": 20, "scale": 8}, "decimal(20,8)", {"data_type": "decimal", "precision": 20, "scale": 8}),
        ({"scale": 2}, "decimal(38,2)", {"data_type": "decimal", "precision": None, "scale": 2}),
        ({"precision": 10, "scale": 0}, "decimal(10,0)", {"data_type": "decimal", "precision": 10, "scale": 0}),
    ],
    ids=["absent", "null", "precise", "scale-only", "zero-scale"],
)
def test_facades_end_to_end_pass_dlt_spark_precision(
    hints: dict[str, Any], spark_type: str, json_hints: dict[str, Any]
) -> None:
    """Direct dlt-to-Spark composition preserves numeric precision instead of decimal string encoding."""
    source = {"records": {"columns": {"amount": {"data_type": "decimal", "nullable": False, **hints}}}}
    document = DltReader().read(source)["records"]
    assert PySparkEmitter().emit(document).jsonValue() == {
        "type": "struct",
        "fields": [{"name": "amount", "type": spark_type, "nullable": False, "metadata": {}}],
    }
    assert DltToJSONSchema().convert(source)["records"] == {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": "urn:dlt:records",
        "title": "records",
        "type": "object",
        "properties": {
            "amount": {
                "type": "string",
                "pattern": r"^-?\d+(\.\d+)?$",
                "x-dlt": json_hints,
            }
        },
        "required": ["amount"],
        "additionalProperties": False,
    }


@pytest.mark.parametrize("preserve", [True, False], ids=["keep", "drop"])
def test_facades_end_to_end_pass_registered_dlt_hints(preserve: bool) -> None:
    """Custom dlt hints need registration even when emission filters them out."""
    source = {"records": {"columns": {"id": {"data_type": "bigint", "nullable": False, "x-owner": "lab"}}}}
    with pytest.raises(DltToJSONSchemaError, match="Invalid x-dlt") as caught:
        DltToJSONSchema(preserve_unknown_hints=preserve).convert(source)
    assert isinstance(caught.value.__cause__, ExtensionError)
    assert caught.value.__cause__.__cause__ is not None
    registry = DEFAULT_EXTENSIONS.extend_payload("x-dlt", {"x-owner": {"type": "string", "enum": ["lab"]}})
    result = DltToJSONSchema(preserve_unknown_hints=preserve, extension_registry=registry).convert(source)
    expected = {"type": "integer", **({"x-dlt": {"x-owner": "lab"}} if preserve else {})}
    assert result["records"]["properties"] == {"id": expected}
    source["records"]["columns"]["id"]["x-owner"] = "unregistered-value"
    with pytest.raises(DltToJSONSchemaError, match="Invalid x-dlt"):
        DltToJSONSchema(preserve_unknown_hints=preserve, extension_registry=registry).convert(source)


@pytest.mark.parametrize("route", ["dlt", "spark"], ids=["json-dlt", "json-spark"])
def test_facades_end_to_end_pass_registered_namespace(route: str) -> None:
    """Facade constructors retain immutable explicit contracts and reject arbitrary keys and values."""
    contract: dict[str, Any] = {
        "type": "object",
        "properties": {"classification": {"enum": ["public"]}},
        "required": ["classification"],
        "additionalProperties": False,
    }
    registry = DEFAULT_EXTENSIONS.register(ExtensionSpec("x-vendor", contract))
    contract["properties"].clear()
    converter = (
        JSONSchemaToDlt(schema_name="records", extension_registry=registry)
        if route == "dlt"
        else JSONSchemaToPySpark(extension_registry=registry, extra_metadata_keywords=frozenset({"x-vendor"}))
    )
    source: dict[str, Any] = {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "object",
        "properties": {
            "name": {"type": "string", "x-vendor": {"classification": "public"}},
        },
    }
    result = converter.convert(source)
    if route == "dlt":
        assert result["tables"]["records"]["columns"] == {
            "name": {"name": "name", "data_type": "text", "nullable": True},
        }
    else:
        assert result.jsonValue() == {
            "type": "struct",
            "fields": [
                {
                    "name": "name",
                    "type": "string",
                    "nullable": True,
                    "metadata": {"jsonschema": {"x-vendor": {"classification": "public"}}},
                }
            ],
        }
    assert converter.extension_registry is registry
    with pytest.raises(FrozenInstanceError):
        registry.specs = ()
    with pytest.raises(ValidationError, match="frozen"):
        converter.extension_registry = DEFAULT_EXTENSIONS
    default = JSONSchemaToDlt(schema_name="records") if route == "dlt" else JSONSchemaToPySpark()
    with pytest.raises(ConversionError, match="Unregistered extension") as caught:
        default.convert(source)
    assert isinstance(caught.value.__cause__, ExtensionError)
    for payload in ({"classification": "private"}, {"classification": "public", "extra": True}):
        source["properties"]["name"]["x-vendor"] = payload
        with pytest.raises(ConversionError, match="Invalid x-vendor") as caught:
            converter.convert(source)
        assert isinstance(caught.value.__cause__, ExtensionError)
        assert caught.value.__cause__.__cause__ is not None


@pytest.mark.parametrize(
    "extensions",
    [{"x-pii": "yes"}, {"x-vendor": True}, {"x-dlt": {"custom": True}}, {"x-xsv-config": {"typo": True}}],
    ids=["invalid-builtin", "unregistered", "unknown-dlt-hint", "invalid-xsv"],
)
def test_facades_end_to_end_fail_invalid_extensions(extensions: dict[str, Any]) -> None:
    """Both forward facades validate extension data even when metadata is not selected."""
    source = {"$schema": JSON_SCHEMA_DIALECT, "properties": {"value": {"type": "string", **extensions}}}
    for converter, error_type in (
        (JSONSchemaToDlt(schema_name="records"), JSONSchemaToDltError),
        (JSONSchemaToPySpark(), JSONSchemaToPySparkError),
    ):
        with pytest.raises(error_type) as caught:
            converter.convert(source)
        assert isinstance(caught.value.__cause__, ExtensionError)


@pytest.mark.parametrize("route", ["dlt", "spark"], ids=["dlt", "spark"])
def test_facades_end_to_end_fail_exact_root_guards(route: str) -> None:
    """Public root messages and missing-dialect exception identities remain stable."""
    converter = JSONSchemaToDlt(schema_name="records") if route == "dlt" else JSONSchemaToPySpark()
    error_type = JSONSchemaToDltError if route == "dlt" else JSONSchemaToPySparkError
    missing_type = InvalidDltSchemaError if route == "dlt" else InvalidSparkSchemaError
    target = "a dlt table" if route == "dlt" else "a StructType"
    name = "JSONSchemaToDlt" if route == "dlt" else "JSONSchemaToPySpark"
    with pytest.raises(missing_type) as caught:
        converter.convert({"type": "object"})
    assert str(caught.value) == (
        f"Input JSON Schema is missing a '$schema' keyword. {name} requires schemas "
        "to explicitly declare their dialect via '$schema'; it will not assume a default."
    )
    with pytest.raises(error_type) as caught:
        converter.convert({"$schema": JSON_SCHEMA_DIALECT, "type": ["object", "null"], "x-unregistered": True})
    assert str(caught.value) == (f"Root schema must be of type 'object' to map to {target}, got: ['object', 'null']")
    with pytest.raises(json.JSONDecodeError):
        converter.convert_from_string("type: object")


@pytest.mark.parametrize("keyword", ["$ref", "$dynamicRef", "$recursiveRef"], ids=["ref", "dynamic", "recursive"])
def test_facades_end_to_end_fail_references_in_ignored_branches(keyword: str) -> None:
    """Reader validation rejects unresolved schema branches even when target dispatch would ignore them."""
    source = {
        "$schema": JSON_SCHEMA_DIALECT,
        "properties": {
            "value": {
                "type": "string",
                "oneOf": [{"type": "string"}, {keyword: "#/$defs/hidden"}],
            }
        },
    }
    for converter, error_type in (
        (JSONSchemaToDlt(schema_name="records"), JSONSchemaToDltError),
        (JSONSchemaToPySpark(), JSONSchemaToPySparkError),
    ):
        with pytest.raises(error_type, match=r"[Uu]nresolved"):
            converter.convert(source)


@pytest.mark.parametrize("required", [True, False], ids=["required", "optional"])
@pytest.mark.parametrize("default", [None, Decimal("0.1000000000000000001")], ids=["null-default", "exact-default"])
def test_facades_end_to_end_pass_iceberg_fragments(required: bool, default: Decimal | None) -> None:
    """Public fragment APIs preserve exact decimal metadata and optional wrapper placement."""
    source = NestedField(1, "amount", DecimalType(38, 19), required=required, initial_default=default, doc="amount")
    scalar = {
        "type": "string",
        "pattern": r"^-?\d{1,19}(\.\d{1,19})?$",
        "x-iceberg": {"logical_type": "decimal", "precision": 38, "scale": 19},
    }
    metadata: dict[str, Any] = {"field_id": 1, "required": required}
    if default is not None:
        metadata["initial_default"] = default
    expected = (
        {**scalar, "x-iceberg": {**scalar["x-iceberg"], **metadata}}
        if required
        else {
            "anyOf": [scalar, {"type": "null"}],
            "x-iceberg": metadata,
        }
    )
    expected["description"] = "amount"
    assert convert_type(source.field_type) == scalar
    assert JsonSchemaEmitter().emit_node(IcebergReader().read_node(source.field_type)) == scalar
    assert convert_field(source) == expected
    assert JsonSchemaEmitter().emit_field(IcebergReader().read_field(source)) == expected


def test_facades_end_to_end_pass_fragment_sources_and_aliases() -> None:
    """Fragment APIs select source policy and facade helper exports retain emitter ownership."""
    assert spark_facade.ConversionContext is spark_emitter.ConversionContext
    emitter = JsonSchemaEmitter()
    for source in (True, False, {"type": "integer"}):
        node = JsonSchemaReader().read_node(source)
        assert emitter.emit_node(node) == source
        assert emitter.emit_field(Field("value", node)) == source
    node = TypedNode(type="string")
    assert emitter.emit_node(node) == {"type": "string"}
    assert emitter.emit_field(Field("value", node)) == {"type": "string"}
    field = (
        DltReader()
        .read({"records": {"columns": {"id": {"data_type": "bigint", "nullable": False}}}})["records"]
        .root.properties[0]
    )
    assert emitter.emit_node(field.node) == {"type": "integer"}
    assert emitter.emit_field(field) == {"type": "integer"}
    context = spark_facade.ConversionContext(Draft202012Validator)
    assert JSONSchemaToPySpark()._convert_property({"type": "string", "x-pii": True}, context) == (StringType(), {})  # noqa: SLF001


@pytest.mark.parametrize("value", [None, True, "registry", {"specs": []}], ids=["null", "boolean", "string", "mapping"])
def test_facades_end_to_end_fail_registry_constructor(value: object) -> None:
    """Facade constructors must not silently accept arbitrary registry substitutes."""
    for facade in (JSONSchemaToDlt, JSONSchemaToPySpark, DltToJSONSchema):
        with pytest.raises(ValidationError):
            facade(**({"schema_name": "records"} if facade is JSONSchemaToDlt else {}), extension_registry=value)


@pytest.mark.parametrize(
    ("source", "message"),
    [({"type": "string"}, None), ({"$ref": "#/$defs/value"}, r"unresolved \$ref"), ({"allOf": []}, "unmerged 'allOf'")],
    ids=["resolved", "unresolved-ref", "unmerged-allof"],
)
def test_facades_end_to_end_pass_dlt_reference_bridge(source: dict[str, Any], message: str | None) -> None:
    """The retained dlt reference guard preserves its exception subtype and diagnostics."""
    if message is None:
        assert JSONSchemaToDlt._reject_unresolved_references(source) is None  # noqa: SLF001
    else:
        with pytest.raises(JSONSchemaToDltError, match=message):
            JSONSchemaToDlt._reject_unresolved_references(source)  # noqa: SLF001
