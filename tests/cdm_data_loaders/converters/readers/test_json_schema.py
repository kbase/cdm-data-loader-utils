"""JSON Schema reader facts remain independent of target policy."""

from decimal import Decimal
from typing import Any

import pytest
from jsonschema import (
    Draft3Validator,
    Draft4Validator,
    Draft6Validator,
    Draft7Validator,
    Draft201909Validator,
    Draft202012Validator,
)

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionSpec
from cdm_data_loaders.converters.core.ir import Field, Provenance, TypedNode
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader

DIALECT = "https://json-schema.org/draft/2020-12/schema"


@pytest.mark.parametrize(
    ("dialect", "validator_cls"),
    [
        ("http://json-schema.org/draft-03/schema#", Draft3Validator),
        ("http://json-schema.org/draft-04/schema#", Draft4Validator),
        ("http://json-schema.org/draft-06/schema#", Draft6Validator),
        ("http://json-schema.org/draft-07/schema#", Draft7Validator),
        ("https://json-schema.org/draft/2019-09/schema", Draft201909Validator),
        (DIALECT, Draft202012Validator),
    ],
    ids=["draft3", "draft4", "draft6", "draft7", "draft2019", "draft2020"],
)
def test_read_pass_detected_validator(dialect: str, validator_cls: type) -> None:
    """Retain the detected validator alongside the original dialect URI."""
    document = JsonSchemaReader().read({"$schema": dialect, "type": "object"})
    assert document.dialect == dialect
    assert document.jsonschema_validator_cls is validator_cls


def test_read_pass_exact_tree() -> None:
    """Required presence and nullable type declarations are separate ordered facts."""
    document = JsonSchemaReader().read(
        {
            "$schema": DIALECT,
            "type": "object",
            "required": ["a", "not-declared"],
            "properties": {"a": {"type": ["string", "null"]}, "b": False},
            "additionalProperties": True,
        }
    )
    assert document.root == TypedNode(
        type="object",
        nullable=False,
        declared_type="object",
        inferred_type="object",
        properties=(
            Field(
                "a",
                TypedNode(
                    type="unknown",
                    nullable=True,
                    declared_type=("string", "null"),
                    source_keywords=("type",),
                    provenance=Provenance("json-schema", ("properties", "a")),
                ),
                required=True,
            ),
            Field(
                "b", TypedNode(type="never", nullable=False, provenance=Provenance("json-schema", ("properties", "b")))
            ),
        ),
        required_names=("a", "not-declared"),
        additional_properties=TypedNode(
            type="any", nullable=True, provenance=Provenance("json-schema", ("additionalProperties",))
        ),
        annotations={"$schema": DIALECT},
        provenance=Provenance("json-schema"),
        source_keywords=("$schema", "type", "required", "properties", "additionalProperties"),
    )


def test_read_pass_constraints_branches_and_literal_data() -> None:
    """No schema traversal occurs inside enum, const, defaults, examples or extensions."""
    literal = {"$ref": "literal", "allOf": ["data"], "x-unregistered": None}
    constraints = {
        "minimum": Decimal("-2.5"),
        "maximum": 4,
        "exclusiveMinimum": -3,
        "exclusiveMaximum": 5,
        "multipleOf": Decimal("0.0000000000000000001"),
        "enum": [literal, None],
        "const": literal,
        "pattern": "^[a-z]$",
        "format": "custom",
        "contentEncoding": "base64",
        "minLength": 1,
    }
    schema = {
        "$schema": DIALECT,
        "properties": {
            "data": {
                **constraints,
                "default": literal,
                "examples": [literal],
                "anyOf": [True, {"type": "null"}],
                "oneOf": [{"type": "integer"}, False],
            }
        },
    }
    node = JsonSchemaReader().read(schema).root.properties[0].node
    assert node.type == "unknown"
    assert node.inferred_type == "string"
    assert node.constraints["multipleOf"] == Decimal("0.0000000000000000001")
    assert node.constraints["enum"] == ({"$ref": "literal", "allOf": ("data",), "x-unregistered": None}, None)
    assert node.any_of[0].type == "any"
    assert node.any_of[1].declared_type == "null"
    assert tuple(branch.type for branch in node.one_of) == ("integer", "never")
    assert node.annotations["default"] == node.constraints["const"]


def test_read_pass_tuple_dynamic_and_schema_keywords() -> None:
    """Tuple forms, dynamic properties and schema-valued keywords remain structured."""
    schema = {
        "$schema": DIALECT,
        "required": [],
        "properties": {
            "tuple": {
                "items": [True, False],
                "prefixItems": [{"type": "integer"}],
                "additionalItems": False,
                "contains": {"type": "string"},
            }
        },
        "patternProperties": {"^first": True, "^second": {"type": "number"}},
        "additionalProperties": {"type": "boolean"},
        "dependentSchemas": {"trigger": {"required": ["other"]}},
        "dependencies": {"literal": ["other"], "schema": {"properties": {"nested": True}}},
        "$defs": {"saved": {"type": "null"}},
        "not": False,
    }
    root = JsonSchemaReader().read(schema).root
    assert root.required_names == ()
    assert tuple(root.pattern_properties) == ("^first", "^second")
    assert root.additional_properties.type == "boolean"
    node = root.properties[0].node
    assert tuple(child.type for child in node.items) == ("any", "never")
    assert node.prefix_items[0].type == "integer"
    assert node.schema_keywords["additionalItems"].type == "never"
    assert node.schema_keywords["contains"].type == "string"
    assert root.schema_maps["dependencies"]["schema"].properties[0].name == "nested"
    assert root.constraints["dependencies"] == {"literal": ("other",)}
    assert root.schema_maps["dependentSchemas"]["trigger"].required_names == ("other",)
    assert root.schema_maps["$defs"]["saved"].type == "null"
    assert root.schema_keywords["not"].type == "never"


@pytest.mark.parametrize(
    "fragment",
    [
        {"$ref": "#/somewhere"},
        {"allOf": []},
        {"$dynamicRef": "#node"},
        {"$recursiveRef": "#node"},
        {"items": 3},
        {"not": {"$ref": "bad"}},
    ],
    ids=["ref", "all-of", "dynamic-ref", "recursive-ref", "invalid-items", "nested-ref"],
)
def test_read_fail_unresolved_or_invalid(fragment: dict[str, Any]) -> None:
    """Defensive schema guards reject unresolved and malformed schema children."""
    with pytest.raises(ConversionError):
        JsonSchemaReader().read({"$schema": DIALECT, "properties": {"bad": fragment}})


@pytest.mark.parametrize(
    "schema", [{}, {"$schema": DIALECT, "type": "array"}, True], ids=["dialect", "root-type", "boolean-root"]
)
def test_read_fail_document(schema: dict[str, Any] | bool) -> None:
    """Document inputs require a dialect and an object-compatible root."""
    with pytest.raises(ConversionError):
        JsonSchemaReader().read(schema)


def test_read_pass_custom_registry_and_empty_schemas() -> None:
    """Custom extensions and empty object schemas need no nonempty required list."""
    registry = DEFAULT_EXTENSIONS.register(ExtensionSpec("x-owned", {"type": ["array", "null"]}))
    root = (
        JsonSchemaReader(registry)
        .read(
            {
                "$schema": DIALECT,
                "x-owned": None,
                "properties": {
                    "any": True,
                    "unknown": {},
                    "empty": {"anyOf": [], "prefixItems": []},
                },
            }
        )
        .root
    )
    assert root.extensions["x-owned"] is None
    assert root.properties[1].node.type == "unknown"
    assert root.properties[2].node.any_of == ()
    assert root.properties[2].node.prefix_items == ()
    assert JsonSchemaReader().read({"$schema": DIALECT, "type": ["object", "null"]}).root.nullable is True
