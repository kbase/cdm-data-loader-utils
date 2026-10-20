"""Target-independent schema model contracts."""

from collections.abc import Callable
from dataclasses import FrozenInstanceError
from decimal import Decimal
from functools import partial
from typing import Any, Final
from unittest.mock import patch

import pytest

from cdm_data_loaders.converters.core import ir
from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.ir import Field, NodeHints, Provenance, SchemaDocument, TypedNode
from cdm_data_loaders.converters.core.ir_values import freeze_mapping, freeze_value, mutable_value

_METADATA_POSITIONS: Final = [
    pytest.param(partial(TypedNode, type="any"), "constraints", id="type-constraints"),
    pytest.param(partial(TypedNode, type="any"), "annotations", id="type-annotations"),
    pytest.param(partial(Field, "value", TypedNode(type="any")), "annotations", id="field-annotations"),
    pytest.param(partial(SchemaDocument, TypedNode(type="any")), "annotations", id="document-annotations"),
]


@pytest.mark.parametrize(("constructor", "keyword"), _METADATA_POSITIONS)
@pytest.mark.parametrize(
    ("namespace", "value"),
    [("x-pii", "bad"), ("x-unregistered", "bad"), ("x-pii", True)],
    ids=["invalid-known", "unknown", "valid-known"],
)
def test_ir_metadata_fail_extension_namespace(
    constructor: Callable[..., TypedNode | Field | SchemaDocument], keyword: str, namespace: str, value: str | bool
) -> None:
    """Schema-level extension keys must use the validated extensions mapping."""
    with pytest.raises(ConversionError, match=rf"{namespace}.*extensions"):
        constructor(**{keyword: {namespace: value}})


@pytest.mark.parametrize(("constructor", "keyword"), _METADATA_POSITIONS)
def test_ir_metadata_pass_literal_extension_keys(
    constructor: Callable[..., TypedNode | Field | SchemaDocument], keyword: str
) -> None:
    """Literal data and unknown non-extension schema keywords remain unrestricted."""
    literal = {"x-pii": "bad", "nested": [{"x-unregistered": True}]}
    metadata = {"examples": [literal], "enum": [literal], "const": literal, "default": literal, "custom": literal}
    model = constructor(**{keyword: metadata}, extensions={"x-pii": True})
    assert mutable_value(getattr(model, keyword)) == metadata
    assert dict(model.extensions) == {"x-pii": True}


def test_schema_document_pass_owned_values_and_presence() -> None:
    """Field presence, nullable values and nested metadata remain independent."""
    values = {"enum": [{"amount": Decimal("1.2300")}], "const": None}
    node = TypedNode(type="number", nullable=True, constraints=values)
    document = SchemaDocument(root=TypedNode(type="object", properties=(Field("amount", node, required=True),)))
    values["enum"].clear()
    assert node.constraints == {"enum": ({"amount": Decimal("1.2300")},), "const": None}
    assert document.root.properties == (Field("amount", node, required=True),)
    with pytest.raises(TypeError):
        node.constraints["const"] = 1
    with pytest.raises(FrozenInstanceError):
        node.type = "string"


def test_schema_document_pass_post_init_freezes_annotations_once() -> None:
    """__post_init__ freezes document annotations exactly once, not once per identity field."""
    root = TypedNode(type="any")
    with patch("cdm_data_loaders.converters.core.ir._schema_metadata", wraps=ir._schema_metadata) as spy:  # noqa: SLF001
        SchemaDocument(root=root, name="n", dialect="d", identifier="i", annotations={"a": 1})
    assert spy.call_count == 1


@pytest.mark.parametrize(
    "value", [object(), {1: "bad"}, float("inf"), Decimal("NaN")], ids=["object", "key", "infinity", "nan"]
)
def test_freeze_value_fail_invalid(value: object) -> None:
    """Unsupported objects and non-finite numbers fail during construction."""
    with pytest.raises((TypeError, ValueError)):
        freeze_value(value)


def test_provenance_pass_owned_metadata() -> None:
    """Metadata preserves absent versus explicit null values."""
    metadata = {"description": None, "tags": ["one"]}
    provenance = Provenance("dlt", ("tables", "entries"), metadata)
    metadata["tags"].append("two")
    assert provenance.metadata == {"description": None, "tags": ("one",)}
    assert "default" not in provenance.metadata


@pytest.mark.parametrize(
    "kwargs",
    [
        {"type": "bogus"},
        {"nullable": "yes"},
        {"inferred_type": []},
        {"hints": {}},
        {"declared_type": [1]},
        {"declared_type": 1},
        {"properties": [1]},
        {"properties": [Field("a", TypedNode(type="any")), Field("a", TypedNode(type="any"))]},
        {"required_names": [1]},
        {"items": 1},
        {"items": [1]},
        {"any_of": [1]},
        {"additional_properties": False},
        {"pattern_properties": {"name": {}}},
        {"schema_maps": {1: {}}},
        {"schema_maps": {"definitions": {1: TypedNode(type="any")}}},
        {"provenance": {}},
        {"source_keywords": [1]},
    ],
    ids=[
        "kind",
        "nullable",
        "inferred",
        "hints",
        "union-entry",
        "declaration",
        "property",
        "duplicate",
        "required",
        "item",
        "tuple-item",
        "branch",
        "additional",
        "pattern",
        "map-key",
        "map-value",
        "provenance",
        "keywords",
    ],
)
def test_typed_node_fail_invalid_fields(kwargs: dict[str, Any]) -> None:
    """Normal constructors reject malformed types and mutable model substitutes."""
    with pytest.raises((TypeError, ValueError)):
        TypedNode(**{"type": "any", **kwargs})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"logical_type": 1},
        {"precision": -1},
        {"scale": True},
        {"timezone": 1},
    ],
    ids=["logical-type", "negative-precision", "boolean-scale", "timezone"],
)
def test_node_hints_fail_invalid(kwargs: dict[str, Any]) -> None:
    """Logical hints validate values rather than coercing invalid inputs."""
    with pytest.raises((TypeError, ValueError)):
        NodeHints(**kwargs)


def test_models_fail_invalid_envelopes_and_metadata() -> None:
    """Envelope, field and source validation run during normal construction."""
    with pytest.raises(ValueError, match="source"):
        Provenance("bogus")
    with pytest.raises(TypeError, match="paths"):
        Provenance("dlt", (True,))
    with pytest.raises(TypeError, match="Fields require"):
        Field("name", {})
    with pytest.raises(TypeError, match="root"):
        SchemaDocument(root={})
    with pytest.raises(TypeError, match="name"):
        SchemaDocument(root=TypedNode(type="any"), name=1)
    with pytest.raises(TypeError, match="mapping"):
        freeze_mapping([])


def test_typed_node_pass_all_containers_owned() -> None:
    """All child sequences, schema maps and nested literals are deeply isolated."""
    child = TypedNode(type="null", nullable=True)
    children = [child]
    nested = {"child": child}
    node = TypedNode(
        type="array",
        items=children,
        prefix_items=children,
        any_of=children,
        one_of=children,
        pattern_properties=nested,
        schema_maps={"$defs": nested},
        schema_keywords={"not": child},
        annotations={"default": [1, {"value": Decimal("1.00")}], "null": None},
    )
    children.clear()
    nested.clear()
    assert node.items == node.prefix_items == node.any_of == node.one_of == (child,)
    assert node.schema_maps == {"$defs": {"child": child}}
    assert node.pattern_properties == {"child": child}
    assert node.schema_keywords == {"not": child}
    with pytest.raises(TypeError):
        node.schema_maps["$defs"]["child"] = child
    value = mutable_value(node.annotations)
    assert value == {"default": [1, {"value": Decimal("1.00")}], "null": None}
    value["default"].clear()
    assert node.annotations["default"] == (1, {"value": Decimal("1.00")})
    assert (freeze_value(1.5),) == (1.5,)
