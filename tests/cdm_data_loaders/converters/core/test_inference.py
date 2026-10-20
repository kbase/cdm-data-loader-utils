"""Tests for target-independent type inference."""

from decimal import Decimal
from typing import Any

import pytest

from cdm_data_loaders.converters.core.inference import (
    IMPLICIT_ARRAY_KEYWORDS,
    IMPLICIT_NUMBER_KEYWORDS,
    IMPLICIT_OBJECT_KEYWORDS,
    IMPLICIT_STRING_KEYWORDS,
    decimal_pattern,
    decimal_places,
    infer_implicit_type,
    is_unconstrained_node,
    json_type_from_enum,
    resolve_union_type,
)
from cdm_data_loaders.converters.core.ir import TypedNode


@pytest.mark.parametrize(
    ("node", "expected"),
    [
        pytest.param(TypedNode(type="never"), True, id="never"),
        pytest.param(TypedNode(type="unknown"), True, id="bare-unknown"),
        pytest.param(TypedNode(type="any"), False, id="any"),
        pytest.param(TypedNode(type="string"), False, id="string"),
        pytest.param(TypedNode(type="unknown", declared_type="alien"), False, id="declared-type"),
        pytest.param(TypedNode(type="unknown", inferred_type="string"), False, id="inferred-type"),
        pytest.param(TypedNode(type="unknown", source_keywords=("type",)), False, id="source-keyword"),
        pytest.param(TypedNode(type="unknown", constraints={"const": 1}), False, id="constraint"),
        pytest.param(TypedNode(type="unknown", annotations={"description": "d"}), False, id="annotation"),
        pytest.param(TypedNode(type="unknown", one_of=()), False, id="empty-one-of"),
        pytest.param(TypedNode(type="unknown", any_of=()), False, id="empty-any-of"),
        pytest.param(
            TypedNode(type="unknown", schema_keywords={"not": TypedNode(type="never")}), False, id="schema-keyword"
        ),
        pytest.param(
            TypedNode(type="unknown", schema_maps={"$defs": {"a": TypedNode(type="string")}}), False, id="schema-map"
        ),
    ],
)
def test_is_unconstrained_node_pass(node: TypedNode, expected: bool) -> None:
    """Only a never node or a keyword-free unknown node counts as unconstrained."""
    assert is_unconstrained_node(node) is expected


@pytest.mark.parametrize(
    ("keyword", "expected"),
    [
        pytest.param(keyword, json_type, id=f"{json_type}-{keyword}")
        for json_type, keywords in (
            ("object", IMPLICIT_OBJECT_KEYWORDS),
            ("array", IMPLICIT_ARRAY_KEYWORDS),
            ("string", IMPLICIT_STRING_KEYWORDS),
            ("number", IMPLICIT_NUMBER_KEYWORDS),
        )
        for keyword in sorted(keywords)
    ],
)
def test_infer_implicit_type_pass_keywords(keyword: str, expected: str) -> None:
    """Infer each keyword family from presence, independent of its value."""
    assert infer_implicit_type({keyword: None}) == expected


@pytest.mark.parametrize(
    ("schema", "expected"),
    [
        ({}, None),
        ({"title": "label", "description": "text"}, None),
        ({"const": 1}, None),
        ({"type": "integer"}, None),
        ({"required": [], "items": {}, "pattern": "", "minimum": 0}, "object"),
        ({"items": {}, "pattern": "", "minimum": 0}, "array"),
        ({"pattern": "", "minimum": 0}, "string"),
        ({"multipleOf": 1}, "number"),
    ],
    ids=[
        "empty",
        "annotations",
        "const",
        "explicit-type-only",
        "object-first",
        "array-first",
        "string-first",
        "numeric",
    ],
)
def test_infer_implicit_type_pass_precedence(schema: dict[str, Any], expected: str | None) -> None:
    """Preserve keyword precedence and leave unconstrained fragments uninferred."""
    original = schema.copy()
    assert infer_implicit_type(schema) == expected
    assert schema == original


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([], "string"),
        ([True, False], "boolean"),
        ([-1, 0, 2], "integer"),
        ([0.0, 1.5], "number"),
        ([1, 2.5], "number"),
        ([1.0, 2.0], "number"),
        (["", "text"], "string"),
        ([True, 1], "string"),
        ([False, 1.5], "string"),
        ([1, "text"], "string"),
        ([None], "string"),
        ([1, None], "string"),
        ([{}], "string"),
        ([[]], "string"),
    ],
    ids=[
        "empty",
        "booleans",
        "integers",
        "floats",
        "mixed-numbers",
        "integral-floats",
        "strings",
        "bool-is-not-int",
        "bool-is-not-number",
        "heterogeneous",
        "null",
        "nullable-int",
        "objects",
        "arrays",
    ],
)
def test_json_type_from_enum_pass_scalar_inference(values: list[Any], expected: str) -> None:
    """Infer scalar enum types without treating booleans as integers or numbers."""
    assert json_type_from_enum(values) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (0, 0),
        (100, 0),
        (-10, 0),
        (0.01, 2),
        (-0.125, 3),
        (1e-7, 7),
        (Decimal("1.2300"), 4),
        (Decimal("1E+3"), 0),
        (Decimal("0.000"), 3),
        (Decimal("1.0000000000000000000000000001"), 28),
        (float("inf"), 0),
        (float("-inf"), 0),
        (float("nan"), 0),
        (Decimal("sNaN"), 0),
    ],
    ids=[
        "zero",
        "integer",
        "negative-integer",
        "hundredth",
        "negative-fraction",
        "scientific",
        "trailing-zeros",
        "positive-exponent",
        "scaled-zero",
        "exact-decimal",
        "infinity",
        "negative-infinity",
        "nan",
        "signaling-nan",
    ],
)
def test_decimal_places_pass_scale(value: float | Decimal, expected: int) -> None:
    """Count decimal scale without converting exact decimals through binary floats."""
    assert decimal_places(value) == expected


@pytest.mark.parametrize(
    ("precision", "scale", "expected"),
    [
        (10, 2, r"^-?\d{1,8}(\.\d{1,2})?$"),
        (2, 2, r"^-?0(\.\d{1,2})?$"),
        (10, 0, r"^-?\d{1,10}$"),
    ],
    ids=["standard", "zero-before-decimal", "no-scale"],
)
def test_decimal_pattern_pass_regex(precision: int, scale: int, expected: str) -> None:
    """Build precise regex patterns for decimal constraints."""
    assert decimal_pattern(precision, scale) == expected


@pytest.mark.parametrize(
    ("node", "fallback", "expected"),
    [
        (TypedNode(type="string"), "text", "string"),
        (TypedNode(type="integer", declared_type=("integer", "null")), "text", "integer"),
        (TypedNode(type="unknown", constraints={"enum": (1, 2)}), "text", "integer"),
        (TypedNode(type="unknown", inferred_type="number"), "text", "number"),
        (TypedNode(type="unknown", one_of=(TypedNode(type="boolean"),)), "text", "boolean"),
        (TypedNode(type="unknown", any_of=(TypedNode(type="integer"),)), "text", "integer"),
        (TypedNode(type="unknown"), "text", "text"),
        (TypedNode(type="unknown"), "json", "json"),
    ],
    ids=[
        "explicit",
        "nullable-union",
        "enum",
        "inferred",
        "one-of",
        "any-of",
        "fallback-text",
        "fallback-custom",
    ],
)
def test_resolve_union_type_pass_precedence(node: TypedNode, fallback: str, expected: str) -> None:
    """Resolve a single type from multiple possibilities using fixed precedence."""
    assert resolve_union_type(node, fallback) == expected
