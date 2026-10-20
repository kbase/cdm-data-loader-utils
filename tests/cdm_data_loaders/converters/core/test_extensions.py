"""Closed extension contracts and custom registry isolation."""

from collections.abc import Callable
from decimal import Decimal, localcontext
from typing import Any

import pytest
from jsonschema.exceptions import SchemaError, ValidationError

from cdm_data_loaders.converters.core.extensions import (
    DEFAULT_EXTENSIONS,
    ExtensionError,
    ExtensionRegistry,
    Extensions,
    ExtensionSpec,
)
from cdm_data_loaders.converters.core.ir import TypedNode
from cdm_data_loaders.converters.core.ir_values import Value, freeze_value
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.custom_metaschema import X_XSV_CONFIG_SCHEMA


@pytest.mark.parametrize(
    ("amount", "divisor", "invalid"),
    [
        (Decimal("0.3"), 0.1, Decimal("0.31")),
        (0.3, Decimal("0.1"), 0.31),
        (0.3, 0.1, 0.31),
        (Decimal("0.3"), Decimal("0.1"), Decimal("0.31")),
        (1, 0.1, Decimal("0.31")),
        (0, Decimal("0.1"), Decimal("0.31")),
        (Decimal("-0.3"), 0.1, Decimal("-0.31")),
        (Decimal("123456789012345678901234567890.3"), 0.1, Decimal("123456789012345678901234567890.31")),
        (Decimal("1e1000"), Decimal("1e-1000"), Decimal("1e-1001")),
    ],
    ids=[
        "decimal-float",
        "float-decimal",
        "floats",
        "decimals",
        "integer",
        "zero",
        "negative",
        "many-digits",
        "exponents",
    ],
)
def test_extensions_pass_decimal_multiple_of_float(
    amount: float | Decimal, divisor: float | Decimal, invalid: float | Decimal
) -> None:
    """Fractional validation stays exact regardless of Decimal context precision."""
    registry = ExtensionRegistry((ExtensionSpec("x-amount", {"type": "number", "multipleOf": divisor}),))
    with localcontext() as context:
        context.prec = 2
        extensions = Extensions({"x-amount": amount}, registry)
        assert extensions["x-amount"] is amount
        with pytest.raises(ExtensionError, match=r"Invalid x-amount at \[\]: .*not a multiple"):
            Extensions({"x-amount": invalid}, registry)


@pytest.mark.parametrize(
    "instance_type", [int, float, Decimal], ids=["integer-instance", "float-instance", "decimal-instance"]
)
@pytest.mark.parametrize("schema_type", [int, float, Decimal], ids=["integer-schema", "float-schema", "decimal-schema"])
@pytest.mark.parametrize(
    ("keyword", "accepted", "rejected"),
    [
        ("multipleOf", "4", "3"),
        ("minimum", "2", "1"),
        ("maximum", "2", "3"),
        ("exclusiveMinimum", "3", "2"),
        ("exclusiveMaximum", "1", "2"),
    ],
    ids=["multiple-of", "minimum", "maximum", "exclusive-minimum", "exclusive-maximum"],
)
def test_extensions_pass_mixed_numeric_keywords(
    instance_type: Callable[[str], float | Decimal],
    schema_type: Callable[[str], float | Decimal],
    keyword: str,
    accepted: str,
    rejected: str,
) -> None:
    """All numeric keywords accept and reject mixed representations consistently."""
    registry = ExtensionRegistry((ExtensionSpec("x-amount", {"type": "number", keyword: schema_type("2")}),))
    amount = instance_type(accepted)
    assert Extensions({"x-amount": amount}, registry)["x-amount"] is amount
    with pytest.raises(ExtensionError) as error:
        Extensions({"x-amount": instance_type(rejected)}, registry)
    assert isinstance(error.value.__cause__, ValidationError)
    assert error.value.__cause__.validator == keyword


@pytest.mark.parametrize(
    ("amount", "boundary"),
    [(Decimal("0.1"), 0.1), (0.1, Decimal("0.1"))],
    ids=["decimal-instance", "float-instance"],
)
@pytest.mark.parametrize(
    "keyword",
    ["minimum", "maximum", "exclusiveMinimum", "exclusiveMaximum"],
    ids=["minimum", "maximum", "exclusive-minimum", "exclusive-maximum"],
)
def test_extensions_pass_exact_fractional_bounds(
    amount: float | Decimal, boundary: float | Decimal, keyword: str
) -> None:
    """Decimal and JSON float boundaries compare as equal without binary rounding."""
    registry = ExtensionRegistry((ExtensionSpec("x-amount", {"type": "number", keyword: boundary}),))
    if keyword.startswith("exclusive"):
        with pytest.raises(ExtensionError):
            Extensions({"x-amount": amount}, registry)
    else:
        assert Extensions({"x-amount": amount}, registry)["x-amount"] is amount


def test_extensions_pass_nested_numeric_refs() -> None:
    """Local refs retain exact arithmetic and nested validation error paths."""
    schema = {
        "$defs": {
            "amount": {
                "$schema": "https://json-schema.org/draft/2020-12/schema",
                "type": "number",
                "multipleOf": 0.1,
                "minimum": 0.1,
                "maximum": Decimal("0.3"),
            }
        },
        "type": "array",
        "items": {"type": "object", "properties": {"amount": {"$ref": "#/$defs/amount"}}},
    }
    registry = ExtensionRegistry((ExtensionSpec("x-amount", schema),))
    amount = Decimal("0.3")
    extensions = Extensions({"x-amount": [{"amount": amount}]}, registry)
    assert extensions["x-amount"] == ({"amount": amount},)
    assert extensions["x-amount"][0]["amount"] is amount
    with pytest.raises(ExtensionError, match=r"Invalid x-amount at \[0, 'amount'\]: .*not a multiple") as error:
        Extensions({"x-amount": [{"amount": Decimal("0.31")}]}, registry)
    assert isinstance(error.value.__cause__, ValidationError)
    assert list(error.value.__cause__.path) == [0, "amount"]


@pytest.mark.parametrize(
    "value",
    [Decimal("0.3"), 1.0, True, "text", None, ({"amount": Decimal("0.3")},)],
    ids=["decimal", "float", "boolean", "string", "null", "nested"],
)
def test_extensions_pass_boolean_schemas(value: Value) -> None:
    """True and false schemas remain unconditional for all finite payloads."""
    allowed = ExtensionRegistry((ExtensionSpec("x-value", schema=True),))
    denied = ExtensionRegistry((ExtensionSpec("x-value", schema=False),))
    assert Extensions({"x-value": value}, allowed)["x-value"] == freeze_value(value)
    with pytest.raises(ExtensionError, match="False schema"):
        Extensions({"x-value": value}, denied)


@pytest.mark.parametrize(
    ("schema", "value", "valid"),
    [
        ({"type": "number"}, Decimal("0.3"), True),
        ({"type": "number"}, True, False),
        ({"type": "integer"}, False, False),
        ({"type": "integer"}, 2.0, True),
        ({"type": "integer"}, 0.3, False),
        ({"type": "integer"}, Decimal(2), False),
        ({"type": "boolean"}, 1, False),
        ({"multipleOf": 0.1, "minimum": 1, "maximum": 2}, True, True),
        ({"multipleOf": 0.1, "minimum": 1, "maximum": 2}, "text", True),
        ({"type": "number"}, "0.3", False),
        ({"enum": [0.1]}, Decimal("0.1"), True),
        ({"const": {"amount": Decimal("0.1")}}, {"amount": 0.1}, True),
    ],
    ids=[
        "decimal-number",
        "boolean-not-number",
        "boolean-not-integer",
        "integral-float",
        "fractional-float",
        "decimal-integer-unchanged",
        "integer-not-boolean",
        "boolean-ignores-numeric",
        "string-ignores-numeric",
        "string-not-number",
        "numeric-enum",
        "nested-numeric-const",
    ],
)
def test_extensions_pass_custom_schema_types(schema: dict[str, Any], value: Value, valid: bool) -> None:
    """Numeric normalization preserves type checks and literal schema comparisons."""
    registry = ExtensionRegistry((ExtensionSpec("x-value", schema),))
    if valid:
        assert Extensions({"x-value": value}, registry)["x-value"] == freeze_value(value)
    else:
        with pytest.raises(ExtensionError):
            Extensions({"x-value": value}, registry)


@pytest.mark.parametrize(
    "value",
    [float("nan"), float("inf"), Decimal("NaN"), Decimal("Infinity")],
    ids=["float-nan", "float-infinity", "decimal-nan", "decimal-infinity"],
)
def test_extensions_fail_nonfinite_numbers(value: float | Decimal) -> None:
    """Non-finite numbers fail before numeric normalization in schemas and payloads."""
    with pytest.raises(ValueError, match="finite"):
        ExtensionSpec("x-amount", {"minimum": value})
    with pytest.raises(ValueError, match="finite"):
        Extensions({"x-amount": value}, ExtensionRegistry((ExtensionSpec("x-amount", schema=True),)))


@pytest.mark.parametrize(
    "values",
    [
        {"x-dlt": {"data_type": "decimal", "precision": 20, "scale": 2, "timezone": None}},
        {"x-iceberg": {"field_id": 1, "required": False, "initial_default": [None, Decimal("1.20")]}},
        {"x-xsv-config": {"x-delimiter": "\t", "x-null-cols": ["a"]}},
        {"x-file-glob": "*.tsv", "x-delimiter": ",", "x-dlt-prefix": "id:", "x-dlt-split": ";", "x-pii": True},
        {},
    ],
    ids=["dlt", "iceberg-default", "xsv", "discovered", "empty"],
)
def test_extensions_pass_builtin(values: dict[str, Any]) -> None:
    """Builtin namespaces validate exact internal key spellings."""
    extensions = Extensions(values)
    assert list(extensions) == list(values)
    assert len(extensions) == len(values)
    assert TypedNode(type="any", extensions=extensions).extensions is extensions


@pytest.mark.parametrize(
    "values",
    [
        {"title": "bad"},
        {"x-unknown": None},
        {"x-pii": "yes"},
        {"x-dlt": {"bogus": True}},
        {"x-iceberg": {"bogus": 1}},
        {"x-xsv-config": {}},
        {"x-delimiter": "::"},
        {"x-dlt": []},
    ],
    ids=["namespace", "unknown", "wrong-type", "dlt-key", "iceberg-key", "empty-xsv", "delimiter", "array"],
)
def test_extensions_fail_invalid(values: dict[str, Any]) -> None:
    """Unknown keys and invalid payloads fail on ordinary node construction."""
    with pytest.raises(ExtensionError):
        TypedNode(type="any", extensions=values)


def test_extension_registry_pass_custom_copy_isolation() -> None:
    """Custom arrays and null are allowed only under an explicit schema."""
    schema = {"type": "array", "items": {"type": ["string", "null"]}}
    registry = DEFAULT_EXTENSIONS.register(ExtensionSpec("x-custom", schema))
    values = {"x-custom": ["data", None]}
    extensions = Extensions(values, registry)
    schema["items"]["type"].append("number")
    values["x-custom"].append("changed")
    assert extensions["x-custom"] == ("data", None)
    with pytest.raises(ExtensionError):
        registry.validate({"x-custom": [1]})
    with pytest.raises(TypeError):
        extensions.payload["x-custom"] = ()
    custom_dlt = registry.extend_payload("x-dlt", {"x-custom": {"type": ["array", "null"]}})
    assert custom_dlt.validate({"x-dlt": {"x-custom": None}}) == {"x-dlt": {"x-custom": None}}
    with pytest.raises(ExtensionError):
        DEFAULT_EXTENSIONS.validate({"x-dlt": {"x-custom": None}})


@pytest.mark.parametrize("namespace", ["", "x-", "ordinary"], ids=["empty", "prefix-only", "no-prefix"])
def test_extension_spec_fail_namespace(namespace: str) -> None:
    """Only named x-prefixed namespaces can be registered."""
    with pytest.raises(ExtensionError):
        ExtensionSpec(namespace, schema=True)


def test_extension_registry_fail_duplicates() -> None:
    """Duplicate declarations require explicit replacement."""
    spec = ExtensionSpec("x-test", schema=True)
    registry = ExtensionRegistry((spec,))
    with pytest.raises(ExtensionError):
        ExtensionRegistry((spec, spec))
    with pytest.raises(ExtensionError):
        registry.register(spec)
    assert registry.register(spec, replace=True) == registry


def test_extension_spec_fail_invalid_schemas_and_payload_extensions() -> None:
    """Invalid schemas, wrong registry entries and payload collisions fail early."""
    with pytest.raises(ExtensionError, match="validation schemas"):
        ExtensionSpec("x-bad", 1)
    with pytest.raises(SchemaError):
        ExtensionSpec("x-bad", {"type": "bogus"})
    with pytest.raises(TypeError, match="entries"):
        ExtensionRegistry(({},))
    with pytest.raises(ExtensionError, match="Not an object"):
        DEFAULT_EXTENSIONS.extend_payload("x-pii", {"extra": True})
    with pytest.raises(ExtensionError, match="already registered"):
        DEFAULT_EXTENSIONS.extend_payload("x-dlt", {"data_type": True})
    with pytest.raises(ExtensionError, match="Unregistered"):
        DEFAULT_EXTENSIONS.extend_payload("x-missing", {"extra": True})


def test_extensions_pass_exact_xsv_schema_and_mapping_api() -> None:
    """The exact existing XSV schema is reused, and mapping methods remain callable."""
    spec = next(spec for spec in DEFAULT_EXTENSIONS.specs if spec.namespace == "x-xsv-config")
    extensions = Extensions({"x-pii": True})
    assert list(extensions.values()) == [True]
    assert list(extensions.items()) == [("x-pii", True)]
    assert spec.schema["additionalProperties"] == X_XSV_CONFIG_SCHEMA["additionalProperties"]
    assert tuple(spec.schema["properties"]) == tuple(X_XSV_CONFIG_SCHEMA["properties"])
    with pytest.raises(KeyError):
        extensions["x-missing"]
    with pytest.raises(TypeError):
        spec.schema["properties"]["x-delimiter"]["type"] = "integer"


@pytest.mark.parametrize(
    "hint",
    [
        "primary_key",
        "unique",
        "foreign_key",
        "sort",
        "cluster",
        "partition",
        "merge_key",
        "row_key",
        "root_key",
        "variant",
        "timezone",
    ],
    ids=["primary", "unique", "foreign", "sort", "cluster", "partition", "merge", "row", "root", "variant", "timezone"],
)
def test_extensions_pass_dlt_boolean_hints(hint: str) -> None:
    """Each observed dlt column hint is explicitly registered as a boolean."""
    assert Extensions({"x-dlt": {hint: True}})["x-dlt"] == {hint: True}
    with pytest.raises(ExtensionError):
        Extensions({"x-dlt": {hint: "yes"}})
