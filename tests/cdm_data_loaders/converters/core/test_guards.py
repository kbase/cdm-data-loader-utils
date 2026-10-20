"""Tests for shared JSON Schema guards and caller-selected errors."""

from typing import Any

import pytest

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.guards import (
    reject_unresolved_references,
    require_object_root,
    require_schema_keyword,
)
from cdm_data_loaders.converters.jsonschema_to_dlt import JSONSchemaToDlt, JSONSchemaToDltError
from cdm_data_loaders.converters.jsonschema_to_pyspark import JSONSchemaToPySpark, JSONSchemaToPySparkError


@pytest.mark.parametrize("error", [JSONSchemaToDltError, JSONSchemaToPySparkError], ids=["dlt", "pyspark"])
@pytest.mark.parametrize("schema", [{}, {"$schema": ""}, {"$schema": None}], ids=["missing", "empty", "null"])
def test_require_schema_keyword_fail_missing_dialect(schema: dict[str, Any], error: type[ConversionError]) -> None:
    """Reject absent dialect declarations with the exact caller-selected exception."""
    with pytest.raises(error, match="missing a '\\$schema'") as caught:
        require_schema_keyword(schema, error)
    assert type(caught.value) is error


@pytest.mark.parametrize(
    "dialect", ["https://json-schema.org/draft/2020-12/schema", "custom"], ids=["standard", "custom"]
)
def test_require_schema_keyword_pass_declared_dialect(dialect: str) -> None:
    """Require a declaration without validating the dialect URI."""
    assert require_schema_keyword({"$schema": dialect}, ConversionError) is None


@pytest.mark.parametrize("schema", [{}, {"type": None}, {"type": "object"}], ids=["implicit", "unspecified", "object"])
def test_require_object_root_pass_object_or_unspecified(schema: dict[str, Any]) -> None:
    """Allow an object or an unspecified root type."""
    assert require_object_root(schema, ConversionError, target="a table") is None


@pytest.mark.parametrize(
    "root_type",
    ["string", "array", "number", "integer", "boolean", "null", ["object", "null"], ""],
    ids=["string", "array", "number", "integer", "boolean", "null", "union", "empty"],
)
@pytest.mark.parametrize(
    ("error", "target"),
    [(JSONSchemaToDltError, "a dlt table"), (JSONSchemaToPySparkError, "a StructType")],
    ids=["dlt", "pyspark"],
)
def test_require_object_root_fail_non_object(
    root_type: str | list[str],
    error: type[ConversionError],
    target: str,
) -> None:
    """Include the target and declared root type in direction-specific errors."""
    with pytest.raises(error) as caught:
        require_object_root({"type": root_type}, error, target=target)
    assert type(caught.value) is error
    assert str(caught.value) == f"Root schema must be of type 'object' to map to {target}, got: {root_type!r}"


@pytest.mark.parametrize("schema", [{}, {"type": "string"}], ids=["empty", "scalar"])
def test_reject_unresolved_references_pass_resolved_schema(schema: dict[str, Any]) -> None:
    """Accept fragments without unresolved references or allOf."""
    assert reject_unresolved_references(schema, ConversionError) is None


@pytest.mark.parametrize("error", [JSONSchemaToDltError, JSONSchemaToPySparkError], ids=["dlt", "pyspark"])
@pytest.mark.parametrize(
    ("schema", "message"),
    [
        ({"$ref": "#/$defs/value"}, "unresolved \\$ref"),
        ({"$ref": "https://example.org/schema.json"}, "unresolved \\$ref"),
        ({"$ref": ""}, "unresolved \\$ref"),
        ({"allOf": []}, "unmerged 'allOf'"),
        ({"$ref": "#/$defs/value", "allOf": []}, "unresolved \\$ref"),
    ],
    ids=["local-ref", "external-ref", "empty-ref", "all-of", "ref-before-all-of"],
)
def test_reject_unresolved_references_fail_unresolved(
    schema: dict[str, Any],
    message: str,
    error: type[ConversionError],
) -> None:
    """Reject unresolved keywords by presence, checking references before allOf."""
    with pytest.raises(error, match=message) as caught:
        reject_unresolved_references(schema, error)
    assert type(caught.value) is error


@pytest.mark.parametrize(
    ("converter", "error", "module"),
    [
        (JSONSchemaToDlt(schema_name="test"), JSONSchemaToDltError, "jsonschema_to_dlt"),
        (JSONSchemaToPySpark(), JSONSchemaToPySparkError, "jsonschema_to_pyspark"),
    ],
    ids=["dlt", "pyspark"],
)
@pytest.mark.parametrize(
    "keyword", ["$schema", "$ref", "allOf"], ids=["missing-dialect", "unresolved-ref", "unmerged-all-of"]
)
def test_convert_fail_direction_specific_diagnostics(
    converter: JSONSchemaToDlt | JSONSchemaToPySpark, error: type[ConversionError], module: str, keyword: str
) -> None:
    """Preserve converter names, dereferencing paths and complete guard messages."""
    name = type(converter).__name__
    schema: dict[str, Any] = {"$schema": "https://json-schema.org/draft/2020-12/schema", "type": "object"}
    if keyword == "$schema":
        del schema["$schema"]
        expected = (
            f"Input JSON Schema is missing a '$schema' keyword. {name} requires schemas "
            "to explicitly declare their dialect via '$schema'; it will not assume a default."
        )
    elif keyword == "$ref":
        schema["$ref"] = "#/$defs/value"
        expected = (
            f"Encountered an unresolved $ref '#/$defs/value'. {name} requires a fully "
            "dereferenced schema -- this includes references to external JSON Schema documents. Use "
            f"`{module}.dereferencing.dereference_schema()` to resolve all $refs before calling `convert()`."
        )
    else:
        schema["allOf"] = []
        expected = (
            f"Encountered an unmerged 'allOf'. {name} requires a fully dereferenced schema "
            f"with 'allOf' already merged. Use `{module}.dereferencing.dereference_schema()` before calling `convert()`."
        )
    with pytest.raises(error) as caught:
        converter.convert(schema)
    assert str(caught.value) == expected
