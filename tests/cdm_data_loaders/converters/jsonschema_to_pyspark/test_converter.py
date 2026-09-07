"""Unit tests for jsonschema_to_pyspark.converter."""

# ruff: noqa: SLF001
import json
import logging
from pathlib import Path
from typing import Any

import pytest
import yaml
from jsonschema import Draft7Validator, Draft202012Validator
from pydantic import ValidationError
from pyspark.sql.types import (
    ArrayType,
    BooleanType,
    DateType,
    DecimalType,
    DoubleType,
    IntegerType,
    LongType,
    MapType,
    NullType,
    StringType,
    StructField,
    StructType,
    TimestampType,
)

from cdm_data_loaders.converters.jsonschema_to_pyspark.converter import (
    ConversionContext,
    InvalidJSONSchemaError,
    JSONSchemaToPySpark,
    JSONSchemaToPySparkError,
    _decimal_places,
    _infer_implicit_type,
    _infer_type_from_enum,
    _metadata_keys_for,
    get_known_jsonschema_keywords,
)
from cdm_data_loaders.converters.jsonschema_to_pyspark.dereferencer import dereference_schema
from tests.cdm_data_loaders.converters.jsonschema_to_pyspark.conftest import base_object_schema


def test_json_schema_to_pyspark_fail_instance_is_frozen(converter: JSONSchemaToPySpark) -> None:
    """Assigning to a JSONSchemaToPySpark instance attribute raises ValidationError (frozen model)."""
    with pytest.raises(ValidationError, match="Instance is frozen"):
        converter.treat_unknown_as_string = False


"""convert"""


def test_convert_fail_missing_schema_keyword(converter: JSONSchemaToPySpark) -> None:
    """convert() rejects a schema with no top-level '$schema' keyword."""
    schema = {"type": "object", "properties": {}}
    with pytest.raises(InvalidJSONSchemaError, match="missing a '\\$schema'"):
        converter.convert(schema)


@pytest.mark.parametrize(
    "root_type",
    ["string", "array", "integer", "number", "boolean", "null"],
)
def test_convert_fail_non_object_root_type(converter: JSONSchemaToPySpark, root_type: str) -> None:
    """convert() rejects any root schema whose declared 'type' isn't 'object'."""
    schema = base_object_schema(type=root_type)
    with pytest.raises(JSONSchemaToPySparkError, match="must be of type 'object'"):
        converter.convert(schema)


def test_convert_fail_unresolved_ref_at_root(converter: JSONSchemaToPySpark) -> None:
    """convert() refuses a root schema containing an un-dereferenced $ref."""
    schema = base_object_schema(**{"$ref": "#/$defs/Foo"})
    del schema["type"]
    with pytest.raises(JSONSchemaToPySparkError, match="unresolved \\$ref"):
        converter.convert(schema)


def test_convert_fail_root_resolves_to_maptype(converter: JSONSchemaToPySpark) -> None:
    """convert() rejects a root schema that has no 'properties' and resolves to a MapType."""
    schema = base_object_schema(patternProperties={"^x-": {"type": "string"}})
    del schema["properties"]
    with pytest.raises(JSONSchemaToPySparkError, match="did not resolve to a StructType"):
        converter.convert(schema)


def test_convert_pass_simple_object(converter: JSONSchemaToPySpark) -> None:
    """convert() produces a StructType with correct nullability and empty metadata for a flat object schema."""
    schema = base_object_schema(
        properties={
            "name": {"type": "string"},
            "age": {"type": "integer"},
        },
        required=["name"],
    )
    result = converter.convert(schema)
    assert result == StructType(
        [
            StructField("name", StringType(), nullable=False, metadata={}),
            StructField("age", LongType(), nullable=True, metadata={}),
        ]
    )


def test_convert_pass_root_schema_without_type_keyword_infers_object(converter: JSONSchemaToPySpark) -> None:
    """convert() succeeds for a root schema that omits 'type' entirely but declares 'properties' (valid JSON Schema)."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "properties": {"a": {"type": "string"}},
    }
    result = converter.convert(schema)
    assert result == StructType([StructField("a", StringType(), nullable=True, metadata={})])


def test_convert_pass_nullable_object_root_type_accepted(converter: JSONSchemaToPySpark) -> None:
    """convert() accepts a root schema whose 'type' is a list consisting solely of 'object'/'null'-style is rejected as non-'object'."""
    schema = base_object_schema(type=["object", "null"])
    with pytest.raises(JSONSchemaToPySparkError, match="must be of type 'object'"):
        converter.convert(schema)


# ConversionContext.metadata_keys


def test_conversion_context_pass_metadata_keys_delegates_to_metadata_keys_for() -> None:
    """ConversionContext.metadata_keys returns the same value as _metadata_keys_for() for its validator_cls."""
    context = ConversionContext(validator_cls=Draft202012Validator)
    assert context.metadata_keys == _metadata_keys_for(Draft202012Validator)


"""convert_from_string / convert_from_file"""


def test_convert_from_string_pass_valid_json(converter: JSONSchemaToPySpark) -> None:
    """convert_from_string() parses JSON text and converts it identically to convert()."""
    schema = base_object_schema(properties={"id": {"type": "integer"}})
    result = converter.convert_from_string(json.dumps(schema))
    assert result == converter.convert(schema)


def test_convert_from_string_fail_invalid_json(converter: JSONSchemaToPySpark) -> None:
    """convert_from_string() propagates a JSONDecodeError for malformed JSON text."""
    with pytest.raises(json.JSONDecodeError):
        converter.convert_from_string("{not valid json")


@pytest.mark.parametrize(
    ("suffix", "dumper"),
    [(".json", json.dumps), (".yaml", yaml.safe_dump)],
)
def test_convert_from_file_pass_supported_extensions(
    converter: JSONSchemaToPySpark, tmp_path: Path, suffix: str, dumper: Any
) -> None:
    """convert_from_file() loads both .json and non-.json (YAML) files correctly."""
    schema = base_object_schema(properties={"id": {"type": "integer"}})
    path = tmp_path / f"schema{suffix}"
    path.write_text(dumper(schema))
    result = converter.convert_from_file(str(path))
    assert result == converter.convert(schema)


def test_convert_from_file_fail_missing_file(converter: JSONSchemaToPySpark, tmp_path: Path) -> None:
    """convert_from_file() raises FileNotFoundError for a nonexistent path."""
    missing = tmp_path / "does_not_exist.json"
    with pytest.raises(FileNotFoundError):
        converter.convert_from_file(str(missing))


"""_reject_unresolved_references"""


@pytest.mark.parametrize(
    ("schema", "match"),
    [
        ({"$ref": "#/$defs/Foo"}, "unresolved \\$ref"),
        ({"allOf": [{"type": "string"}]}, "unmerged 'allOf'"),
    ],
)
def test_reject_unresolved_references_fail_ref_or_all_of_present(schema: dict[str, Any], match: str) -> None:
    """_reject_unresolved_references() raises for schemas still containing $ref or allOf."""
    with pytest.raises(JSONSchemaToPySparkError, match=match):
        JSONSchemaToPySpark._reject_unresolved_references(schema)


def test_reject_unresolved_references_pass_clean_schema() -> None:
    """_reject_unresolved_references() is a no-op for a schema with no $ref/allOf."""
    assert JSONSchemaToPySpark._reject_unresolved_references({"type": "string"}) is None


def test_reject_unresolved_references_fail_ref_checked_before_all_of() -> None:
    """_reject_unresolved_references() reports the $ref error first when both $ref and allOf are present."""
    schema = {"$ref": "#/$defs/Foo", "allOf": [{"type": "string"}]}
    with pytest.raises(JSONSchemaToPySparkError, match="unresolved \\$ref"):
        JSONSchemaToPySpark._reject_unresolved_references(schema)


"""_convert_boolean_schema"""


@pytest.mark.parametrize("bool_schema", [True, False])
def test_convert_boolean_schema_pass_falls_back_to_string_by_default(
    converter: JSONSchemaToPySpark, bool_schema: bool
) -> None:
    """_convert_boolean_schema() maps both `true` and `false` schemas to StringType when treat_unknown_as_string=True."""
    assert converter._convert_boolean_schema(bool_schema) == StringType()


@pytest.mark.parametrize("bool_schema", [True, False])
def test_convert_boolean_schema_fail_raises_when_disallowed(
    strict_converter: JSONSchemaToPySpark, bool_schema: bool
) -> None:
    """_convert_boolean_schema() raises for both `true` and `false` schemas when treat_unknown_as_string=False."""
    with pytest.raises(JSONSchemaToPySparkError, match="no PySpark equivalent"):
        strict_converter._convert_boolean_schema(bool_schema)


"""_build_metadata"""


def test_build_metadata_pass_description_becomes_comment_only(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_build_metadata() surfaces 'description' only as top-level 'comment', not inside 'jsonschema'."""
    metadata = converter._build_metadata({"type": "string", "description": "A name"}, ctx)
    assert metadata == {"comment": "A name"}


def test_build_metadata_pass_title_goes_into_jsonschema_dict(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_build_metadata() copies 'title' into the nested 'jsonschema' metadata dict."""
    metadata = converter._build_metadata({"type": "string", "title": "Name"}, ctx)
    assert metadata == {"jsonschema": {"title": "Name"}}


def test_build_metadata_pass_no_matching_keywords_returns_empty_dict(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_build_metadata() returns an entirely empty dict when there's no title/description/extra data."""
    metadata = converter._build_metadata({"type": "string", "pattern": "^a"}, ctx)
    assert metadata == {}


def test_build_metadata_pass_excludes_structural_keywords_even_as_allowed_extra() -> None:
    """_build_metadata() never copies structural/compositional keywords, even if listed in extra_metadata_keywords."""
    # 'properties'/'type' can never actually reach ctx.allowed_extra_metadata_keywords in practice
    # (they're excluded from ctx.metadata_keys), but this guards the invariant directly.
    converter = JSONSchemaToPySpark(extra_metadata_keywords={"properties", "type"})  # pyright: ignore[reportArgumentType]
    ctx = converter._build_context(Draft202012Validator)
    schema = {"type": "object", "properties": {"a": {"type": "string"}}}
    metadata = converter._build_metadata(schema, ctx)
    assert metadata == {}


@pytest.mark.parametrize(
    ("extra_keywords", "schema_extra_key", "should_appear"),
    [
        ({"pattern"}, "pattern", True),  # known standard keyword -> allowed
        ({"x-pii"}, "x-pii", True),  # vendor 'x-' prefix -> allowed
        ({"totally-made-up"}, "totally-made-up", False),  # neither -> filtered out
    ],
)
def test_build_metadata_pass_extra_metadata_keywords_filtering(
    extra_keywords: set[str],
    schema_extra_key: str,
    should_appear: bool,
) -> None:
    """_build_metadata() only honours ctx.allowed_extra_metadata_keywords (standard keywords or 'x-' prefixed)."""
    converter = JSONSchemaToPySpark(extra_metadata_keywords=extra_keywords)  # pyright: ignore[reportArgumentType]
    ctx = converter._build_context(Draft202012Validator)
    schema = {"type": "string", schema_extra_key: "some-value"}
    metadata = converter._build_metadata(schema, ctx)
    if should_appear:
        assert metadata["jsonschema"][schema_extra_key] == "some-value"
    else:
        assert metadata == {}


def test_build_metadata_pass_title_and_description_together(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_build_metadata() combines 'title' (nested) and 'description' (top-level comment) when both are present."""
    metadata = converter._build_metadata({"type": "string", "title": "Name", "description": "desc"}, ctx)
    assert metadata == {"jsonschema": {"title": "Name"}, "comment": "desc"}


"""_convert_type / _convert_property"""


@pytest.mark.parametrize("bool_schema", [True, False])
def test_convert_type_pass_boolean_schema_is_stringtype(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, bool_schema: bool
) -> None:
    """_convert_type() maps JSON Schema boolean schemas (true/false) to StringType by default."""
    assert converter._convert_type(bool_schema, ctx) == StringType()


@pytest.mark.parametrize("bool_schema", [True, False])
def test_convert_type_fail_boolean_schema_when_disallowed(
    strict_converter: JSONSchemaToPySpark, strict_ctx: ConversionContext, bool_schema: bool
) -> None:
    """_convert_type() raises for boolean schemas when treat_unknown_as_string=False."""
    with pytest.raises(JSONSchemaToPySparkError, match="no PySpark equivalent"):
        strict_converter._convert_type(bool_schema, strict_ctx)


def test_convert_type_fail_unresolved_ref(converter: JSONSchemaToPySpark, ctx: ConversionContext) -> None:
    """_convert_type() raises when handed a schema still containing $ref."""
    with pytest.raises(JSONSchemaToPySparkError, match="unresolved \\$ref"):
        converter._convert_type({"$ref": "#/$defs/Foo"}, ctx)


@pytest.mark.parametrize("bool_schema", [True, False])
def test_convert_property_pass_boolean_schema_returns_string_and_empty_metadata(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, bool_schema: bool
) -> None:
    """_convert_property() maps boolean schemas to (StringType, {}) with no metadata by default."""
    data_type, metadata = converter._convert_property(bool_schema, ctx)
    assert data_type == StringType()
    assert metadata == {}


@pytest.mark.parametrize("bool_schema", [True, False])
def test_convert_property_fail_boolean_schema_when_disallowed(
    strict_converter: JSONSchemaToPySpark, strict_ctx: ConversionContext, bool_schema: bool
) -> None:
    """_convert_property() raises for boolean schemas when treat_unknown_as_string=False."""
    with pytest.raises(JSONSchemaToPySparkError, match="no PySpark equivalent"):
        strict_converter._convert_property(bool_schema, strict_ctx)


def test_convert_property_pass_returns_type_and_metadata(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_property() returns both the converted DataType and built metadata."""
    data_type, metadata = converter._convert_property({"type": "string", "description": "desc"}, ctx)
    assert data_type == StringType()
    assert metadata == {"comment": "desc"}


"""_dispatch_type"""


@pytest.mark.parametrize(
    ("schema", "expected"),
    [
        ({"type": "string"}, StringType()),
        ({"type": "boolean"}, BooleanType()),
        ({"type": "null"}, NullType()),
        ({"type": "integer"}, LongType()),
        ({"type": "number"}, DoubleType()),
        ({"type": "object"}, StructType([])),
        ({"type": "array", "items": {"type": "string"}}, ArrayType(StringType(), containsNull=True)),
    ],
)
def test_dispatch_type_pass_scalar_and_container_types(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, schema: dict[str, Any], expected: type
) -> None:
    """_dispatch_type() routes each JSON Schema 'type' value to the correct PySpark DataType."""
    assert converter._dispatch_type(schema, ctx) == expected


def test_dispatch_type_pass_multi_type_union_single_non_null_type(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() collapses a ['string', 'null'] union to StringType."""
    assert converter._dispatch_type({"type": ["string", "null"]}, ctx) == StringType()


def test_dispatch_type_pass_multi_type_union_all_null(converter: JSONSchemaToPySpark, ctx: ConversionContext) -> None:
    """_dispatch_type() maps a type list containing only 'null' to NullType."""
    assert converter._dispatch_type({"type": ["null"]}, ctx) == NullType()


def test_dispatch_type_pass_multi_type_union_collapses_to_string_by_default(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() collapses a multi-non-null-type union to StringType when treat_unknown_as_string=True."""
    assert converter._dispatch_type({"type": ["string", "integer"]}, ctx) == StringType()


def test_dispatch_type_fail_multi_type_union_when_disallowed(
    strict_converter: JSONSchemaToPySpark, strict_ctx: ConversionContext
) -> None:
    """_dispatch_type() raises on a multi-non-null-type union when treat_unknown_as_string=False."""
    with pytest.raises(JSONSchemaToPySparkError, match="multi-type union"):
        strict_converter._dispatch_type({"type": ["string", "integer"]}, strict_ctx)


@pytest.mark.parametrize(
    ("enum_values", "expected"),
    [
        ([True, False], BooleanType()),
        ([1, 2, 3], LongType()),
        ([1, 2.5], DoubleType()),
        (["a", "b"], StringType()),
        ([], StringType()),
    ],
)
def test_dispatch_type_pass_typeless_enum_infers_type(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, enum_values: list[Any], expected: type
) -> None:
    """_dispatch_type() infers a DataType from 'enum' values when no 'type' keyword is present."""
    assert converter._dispatch_type({"enum": enum_values}, ctx) == expected


@pytest.mark.parametrize("combiner", ["oneOf", "anyOf"])
def test_dispatch_type_pass_combiner_uses_first_branch(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, combiner: str
) -> None:
    """_dispatch_type() approximates 'oneOf'/'anyOf' by converting only the first branch."""
    schema = {combiner: [{"type": "integer"}, {"type": "string"}]}
    assert converter._dispatch_type(schema, ctx) == LongType()


def test_dispatch_type_pass_unknown_construct_falls_back_to_string(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() falls back to StringType for unrecognized schemas when treat_unknown_as_string=True."""
    assert converter._dispatch_type({"not": {"type": "string"}}, ctx) == StringType()


def test_dispatch_type_fail_unknown_construct_when_disallowed(
    strict_converter: JSONSchemaToPySpark, strict_ctx: ConversionContext
) -> None:
    """_dispatch_type() raises for unrecognized schemas when treat_unknown_as_string=False."""
    with pytest.raises(JSONSchemaToPySparkError, match="Unsupported/unknown schema type"):
        strict_converter._dispatch_type({"not": {"type": "string"}}, strict_ctx)


@pytest.mark.parametrize("ignored_keyword", ["not", "if", "then", "else"])
def test_dispatch_type_pass_ignored_conditional_keyword_logs_warning(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, ignored_keyword: str, caplog: pytest.LogCaptureFixture
) -> None:
    """_dispatch_type() logs a warning for each unsupported conditional keyword ('not'/'if'/'then'/'else') it ignores."""
    schema = {ignored_keyword: {"type": "string"}}
    with caplog.at_level(logging.WARNING):
        result = converter._dispatch_type(schema, ctx)
    assert result == StringType()
    assert ignored_keyword in caplog.text


@pytest.mark.parametrize("combiner", ["oneOf", "anyOf"])
def test_dispatch_type_pass_empty_combiner_list_falls_through_to_unknown_handling(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, combiner: str
) -> None:
    """_dispatch_type() treats an empty 'oneOf'/'anyOf' list as absent, falling through to the unknown-construct fallback."""
    assert converter._dispatch_type({combiner: []}, ctx) == StringType()


@pytest.mark.parametrize("combiner", ["oneOf", "anyOf"])
def test_dispatch_type_fail_empty_combiner_list_when_disallowed(
    strict_converter: JSONSchemaToPySpark, strict_ctx: ConversionContext, combiner: str
) -> None:
    """_dispatch_type() raises for an empty 'oneOf'/'anyOf' list when treat_unknown_as_string=False."""
    with pytest.raises(JSONSchemaToPySparkError, match="Unsupported/unknown schema type"):
        strict_converter._dispatch_type({combiner: []}, strict_ctx)


"""_dispatch_type"""


def test_dispatch_type_pass_if_then_else_silently_ignored_when_type_present(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, caplog: pytest.LogCaptureFixture
) -> None:
    """_dispatch_type() drops 'if'/'then'/'else' with no warning when a recognized 'type' is also present."""
    schema = {
        "type": "object",
        "if": {"properties": {"country": {"const": "US"}}},
        "then": {"required": ["zip"]},
        "else": {"required": ["postal_code"]},
        "properties": {"country": {"type": "string"}},
    }
    with caplog.at_level(logging.WARNING):
        result = converter._dispatch_type(schema, ctx)
    assert isinstance(result, StructType)
    assert "if" not in caplog.text


def test_dispatch_type_pass_not_silently_ignored_when_type_present(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, caplog: pytest.LogCaptureFixture
) -> None:
    """_dispatch_type() drops 'not' with no warning when a recognized 'type' is also present."""
    schema = {"type": "string", "not": {"enum": ["forbidden"]}}
    with caplog.at_level(logging.WARNING):
        result = converter._dispatch_type(schema, ctx)
    assert result == StringType()
    assert "not" not in caplog.text.lower() or "unsupported conditional keyword" not in caplog.text


@pytest.mark.parametrize("combiner", ["oneOf", "anyOf"])
def test_dispatch_type_pass_combiner_silently_ignored_when_type_present(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, combiner: str, caplog: pytest.LogCaptureFixture
) -> None:
    """_dispatch_type() drops 'oneOf'/'anyOf' entirely if `type` value is present.

    Base 'type' takes precedence over the union approximation.
    """
    schema = {
        "type": "object",
        "properties": {"a": {"type": "string"}},
        combiner: [{"type": "integer"}],
    }
    with caplog.at_level(logging.WARNING):
        result = converter._dispatch_type(schema, ctx)
    assert isinstance(result, StructType)
    assert result.fieldNames() == ["a"]
    assert "Approximating" not in caplog.text


# _dispatch type with _infer_implicit_type
def test_dispatch_type_pass_infers_object_when_type_omitted(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() infers 'object' and produces a StructType for a schema with 'properties' but no 'type'."""
    schema = {"properties": {"a": {"type": "string"}}, "required": ["a"]}
    result = converter._dispatch_type(schema, ctx)
    assert result == StructType([StructField("a", StringType(), nullable=False, metadata={})])


def test_dispatch_type_pass_infers_array_when_type_omitted(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() infers 'array' and produces an ArrayType for a schema with 'items' but no 'type'."""
    result = converter._dispatch_type({"items": {"type": "integer"}}, ctx)
    assert result == ArrayType(LongType(), containsNull=True)


def test_dispatch_type_pass_infers_string_when_type_omitted(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() infers 'string' and honours 'format' for a schema with string keywords but no 'type'."""
    result = converter._dispatch_type({"format": "date"}, ctx)
    assert result == DateType()


def test_dispatch_type_pass_infers_number_when_type_omitted(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() infers 'number' and produces a DecimalType for a schema with 'multipleOf' but no 'type'."""
    result = converter._dispatch_type({"multipleOf": 0.01}, ctx)
    assert result == DecimalType(38, 2)


def test_dispatch_type_pass_enum_still_takes_priority_over_implicit_inference(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() still infers type from 'enum' first, even if other type-specific keywords are also present."""
    schema = {"enum": [1, 2, 3], "minimum": 0}  # would otherwise imply 'number' -> DoubleType
    assert converter._dispatch_type(schema, ctx) == LongType()


def test_dispatch_type_pass_no_inferrable_keywords_falls_back_to_string(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_dispatch_type() still falls back to StringType for a genuinely unconstrained, type-less schema."""
    assert converter._dispatch_type({"description": "anything goes"}, ctx) == StringType()


def test_dispatch_type_fail_no_inferrable_keywords_when_disallowed(
    strict_converter: JSONSchemaToPySpark, strict_ctx: ConversionContext
) -> None:
    """_dispatch_type() still raises for a genuinely unconstrained, type-less schema when treat_unknown_as_string=False."""
    with pytest.raises(JSONSchemaToPySparkError, match="Unsupported/unknown schema type"):
        strict_converter._dispatch_type({"description": "anything goes"}, strict_ctx)


"""_convert_object"""


def test_convert_object_pass_properties_and_required_control_nullability(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_object() marks properties listed in 'required' as non-nullable, others nullable."""
    schema = {
        "type": "object",
        "properties": {"a": {"type": "string"}, "b": {"type": "string"}},
        "required": ["a"],
    }
    result = converter._convert_object(schema, ctx)
    assert result == StructType(
        [StructField("a", StringType(), nullable=False), StructField("b", StringType(), nullable=True)]
    )


def test_convert_object_pass_pattern_properties_maps_to_maptype(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_object() maps an object with only 'patternProperties' to a MapType."""
    schema = {"type": "object", "patternProperties": {"^x-": {"type": "integer"}}}
    result = converter._convert_object(schema, ctx)
    assert result == MapType(StringType(), LongType(), valueContainsNull=True)


def test_convert_object_pass_additional_properties_schema_maps_to_maptype(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_object() maps an object with schema-valued 'additionalProperties' to a MapType."""
    schema = {"type": "object", "additionalProperties": {"type": "string"}}
    result = converter._convert_object(schema, ctx)
    assert result == MapType(StringType(), StringType(), valueContainsNull=True)


def test_convert_object_pass_no_properties_returns_empty_structtype(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_object() falls back to an empty StructType for a fully open object schema."""
    result = converter._convert_object({"type": "object"}, ctx)
    assert result == StructType([])


def test_convert_object_pass_properties_take_precedence_over_pattern_properties(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_object() prefers fixed 'properties' over 'patternProperties' when both are present."""
    schema = {
        "type": "object",
        "properties": {"a": {"type": "string"}},
        "patternProperties": {"^x-": {"type": "integer"}},
    }
    result = converter._convert_object(schema, ctx)
    assert isinstance(result, StructType)
    assert result.fieldNames() == ["a"]


@pytest.mark.parametrize("additional_properties", [True, False])
def test_convert_object_pass_boolean_additional_properties_falls_back_to_empty_structtype(
    converter: JSONSchemaToPySpark, ctx: ConversionContext, additional_properties: bool
) -> None:
    """_convert_object() treats boolean (not schema-valued) 'additionalProperties' as no dynamic-value schema."""
    schema = {"type": "object", "additionalProperties": additional_properties}
    result = converter._convert_object(schema, ctx)
    assert result == StructType([])


"""_convert_array"""


def test_convert_array_pass_simple_items(converter: JSONSchemaToPySpark, ctx: ConversionContext) -> None:
    """_convert_array() converts a homogeneous 'items' schema to the matching ArrayType."""
    result = converter._convert_array({"type": "array", "items": {"type": "string"}}, ctx)
    assert result == ArrayType(StringType(), containsNull=True)


def test_convert_array_pass_prefix_items_uses_first_slot(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_array() uses only the first 'prefixItems' schema as the element type."""
    schema = {"type": "array", "prefixItems": [{"type": "integer"}, {"type": "string"}]}
    result = converter._convert_array(schema, ctx)
    assert result == ArrayType(LongType(), containsNull=True)


def test_convert_array_pass_tuple_style_items_uses_first_element(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_array() uses only the first schema of a Draft-07 tuple-style 'items' list."""
    schema = {"type": "array", "items": [{"type": "boolean"}, {"type": "string"}]}
    result = converter._convert_array(schema, ctx)
    assert result == ArrayType(BooleanType(), containsNull=True)


def test_convert_array_pass_missing_items_defaults_to_string_element(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_array() defaults the element type to StringType when 'items' is absent."""
    result = converter._convert_array({"type": "array"}, ctx)
    assert result == ArrayType(StringType(), containsNull=True)


def test_convert_array_pass_empty_prefix_items_defaults_to_string_element(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_array() defaults to StringType element when 'prefixItems' is an empty list."""
    result = converter._convert_array({"type": "array", "prefixItems": []}, ctx)
    assert result == ArrayType(StringType(), containsNull=True)


def test_convert_array_pass_empty_tuple_style_items_defaults_to_string_element(
    converter: JSONSchemaToPySpark, ctx: ConversionContext
) -> None:
    """_convert_array() defaults to StringType element when tuple-style 'items' is an empty list."""
    result = converter._convert_array({"type": "array", "items": []}, ctx)
    assert result == ArrayType(StringType(), containsNull=True)


"""_convert_string"""


@pytest.mark.parametrize(
    ("fmt", "expected"),
    [
        ("date", DateType()),
        ("date-time", TimestampType()),
        ("uuid", StringType()),
        ("time", StringType()),
    ],
)
def test_convert_string_pass_known_format_maps_to_expected_type(
    converter: JSONSchemaToPySpark, fmt: str, expected: type
) -> None:
    """_convert_string() maps known 'format' values per DEFAULT_FORMAT_MAP."""
    assert converter._convert_string({"type": "string", "format": fmt}) == expected


def test_convert_string_pass_unknown_format_defaults_to_string(converter: JSONSchemaToPySpark) -> None:
    """_convert_string() falls back to StringType for a 'format' value not in format_map."""
    assert converter._convert_string({"type": "string", "format": "not-a-real-format"}) == StringType()


def test_convert_string_pass_no_format_returns_string(converter: JSONSchemaToPySpark) -> None:
    """_convert_string() returns StringType when no 'format' keyword is present."""
    assert converter._convert_string({"type": "string"}) == StringType()


def test_convert_string_pass_custom_format_map_override_takes_precedence() -> None:
    """A user-supplied format_map entry overrides the corresponding DEFAULT_FORMAT_MAP entry."""
    converter = JSONSchemaToPySpark(format_map={"date-time": StringType()})  # pyright: ignore[reportArgumentType]
    assert converter._convert_string({"type": "string", "format": "date-time"}) == StringType()
    # untouched entries still come from the default map
    assert converter._convert_string({"type": "string", "format": "date"}) == DateType()


"""_convert_integer"""


@pytest.mark.parametrize(
    ("schema", "expected"),
    [
        ({"type": "integer", "minimum": 0, "maximum": 100}, IntegerType()),
        ({"type": "integer", "minimum": -2_147_483_648, "maximum": 2_147_483_647}, IntegerType()),
        ({"type": "integer", "minimum": 0, "maximum": 9_999_999_999}, LongType()),
        ({"type": "integer", "exclusiveMinimum": 0, "exclusiveMaximum": 100}, IntegerType()),
        ({"type": "integer"}, LongType()),
        ({"type": "integer", "minimum": 0}, LongType()),
    ],
)
def test_convert_integer_pass_bounds_determine_int_or_long(schema: dict[str, Any], expected: type) -> None:
    """_convert_integer() picks IntegerType when bounds fit in int32, else LongType."""
    assert JSONSchemaToPySpark._convert_integer(schema) == expected


"""_convert_number"""


@pytest.mark.parametrize(
    ("multiple_of", "expected_scale"),
    [(1, 0), (0.1, 1), (0.01, 2), (0.001, 3)],
)
def test_convert_number_pass_multiple_of_computes_decimal_scale(multiple_of: float, expected_scale: int) -> None:
    """_convert_number() derives DecimalType scale from the 'multipleOf' keyword."""
    result = JSONSchemaToPySpark._convert_number({"type": "number", "multipleOf": multiple_of})
    assert result == DecimalType(precision=38, scale=expected_scale)


def test_convert_number_pass_without_multiple_of_returns_double() -> None:
    """_convert_number() returns DoubleType when 'multipleOf' is absent."""
    assert JSONSchemaToPySpark._convert_number({"type": "number"}) == DoubleType()


"""_decimal_places"""


@pytest.mark.parametrize(
    ("value", "expected"),
    [(1, 0), (1.0, 1), (0.1, 1), (0.25, 2), (0.001, 3), (100, 0), (100.000001, 6)],
)
def test_decimal_places_pass_various_magnitudes(value: float, expected: int) -> None:
    """_decimal_places() returns the number of fractional digits needed to represent `value`."""
    assert _decimal_places(value) == expected


"""_infer_type_from_enum"""


@pytest.mark.parametrize(
    ("values", "expected"),
    [
        ([], StringType()),
        ([True, False], BooleanType()),
        ([1, 2, 3], LongType()),
        ([1, 2.5, 3], DoubleType()),
        (["a", "b", "c"], StringType()),
        ([1, "a"], StringType()),
        ([None], StringType()),
    ],
)
def test_infer_type_from_enum_pass_various_value_sets(values: list[Any], expected: list[Any]) -> None:
    """_infer_type_from_enum() infers a PySpark type consistent with a homogeneous enum value list."""
    assert _infer_type_from_enum(values) == expected


"""_infer_implicit_type"""


@pytest.mark.parametrize(
    ("schema", "expected"),
    [
        ({"properties": {"a": {"type": "string"}}}, "object"),
        ({"patternProperties": {"^x-": {}}}, "object"),
        ({"additionalProperties": {"type": "string"}}, "object"),
        ({"required": ["a"]}, "object"),
        ({"items": {"type": "string"}}, "array"),
        ({"prefixItems": [{"type": "string"}]}, "array"),
        ({"uniqueItems": True}, "array"),
        ({"pattern": "^a"}, "string"),
        ({"format": "date"}, "string"),
        ({"minLength": 1}, "string"),
        ({"minimum": 0}, "number"),
        ({"multipleOf": 0.01}, "number"),
        ({"description": "no constraining keywords"}, None),
        ({}, None),
    ],
)
def test_infer_implicit_type_pass_various_keyword_sets(schema: dict[str, Any], expected: str | None) -> None:
    """_infer_implicit_type() infers object/array/string/number from type-specific keywords, or None if unconstrained."""
    assert _infer_implicit_type(schema) == expected


def test_infer_implicit_type_pass_object_keywords_take_priority_over_array() -> None:
    """_infer_implicit_type() prefers 'object' over 'array' when a schema (unusually) mixes both keyword groups."""
    schema = {"properties": {"a": {"type": "string"}}, "items": {"type": "string"}}
    assert _infer_implicit_type(schema) == "object"


def test_infer_implicit_type_pass_array_keywords_take_priority_over_string() -> None:
    """_infer_implicit_type() prefers 'array' over 'string' when a schema mixes both keyword groups."""
    schema = {"items": {"type": "string"}, "pattern": "^a"}
    assert _infer_implicit_type(schema) == "array"


def test_infer_implicit_type_pass_string_keywords_take_priority_over_number() -> None:
    """_infer_implicit_type() prefers 'string' over 'number' when a schema mixes both keyword groups."""
    schema = {"pattern": "^a", "minimum": 0}
    assert _infer_implicit_type(schema) == "string"


def test_infer_implicit_type_pass_never_infers_integer() -> None:
    """_infer_implicit_type() infers 'number', never 'integer', from numeric keywords alone."""
    assert _infer_implicit_type({"minimum": 0, "maximum": 10}) == "number"


"""get_known_jsonschema_keywords / _metadata_keys_for"""


def test_get_known_jsonschema_keywords_pass_excludes_ref_and_identity_keywords() -> None:
    """get_known_jsonschema_keywords() never includes $ref/$id/$schema-style identity keywords."""
    keywords = get_known_jsonschema_keywords(Draft202012Validator)
    assert "$ref" not in keywords
    assert "$id" not in keywords
    assert "$schema" not in keywords


def test_get_known_jsonschema_keywords_pass_includes_annotation_and_assertion_keywords() -> None:
    """get_known_jsonschema_keywords() includes both annotation and draft-specific assertion keywords."""
    keywords = get_known_jsonschema_keywords(Draft202012Validator)
    assert "description" in keywords
    assert "pattern" in keywords
    assert "unevaluatedProperties" in keywords


def test_metadata_keys_for_pass_excludes_structural_keywords() -> None:
    """_metadata_keys_for() excludes structural/compositional keywords like 'properties'/'type'."""
    keys = _metadata_keys_for(Draft202012Validator)
    assert "properties" not in keys
    assert "type" not in keys
    assert "description" in keys


def test_metadata_keys_for_pass_cached_per_validator_class() -> None:
    """_metadata_keys_for() returns the identical cached object for repeated calls with the same class."""
    assert _metadata_keys_for(Draft7Validator) is _metadata_keys_for(Draft7Validator)


def test_get_known_jsonschema_keywords_pass_draft7_includes_annotation_keywords() -> None:
    """get_known_jsonschema_keywords() includes Draft-07's own annotation/assertion keywords for Draft7Validator."""
    keywords = get_known_jsonschema_keywords(Draft7Validator)
    assert "description" in keywords
    assert "pattern" in keywords
    assert "$ref" not in keywords


"""_merge_with_builtin_format_map"""


def test_merge_with_builtin_format_map_pass_merges_user_overrides_over_defaults() -> None:
    """_merge_with_builtin_format_map() merges a non-empty user mapping on top of DEFAULT_FORMAT_MAP."""
    result = JSONSchemaToPySpark._merge_with_builtin_format_map({"date-time": StringType()})
    assert result["date-time"] == StringType()
    assert result["date"] == DateType()  # untouched default entry survives


@pytest.mark.parametrize("value", [None, {}, "", 0])
def test_merge_with_builtin_format_map_pass_falsy_value_passed_through_unchanged(
    value: None | dict | str | int,
) -> None:
    """_merge_with_builtin_format_map() passes falsy values straight through without merging."""
    assert JSONSchemaToPySpark._merge_with_builtin_format_map(value) == value


def test_merge_with_builtin_format_map_pass_non_dict_value_passed_through_unchanged() -> None:
    """_merge_with_builtin_format_map() passes non-dict/non-frozendict values through unchanged."""
    sentinel = object()
    assert JSONSchemaToPySpark._merge_with_builtin_format_map(sentinel) is sentinel


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ({"a", "b"}, frozenset({"a", "b"})),
        (["a", "b"], frozenset({"a", "b"})),
        (("a", "b"), frozenset({"a", "b"})),
        (frozenset({"a", "b"}), frozenset({"a", "b"})),
    ],
)
def test_coerce_extra_metadata_keywords_pass_coerces_set_list_tuple_to_frozenset(
    value: list | tuple | set | frozenset, expected: frozenset[str]
) -> None:
    """_coerce_extra_metadata_keywords() coerces set/list/tuple inputs (and passes frozenset through) unchanged in value."""
    assert JSONSchemaToPySpark._coerce_extra_metadata_keywords(value) == expected


def test_coerce_extra_metadata_keywords_pass_frozenset_identity_preserved() -> None:
    """_coerce_extra_metadata_keywords() returns the exact same frozenset object when already a frozenset."""
    value = frozenset({"a"})
    assert JSONSchemaToPySpark._coerce_extra_metadata_keywords(value) is value


def test_coerce_extra_metadata_keywords_pass_other_type_passed_through_unchanged() -> None:
    """_coerce_extra_metadata_keywords() passes non-set/list/tuple/frozenset values through unchanged."""
    sentinel = object()
    assert JSONSchemaToPySpark._coerce_extra_metadata_keywords(sentinel) is sentinel


@pytest.mark.parametrize("value", [["a"], ("a",), {"a"}])
def test_json_schema_to_pyspark_pass_extra_metadata_keywords_accepts_list_tuple_set(
    value: list[str] | tuple[str] | set[str],
) -> None:
    """JSONSchemaToPySpark(extra_metadata_keywords=...) accepts list/tuple/set and stores a frozenset."""
    converter = JSONSchemaToPySpark(extra_metadata_keywords=value)  # pyright: ignore[reportArgumentType]
    assert converter.extra_metadata_keywords == frozenset({"a"})


"""End-to-end test covering every JSON Schema construct the converter supports,
in a single schema, asserting the exact resulting StructType.
"""


def _all_types_schema() -> dict[str, Any]:
    """Build a schema covering every JSON Schema construct handled by the converter."""
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "title": "AllTypes",
        "properties": {
            # --- primitives ---
            "a_string": {"type": "string", "title": "A String", "description": "a plain string"},
            "a_formatted_date": {"type": "string", "format": "date"},
            "a_formatted_datetime": {"type": "string", "format": "date-time"},
            "a_formatted_uuid": {"type": "string", "format": "uuid"},
            "a_small_integer": {"type": "integer", "minimum": 0, "maximum": 100},
            "a_large_integer": {"type": "integer", "minimum": 0, "maximum": 99_999_999_999},
            "an_unbounded_integer": {"type": "integer"},
            "a_number": {"type": "number"},
            "a_decimal_number": {"type": "number", "multipleOf": 0.01},
            "a_boolean": {"type": "boolean"},
            "a_null": {"type": "null"},
            # --- enum (typeless) ---
            "a_string_enum": {"enum": ["red", "green", "blue"]},
            "an_integer_enum": {"enum": [1, 2, 3]},
            "a_boolean_enum": {"enum": [True, False]},
            # --- nullable union ---
            "a_nullable_string": {"type": ["string", "null"]},
            # --- array ---
            "a_string_array": {"type": "array", "items": {"type": "string"}},
            "a_nested_array": {
                "type": "array",
                "items": {"type": "array", "items": {"type": "integer"}},
            },
            # --- nested object (fixed properties) ---
            "a_nested_object": {
                "type": "object",
                "properties": {
                    "inner_field": {"type": "string"},
                },
                "required": ["inner_field"],
            },
            # --- object -> MapType (dynamic keys) ---
            "a_string_map": {"type": "object", "additionalProperties": {"type": "string"}},
            "a_pattern_property_map": {
                "type": "object",
                "patternProperties": {"^x-": {"type": "integer"}},
            },
            # --- oneOf / anyOf approximated by first branch ---
            "a_one_of_field": {"oneOf": [{"type": "integer"}, {"type": "string"}]},
            "an_any_of_field": {"anyOf": [{"type": "boolean"}, {"type": "string"}]},
        },
        "required": [
            "a_string",
            "a_small_integer",
        ],
    }


def test_convert_e2e_pass_all_supported_types_produce_expected_structtype() -> None:
    """A schema exercising every supported JSON Schema construct converts to the exact expected StructType."""
    schema = _all_types_schema()
    dereferenced = dereference_schema(schema)
    result = JSONSchemaToPySpark().convert(dereferenced)

    expected = StructType(
        [
            StructField(
                "a_string",
                StringType(),
                nullable=False,
                metadata={"jsonschema": {"title": "A String"}, "comment": "a plain string"},
            ),
            StructField("a_formatted_date", DateType(), nullable=True, metadata={}),
            StructField("a_formatted_datetime", TimestampType(), nullable=True, metadata={}),
            StructField("a_formatted_uuid", StringType(), nullable=True, metadata={}),
            StructField("a_small_integer", IntegerType(), nullable=False, metadata={}),
            StructField("a_large_integer", LongType(), nullable=True, metadata={}),
            StructField("an_unbounded_integer", LongType(), nullable=True, metadata={}),
            StructField("a_number", DoubleType(), nullable=True, metadata={}),
            StructField("a_decimal_number", DecimalType(38, 2), nullable=True, metadata={}),
            StructField("a_boolean", BooleanType(), nullable=True, metadata={}),
            StructField("a_null", NullType(), nullable=True, metadata={}),
            StructField("a_string_enum", StringType(), nullable=True, metadata={}),
            StructField("an_integer_enum", LongType(), nullable=True, metadata={}),
            StructField("a_boolean_enum", BooleanType(), nullable=True, metadata={}),
            StructField("a_nullable_string", StringType(), nullable=True, metadata={}),
            StructField("a_string_array", ArrayType(StringType(), containsNull=True), nullable=True, metadata={}),
            StructField(
                "a_nested_array",
                ArrayType(ArrayType(LongType(), containsNull=True), containsNull=True),
                nullable=True,
                metadata={},
            ),
            StructField(
                "a_nested_object",
                StructType([StructField("inner_field", StringType(), nullable=False, metadata={})]),
                nullable=True,
                metadata={},
            ),
            StructField(
                "a_string_map", MapType(StringType(), StringType(), valueContainsNull=True), nullable=True, metadata={}
            ),
            StructField(
                "a_pattern_property_map",
                MapType(StringType(), LongType(), valueContainsNull=True),
                nullable=True,
                metadata={},
            ),
            StructField("a_one_of_field", LongType(), nullable=True, metadata={}),
            StructField("an_any_of_field", BooleanType(), nullable=True, metadata={}),
        ]
    )
    assert result == expected


def test_convert_e2e_pass_all_supported_types_field_names_match_schema_properties() -> None:
    """The resulting StructType's field names exactly match the schema's 'properties' keys, in order."""
    schema = _all_types_schema()
    result = JSONSchemaToPySpark().convert(schema)
    assert result.fieldNames() == list(schema["properties"].keys())


def test_convert_e2e_pass_all_supported_types_required_fields_are_non_nullable() -> None:
    """Every field listed in the schema's root 'required' array is non-nullable in the result; all others are nullable."""
    schema = _all_types_schema()
    result = JSONSchemaToPySpark().convert(schema)
    required = set(schema["required"])
    for field in result.fields:
        assert field.nullable is (field.name not in required)


def test_build_context_pass_threads_converters_extra_metadata_keywords() -> None:
    """_build_context() produces a ConversionContext whose extra_metadata_keywords match the converter's own."""
    converter = JSONSchemaToPySpark(extra_metadata_keywords={"x-pii", "pattern"})  # pyright: ignore[reportArgumentType]
    ctx = converter._build_context(Draft202012Validator)
    assert ctx.extra_metadata_keywords == converter.extra_metadata_keywords


def test_convert_e2e_pass_extra_metadata_keywords_take_effect_via_convert() -> None:
    """A converter's extra_metadata_keywords reach _build_metadata correctly when driven through convert()."""
    converter = JSONSchemaToPySpark(extra_metadata_keywords={"x-pii"})  # pyright: ignore[reportArgumentType]
    schema = base_object_schema(properties={"ssn": {"type": "string", "x-pii": True}})
    result = converter.convert(schema)
    assert result["ssn"].metadata == {"jsonschema": {"x-pii": True}}


def test_convert_e2e_pass_property_without_type_keyword_infers_object() -> None:
    """A nested property schema that omits 'type' but declares 'properties' still converts to a StructType."""
    schema = {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "properties": {
            "address": {
                "properties": {"city": {"type": "string"}},
                "required": ["city"],
            },
        },
    }
    result = JSONSchemaToPySpark().convert(schema)
    assert result["address"].dataType == StructType([StructField("city", StringType(), nullable=False, metadata={})])
