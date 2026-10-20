"""Unit tests for pyiceberg_to_jsonschema."""

import re

import pytest
from jsonschema import Draft202012Validator, FormatChecker
from pyiceberg.partitioning import PartitionField, PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.transforms import IdentityTransform
from pyiceberg.types import (
    BooleanType,
    DateType,
    DecimalType,
    DoubleType,
    FixedType,
    FloatType,
    IntegerType,
    ListType,
    LongType,
    MapType,
    NestedField,
    StringType,
    StructType,
    TimestampType,
    TimestamptzType,
    TimeType,
    UUIDType,
)

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.emitters.json_schema import JSON_SCHEMA_DIALECT, decimal_pattern
from cdm_data_loaders.converters.pyiceberg_to_jsonschema import (
    convert_field,
    convert_struct,
    convert_type,
    table_to_json_schema,
)
from tests.cdm_data_loaders.converters.conftest import (
    JSON_SCHEMA_KEYWORDS,
    iter_extension_keys,
    iter_schema_keywords,
    make_table,
    to_snake_case,
)


def test_decimal_pattern_with_scale() -> None:
    """A scaled decimal produces an integer-or-fractional pattern."""
    pattern = decimal_pattern(10, 2)
    assert pattern == r"^-?\d{1,8}(\.\d{1,2})?$"
    assert re.match(pattern, "123.45")
    assert re.match(pattern, "-123.45")
    assert not re.match(pattern, "123.456")


def test_decimal_pattern_zero_scale() -> None:
    """A zero-scale decimal produces an integer-only pattern."""
    assert decimal_pattern(10, 0) == r"^-?\d{1,10}$"


def test_decimal_pattern_precision_equals_scale_raises() -> None:
    """A scale equal to precision produces a leading-zero pattern that compiles."""
    pattern = decimal_pattern(2, 2)
    assert pattern == r"^-?0(\.\d{1,2})?$"
    re.compile(pattern)


@pytest.mark.parametrize(
    ("field_type", "expected"),
    [
        pytest.param(BooleanType(), {"type": "boolean"}, id="boolean"),
        pytest.param(
            IntegerType(),
            {"type": "integer", "minimum": -2147483648, "maximum": 2147483647},
            id="integer",
        ),
        pytest.param(
            LongType(),
            {"type": "integer", "minimum": -9223372036854775808, "maximum": 9223372036854775807},
            id="long",
        ),
        pytest.param(FloatType(), {"type": "number"}, id="float"),
        pytest.param(DoubleType(), {"type": "number"}, id="double"),
        pytest.param(StringType(), {"type": "string"}, id="string"),
        pytest.param(DateType(), {"type": "string", "format": "date"}, id="date"),
        pytest.param(TimeType(), {"type": "string", "format": "time"}, id="time"),
        pytest.param(UUIDType(), {"type": "string", "format": "uuid"}, id="uuid"),
        pytest.param(TimestamptzType(), {"type": "string", "format": "date-time"}, id="timestamptz"),
        pytest.param(
            DecimalType(10, 2),
            {
                "type": "string",
                "pattern": r"^-?\d{1,8}(\.\d{1,2})?$",
                "x-iceberg": {"logical_type": "decimal", "precision": 10, "scale": 2},
            },
            id="decimal",
        ),
        pytest.param(
            TimestampType(),
            {
                "type": "string",
                "description": "ISO 8601 timestamp with no timezone attached.",
                "x-iceberg": {"logical_type": "timestamp-without-tz"},
            },
            id="timestamp",
        ),
        pytest.param(
            FixedType(16),
            {
                "type": "string",
                "contentEncoding": "base64",
                "x-iceberg": {"logical_type": "fixed", "length": 16},
            },
            id="fixed",
        ),
    ],
)
def test_convert_type_pass_known_scalars(field_type: object, expected: dict) -> None:
    """Every scalar Iceberg type maps to its expected JSON Schema node."""
    assert convert_type(field_type) == expected  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize(
    ("field_type", "valid_instance", "invalid_instance"),
    [
        pytest.param(BooleanType(), True, "not-a-boolean", id="boolean"),
        pytest.param(IntegerType(), 5, 2.5, id="integer"),
        pytest.param(LongType(), 5, 9223372036854775808, id="long"),
        pytest.param(FloatType(), 1.5, "not-a-number", id="float"),
        pytest.param(DoubleType(), 1.5, "not-a-number", id="double"),
        pytest.param(StringType(), "s", 5, id="string"),
        pytest.param(DateType(), "2026-09-22", "not-a-date", id="date"),
        pytest.param(TimeType(), "12:30:45", "not-a-time", id="time"),
        pytest.param(UUIDType(), "550e8400-e29b-41d4-a716-446655440000", "not-a-uuid", id="uuid"),
        pytest.param(TimestamptzType(), "2026-09-22T00:00:00+00:00", "not-a-datetime", id="timestamptz"),
        pytest.param(DecimalType(10, 2), "123.45", 123.45, id="decimal"),
        pytest.param(FixedType(16), "base64bytes", 5, id="fixed"),
    ],
)
def test_convert_type_pass_scalars_are_valid_draft_2020_12(
    field_type: object, valid_instance: object, invalid_instance: object
) -> None:
    """Every scalar mapping produces a schema that meta-validates and rejects bad instances."""
    doc = convert_struct((NestedField(1, "value", field_type, required=True),))  # pyright: ignore[reportArgumentType]
    Draft202012Validator.check_schema(doc)
    validator = Draft202012Validator(doc, format_checker=FormatChecker())
    assert validator.is_valid({"value": valid_instance})  # pyright: ignore[reportArgumentType]
    assert not validator.is_valid({"value": invalid_instance})  # pyright: ignore[reportArgumentType]


def test_convert_type_fail_unknown_type() -> None:
    """An unknown Iceberg type raises ConversionError."""
    with pytest.raises(ConversionError, match="No JSON Schema mapping"):
        convert_type(object())  # pyright: ignore[reportArgumentType]


def test_convert_type_pass_required_list() -> None:
    """A required-element list maps to a plain array."""
    node = convert_type(ListType(element_id=2, element=StringType(), element_required=True))
    assert node == {"type": "array", "items": {"type": "string"}}


def test_convert_type_pass_optional_element_list() -> None:
    """An optional-element list wraps items in anyOf with null."""
    node = convert_type(ListType(element_id=2, element=StringType(), element_required=False))
    assert node == {"type": "array", "items": {"anyOf": [{"type": "string"}, {"type": "null"}]}}


def test_convert_type_pass_required_value_map() -> None:
    """A required-value map maps to key/value pair objects."""
    node = convert_type(
        MapType(key_id=2, key_type=StringType(), value_id=3, value_type=LongType(), value_required=True)
    )
    assert node == {
        "type": "array",
        "description": "Iceberg map rendered as key/value pairs (JSON object keys must be strings).",
        "items": {
            "type": "object",
            "properties": {
                "key": {"type": "string"},
                "value": {
                    "type": "integer",
                    "minimum": -9223372036854775808,
                    "maximum": 9223372036854775807,
                },
            },
            "required": ["key", "value"],
            "additionalProperties": False,
        },
    }


def test_convert_type_pass_optional_value_map() -> None:
    """An optional-value map wraps the value schema in anyOf with null."""
    node = convert_type(
        MapType(key_id=2, key_type=StringType(), value_id=3, value_type=LongType(), value_required=False)
    )
    assert node["items"]["properties"]["value"] == {
        "anyOf": [
            {
                "type": "integer",
                "minimum": -9223372036854775808,
                "maximum": 9223372036854775807,
            },
            {"type": "null"},
        ]
    }


def test_convert_type_pass_nested_struct() -> None:
    """A struct type recurses into convert_struct."""
    node = convert_type(
        StructType(
            NestedField(1, "inner", StringType(), required=True),
            NestedField(2, "count", LongType(), required=False),
        )
    )
    assert node == {
        "type": "object",
        "properties": {
            "inner": {"type": "string", "x-iceberg": {"field_id": 1, "required": True}},
            "count": {
                "anyOf": [
                    {
                        "type": "integer",
                        "minimum": -9223372036854775808,
                        "maximum": 9223372036854775807,
                    },
                    {"type": "null"},
                ],
                "x-iceberg": {"field_id": 2, "required": False},
            },
        },
        "required": ["inner"],
        "additionalProperties": False,
    }


def test_convert_field_pass_required_with_doc() -> None:
    """A required documented field keeps its schema and carries iceberg metadata."""
    node = convert_field(NestedField(7, "id", LongType(), required=True, doc="primary key"))
    assert node == {
        "type": "integer",
        "minimum": -9223372036854775808,
        "maximum": 9223372036854775807,
        "description": "primary key",
        "x-iceberg": {"field_id": 7, "required": True},
    }


def test_convert_field_pass_optional_wraps_anyof() -> None:
    """An optional field wraps its schema in anyOf with null and keeps iceberg metadata at the top level."""
    node = convert_field(NestedField(3, "note", StringType(), required=False))
    assert node == {
        "anyOf": [{"type": "string"}, {"type": "null"}],
        "x-iceberg": {"field_id": 3, "required": False},
    }


def test_convert_field_pass_defaults_added_to_iceberg_meta() -> None:
    """initial_default and write_default are recorded in the iceberg metadata."""
    node = convert_field(
        NestedField(5, "status", StringType(), required=True, initial_default="new", write_default="new")
    )
    assert node["x-iceberg"] == {"field_id": 5, "required": True, "initial_default": "new", "write_default": "new"}


def test_convert_field_pass_defaults_omitted_when_none() -> None:
    """Defaults that are None are omitted from the iceberg metadata."""
    node = convert_field(NestedField(5, "status", StringType(), required=True))
    assert node["x-iceberg"] == {"field_id": 5, "required": True}


def test_convert_struct_pass_mixed_required() -> None:
    """convert_struct lists only required fields in the required array."""
    node = convert_struct(
        (
            NestedField(1, "a", StringType(), required=True),
            NestedField(2, "b", StringType(), required=False),
            NestedField(3, "c", LongType(), required=True),
        )
    )
    assert node == {
        "type": "object",
        "properties": {
            "a": {"type": "string", "x-iceberg": {"field_id": 1, "required": True}},
            "b": {"anyOf": [{"type": "string"}, {"type": "null"}], "x-iceberg": {"field_id": 2, "required": False}},
            "c": {
                "type": "integer",
                "minimum": -9223372036854775808,
                "maximum": 9223372036854775807,
                "x-iceberg": {"field_id": 3, "required": True},
            },
        },
        "required": ["a", "c"],
        "additionalProperties": False,
    }


def test_convert_struct_pass_no_required() -> None:
    """A struct with no required fields omits the required key entirely."""
    node = convert_struct((NestedField(1, "a", StringType(), required=False),))
    assert node == {
        "type": "object",
        "properties": {
            "a": {"anyOf": [{"type": "string"}, {"type": "null"}], "x-iceberg": {"field_id": 1, "required": False}}
        },
        "additionalProperties": False,
    }


def test_convert_struct_pass_empty() -> None:
    """An empty struct produces an object with no properties or required keys."""
    node = convert_struct(())
    assert node == {"type": "object", "properties": {}, "additionalProperties": False}


def test_table_to_json_schema_pass_minimal_table() -> None:
    """A table without partition fields or comment produces the expected envelope."""
    iceberg_format_version_v2 = 2
    table = make_table(
        ("ns", "tbl"),
        Schema(
            NestedField(1, "id", LongType(), required=True),
            NestedField(2, "name", StringType()),
        ),
    )
    doc = table_to_json_schema(table, ("ns", "tbl"))
    assert doc["$schema"] == JSON_SCHEMA_DIALECT
    assert doc["$id"] == "urn:iceberg:ns.tbl"
    assert doc["title"] == "tbl"
    assert "description" not in doc
    assert doc["x-iceberg"]["identifier"] == ["ns", "tbl"]
    assert doc["x-iceberg"]["partition_spec"] == []
    assert doc["x-iceberg"]["format_version"] == iceberg_format_version_v2
    assert doc["x-iceberg"]["current_snapshot_id"] is None
    assert doc["x-iceberg"]["identifier_field_ids"] == []
    assert doc["x-iceberg"]["properties"] == {}
    assert doc["x-iceberg"]["schema_id"] == table.schema().schema_id


def test_table_to_json_schema_pass_partitioned_with_comment() -> None:
    """Partition fields and table comments are rendered into the envelope."""
    schema = Schema(
        NestedField(1, "id", LongType(), required=True),
        NestedField(2, "day", StringType(), required=True),
    )
    spec = PartitionSpec(PartitionField(source_id=2, field_id=1000, transform=IdentityTransform(), name="day"))
    table = make_table(("ns", "partitioned"), schema, partition_spec=spec, properties={"comment": "test comment"})
    doc = table_to_json_schema(table, ("ns", "partitioned"))
    assert doc["description"] == "test comment"
    assert doc["x-iceberg"]["partition_spec"] == [{"source_id": 2, "name": "day", "transform": "identity"}]
    assert doc["x-iceberg"]["identifier_field_ids"] == []


def test_table_to_json_schema_pass_identifier_fields() -> None:
    """Identifier field ids are sorted into the envelope."""
    schema = Schema(
        NestedField(1, "id", UUIDType(), required=True),
        NestedField(2, "other", UUIDType(), required=True),
        identifier_field_ids=[2, 1],
    )
    table = make_table(("ns", "identified"), schema)
    doc = table_to_json_schema(table, ("ns", "identified"))
    assert doc["x-iceberg"]["identifier_field_ids"] == [1, 2]
    first_field_id = 1
    second_field_id = 2
    assert doc["properties"]["id"]["x-iceberg"]["field_id"] == first_field_id
    assert doc["properties"]["other"]["x-iceberg"]["field_id"] == second_field_id


def test_table_to_json_schema_pass_no_nonstandard_keys_without_x_prefix() -> None:
    """Every key in a dumped document is either a standard JSON Schema keyword or x- prefixed."""
    table = make_table(
        ("ns", "tbl"),
        Schema(
            NestedField(1, "id", LongType(), required=True),
            NestedField(2, "rate", DecimalType(10, 2), required=False),
        ),
        properties={"comment": "a comment"},
    )
    doc = table_to_json_schema(table, ("ns", "tbl"))
    for key in iter_schema_keywords(doc):
        assert key in JSON_SCHEMA_KEYWORDS or key.startswith("x-"), f"non-standard key without x- prefix: {key}"


def test_table_to_json_schema_pass_extension_keys_are_snake_case() -> None:
    """Every x- prefixed key and every key inside extension blocks is snake case."""
    table = make_table(
        ("ns", "tbl"),
        Schema(
            NestedField(1, "id", LongType(), required=True),
            NestedField(2, "rate", DecimalType(10, 2), required=False),
        ),
        properties={"comment": "a comment"},
    )
    doc = table_to_json_schema(table, ("ns", "tbl"))
    for key in iter_extension_keys(doc):
        assert key == to_snake_case(key), f"extension key is not snake case: {key}"


def test_table_to_json_schema_pass_namespace_schema_adds_table_under_defs() -> None:
    """Passing a namespace document adds the table's schema under $defs and returns that document."""
    table = make_table(("ns", "tbl"), Schema(NestedField(1, "id", LongType(), required=True)))
    namespace_doc = {"$schema": JSON_SCHEMA_DIALECT, "$id": "urn:iceberg:ns", "title": "ns"}

    result = table_to_json_schema(table, ("ns", "tbl"), schema=namespace_doc)

    assert result is namespace_doc
    assert set(namespace_doc["$defs"]) == {"tbl"}
    assert "$schema" not in namespace_doc["$defs"]["tbl"]
    assert namespace_doc["$defs"]["tbl"]["title"] == "tbl"


def test_table_to_json_schema_pass_namespace_schema_keeps_prior_tables() -> None:
    """A second call with the same namespace document keeps earlier tables' entries."""
    schema = Schema(NestedField(1, "id", LongType(), required=True))
    first_table = make_table(("ns", "first"), schema)
    second_table = make_table(("ns", "second"), schema)
    namespace_doc = {"$schema": JSON_SCHEMA_DIALECT, "$id": "urn:iceberg:ns", "title": "ns"}

    table_to_json_schema(first_table, ("ns", "first"), schema=namespace_doc)
    table_to_json_schema(second_table, ("ns", "second"), schema=namespace_doc)

    assert set(namespace_doc["$defs"]) == {"first", "second"}
