"""Iceberg source types, scopes and loaded table envelopes."""

from datetime import UTC, date, datetime, time
from decimal import Decimal
from uuid import UUID

import pytest
from pyiceberg.partitioning import PartitionField, PartitionSpec
from pyiceberg.schema import Schema
from pyiceberg.transforms import IdentityTransform
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
    TimestampType,
    TimestamptzType,
    TimeType,
    UnknownType,
    UUIDType,
)

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.readers.iceberg import IcebergReader, iceberg_value
from tests.cdm_data_loaders.converters.conftest import make_table


@pytest.mark.parametrize(
    ("source", "kind", "logical", "width", "timezone"),
    [
        (BooleanType(), "boolean", "boolean", None, None),
        (IntegerType(), "integer", "int", 32, None),
        (LongType(), "integer", "long", 64, None),
        (FloatType(), "number", "float", 32, None),
        (DoubleType(), "number", "double", 64, None),
        (StringType(), "string", "string", None, None),
        (DateType(), "string", "date", None, None),
        (TimeType(), "string", "time", None, None),
        (UUIDType(), "string", "uuid", None, None),
        (BinaryType(), "string", "binary", None, None),
        (TimestampType(), "string", "timestamp-without-tz", None, False),
        (TimestamptzType(), "string", "timestamp-with-tz", None, True),
    ],
    ids=["boolean", "integer", "long", "float", "double", "string", "date", "time", "uuid", "binary", "naive", "aware"],
)
def test_read_pass_primitive_facts(
    source: IcebergType,
    kind: str,
    logical: str,
    width: int | None,
    timezone: bool | None,
) -> None:
    """Physical widths and timestamp timezone distinctions remain available to emitters."""
    node = IcebergReader().read(Schema(NestedField(1, "data", source, required=True))).root.properties[0].node
    assert (node.type, node.hints.logical_type, node.hints.bit_width, node.hints.timezone) == (
        kind,
        logical,
        width,
        timezone,
    )
    assert node.nullable is False


def test_read_pass_nested_types_and_metadata_scopes() -> None:
    """Maps, list elements, structs and defaults carry separate type and field facts."""
    schema = Schema(
        NestedField(1, "amount", DecimalType(20, 4), initial_default=Decimal("1.2300"), write_default=None, doc=None),
        NestedField(2, "fixed", FixedType(12)),
        NestedField(
            3, "entries", ListType(4, StructType(NestedField(5, "name", StringType())), element_required=False)
        ),
        NestedField(6, "map", MapType(7, IntegerType(), 8, StringType(), value_required=False)),
    )
    document = IcebergReader().read(schema)
    amount, fixed, entries, mapping = document.root.properties
    assert amount.node.type == "number"
    assert (amount.node.hints.precision, amount.node.hints.scale) == (20, 4)
    assert amount.node.constraints == {}
    assert amount.extensions["x-iceberg"] == {
        "field_id": 1,
        "required": False,
        "initial_default": Decimal("1.2300"),
        "write_default": None,
    }
    assert amount.annotations == {"description": None}
    assert not amount.node.extensions
    assert (fixed.node.hints.length,) == (12,)
    assert fixed.extensions["x-iceberg"]["initial_default"] is None
    assert entries.node.items.type == "object"
    assert entries.node.items.nullable is True
    assert entries.node.items.properties[0].name == "name"
    assert entries.node.extensions["x-iceberg"] == {"element_id": 4, "element_required": False}
    assert mapping.node.type == "map"
    assert mapping.node.key_type.type == "integer"
    assert mapping.node.key_type.nullable is False
    assert mapping.node.value_type.nullable is True
    assert mapping.node.extensions["x-iceberg"] == {"key_id": 7, "value_id": 8, "value_required": False}


def test_read_table_pass_envelope_without_io() -> None:
    """Real in-memory table metadata includes partition, property and snapshot facts."""
    schema = Schema(NestedField(1, "id", LongType(), required=True), identifier_field_ids=[1])
    table = make_table(
        ("ns", "table"),
        schema,
        PartitionSpec(PartitionField(1, 1000, IdentityTransform(), "id")),
        {"comment": "records"},
    )
    document = IcebergReader().read_table(table, ("ns", "table"))
    assert document.identifier == "urn:iceberg:ns.table"
    assert document.name == "table"
    assert document.annotations == {"description": "records"}
    assert document.extensions["x-iceberg"] == {
        "identifier": ("ns", "table"),
        "schema_id": table.schema().schema_id,
        "format_version": table.metadata.format_version,
        "current_snapshot_id": None,
        "location": table.metadata.location,
        "properties": {"comment": "records"},
        "partition_spec": ({"source_id": 1, "name": "id", "transform": "identity"},),
        "identifier_field_ids": (1,),
    }


def test_read_fail_invalid_source_and_identifier() -> None:
    """Invalid source objects and empty identifiers fail at the public boundary with ConversionError."""
    with pytest.raises(ConversionError):
        IcebergReader().read({})
    with pytest.raises(ConversionError, match="identifier"):
        IcebergReader().read_table(make_table(("ns", "table"), Schema()), ())
    assert IcebergReader().read(Schema()).root.properties == ()


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (Decimal("1.2300"), Decimal("1.2300")),
        (b"\x00\xff", "AP8="),
        (UUID(int=0), "00000000-0000-0000-0000-000000000000"),
        (datetime(2026, 1, 1, tzinfo=UTC), "2026-01-01T00:00:00+00:00"),
        (date(2026, 1, 1), "2026-01-01"),
        (time(1, 2, tzinfo=UTC), "01:02:00+00:00"),
        ({"nested": [None, Decimal("0.0001")]}, {"nested": (None, Decimal("0.0001"))}),
    ],
    ids=["decimal", "binary", "uuid", "timestamp", "date", "time", "nested"],
)
def test_iceberg_value_pass_precise_serialization(value: object, expected: object) -> None:
    """Typed metadata serializes explicitly without a Decimal-to-float conversion."""
    assert iceberg_value(value) == expected


def test_read_fail_unsupported_type_and_metadata() -> None:
    """Unmapped Iceberg types and arbitrary Python defaults are not guessed."""
    with pytest.raises(ConversionError, match="No IR mapping"):
        IcebergReader().read(Schema(NestedField(1, "unknown", UnknownType())))
    with pytest.raises(TypeError, match="Unsupported"):
        iceberg_value(object())
