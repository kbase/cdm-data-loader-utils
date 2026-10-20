"""Unit tests for pyiceberg_to_linkml."""

import pytest
from linkml_runtime.dumpers import yaml_dumper
from linkml_runtime.linkml_model.meta import SchemaDefinition, SlotDefinition
from linkml_runtime.loaders import yaml_loader
from pyiceberg.schema import Schema
from pyiceberg.types import (
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
    TimestampType,
    TimestamptzType,
    TimeType,
    UUIDType,
)

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.pyiceberg_to_linkml import (
    convert_field,
    convert_struct,
    convert_type,
    table_to_linkml,
)
from tests.cdm_data_loaders.converters.conftest import make_table


@pytest.mark.parametrize(
    ("field_type", "expected_range"),
    [
        pytest.param(BooleanType(), "boolean", id="boolean"),
        pytest.param(IntegerType(), "integer", id="integer"),
        pytest.param(LongType(), "integer", id="long"),
        pytest.param(FloatType(), "float", id="float"),
        pytest.param(DoubleType(), "float", id="double"),
        pytest.param(StringType(), "string", id="string"),
        pytest.param(DateType(), "date", id="date"),
        pytest.param(TimeType(), "time", id="time"),
        pytest.param(UUIDType(), "string", id="uuid"),
        pytest.param(TimestamptzType(), "datetime", id="timestamptz"),
        pytest.param(TimestampType(), "datetime", id="timestamp"),
        pytest.param(DecimalType(10, 2), "decimal", id="decimal"),
        pytest.param(FixedType(16), "string", id="fixed"),
    ],
)
def test_convert_type_pass_known_scalars(field_type: IcebergType, expected_range: str) -> None:
    """Every scalar Iceberg type maps to its expected LinkML range."""
    result = convert_type(field_type)
    assert result == SlotDefinition(name="node", range=expected_range)


def test_convert_type_fail_unknown_type() -> None:
    """An unknown Iceberg type raises ConversionError."""
    with pytest.raises(ConversionError, match="No LinkML mapping"):
        convert_type(object())  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("required", [True, False], ids=["required-elements", "optional-elements"])
def test_convert_type_pass_list(required: bool) -> None:
    """Lists map to official multivalued slot definitions."""
    node = convert_type(ListType(element_id=2, element=StringType(), element_required=required))
    assert node == SlotDefinition(name="node", range="string", multivalued=True)


def test_convert_type_pass_required_value_map() -> None:
    """A map maps to a custom class with key/value attributes."""
    result = convert_type(
        MapType(key_id=2, key_type=StringType(), value_id=3, value_type=LongType(), value_required=True)
    )
    assert result == SlotDefinition(name="node", range="Node_node_map", multivalued=True)


@pytest.mark.parametrize("required", [True, False], ids=["required-field", "optional-field"])
def test_convert_field_pass_cardinality_description(required: bool) -> None:
    """Fields preserve names, requiredness, and descriptions in slot models."""
    field = NestedField(1, "name", StringType(), required=required, doc="Display name.")
    result = convert_field(field)
    assert result == SlotDefinition(
        name="name", range="string", required=True if required else None, description="Display name."
    )


def test_convert_struct_pass_class_range() -> None:
    """Struct conversion returns a slot referencing the generated class."""
    result = convert_struct((NestedField(1, "name", StringType()),))
    assert result == SlotDefinition(name="node", range="Node_node")


def test_table_to_linkml_pass_simple() -> None:
    """A simple table produces a LinkML model that round-trips through YAML."""
    schema = Schema(
        identifier_field_ids=[],
        fields=(NestedField(1, "id", StringType(), required=True),),
    )
    table = make_table(
        schema=schema,
        identifier=("namespace", "table"),
    )
    result = table_to_linkml(table, ("namespace", "table"))
    assert isinstance(result, SchemaDefinition)
    assert result.name == "table"
    assert result.id == "urn:iceberg:namespace.table"
    assert set(result.classes) == {"table"}
    assert result.classes["table"].attributes == {"id": SlotDefinition(name="id", range="string", required=True)}
    assert yaml_loader.loads(yaml_dumper.dumps(result), target_class=SchemaDefinition) == result
