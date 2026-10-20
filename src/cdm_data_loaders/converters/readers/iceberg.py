"""Read Iceberg schema facts and loaded table metadata without I/O."""

from base64 import b64encode
from collections.abc import Mapping
from dataclasses import dataclass, replace
from datetime import date, datetime, time
from typing import Final
from uuid import UUID

from pyiceberg.schema import Schema
from pyiceberg.table import Table
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
    UUIDType,
)

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry, Extensions
from cdm_data_loaders.converters.core.ir import Field, NodeHints, NodeType, Provenance, SchemaDocument, TypedNode
from cdm_data_loaders.converters.core.ir_values import Value, freeze_value


def iceberg_value(value: object) -> Value:
    """Serialize typed source metadata explicitly, retaining Decimal as a number.

    UUID and temporal values use their textual forms; binary values use base64.
    Unsupported Python objects are rejected, not coerced through ``str``.
    """
    if isinstance(value, UUID):
        return str(value)
    if isinstance(value, (datetime, date, time)):
        return value.isoformat()
    if isinstance(value, bytes):
        return b64encode(value).decode("ascii")
    if isinstance(value, Mapping):
        return freeze_value({key: iceberg_value(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(iceberg_value(item) for item in value)
    return freeze_value(value)


def _metadata(value: Mapping[str, object]) -> Mapping[str, Value]:
    """Normalize loaded Iceberg metadata to owned IR values."""
    return {key: iceberg_value(item) for key, item in value.items()}


PRIMITIVES: Final[dict[type[IcebergType], tuple[NodeType, str, int | None, bool | None]]] = {
    BooleanType: ("boolean", "boolean", None, None),
    IntegerType: ("integer", "int", 32, None),
    LongType: ("integer", "long", 64, None),
    FloatType: ("number", "float", 32, None),
    DoubleType: ("number", "double", 64, None),
    StringType: ("string", "string", None, None),
    DateType: ("string", "date", None, None),
    TimeType: ("string", "time", None, None),
    UUIDType: ("string", "uuid", None, None),
    TimestampType: ("string", "timestamp-without-tz", None, False),
    TimestamptzType: ("string", "timestamp-with-tz", None, True),
    BinaryType: ("string", "binary", None, None),
}


@dataclass(frozen=True, slots=True)
class IcebergReader:
    """Read physical/logical types separately from JSON rendering conventions."""

    extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS
    mapping_target: str = "IR"

    def read(self, source: Schema) -> SchemaDocument:
        """Build a structural root from an Iceberg Schema."""
        if not isinstance(source, Schema):
            msg = "IcebergReader.read requires a Schema"
            raise ConversionError(msg)
        return SchemaDocument(
            root=self._type(source.as_struct(), ()),
            extensions=Extensions(
                {
                    "x-iceberg": {
                        "schema_id": source.schema_id,
                        "identifier_field_ids": sorted(source.identifier_field_ids),
                    }
                },
                self.extension_registry,
            ),
            provenance=Provenance("iceberg", metadata=_metadata(source.model_dump(mode="python", by_alias=True))),
        )

    def read_table(self, table: Table, identifier: tuple[str, ...]) -> SchemaDocument:
        """Capture a loaded table envelope using only in-memory metadata accessors."""
        if not identifier or any(not isinstance(part, str) or not part for part in identifier):
            msg = "An Iceberg table identifier must contain nonempty strings"
            raise ConversionError(msg)
        document = self.read(table.schema())
        metadata = table.metadata
        iceberg = {
            "identifier": list(identifier),
            "schema_id": table.schema().schema_id,
            "format_version": metadata.format_version,
            "current_snapshot_id": metadata.current_snapshot_id,
            "location": metadata.location,
            "properties": dict(table.properties),
            "partition_spec": [
                {"source_id": part.source_id, "name": part.name, "transform": str(part.transform)}
                for part in table.spec().fields
            ],
            "identifier_field_ids": sorted(table.schema().identifier_field_ids),
        }
        return replace(
            document,
            name=identifier[-1],
            identifier=f"urn:iceberg:{'.'.join(identifier)}",
            annotations={"description": table.properties["comment"]} if "comment" in table.properties else {},
            extensions=Extensions({"x-iceberg": iceberg}, self.extension_registry),
            provenance=Provenance("iceberg", identifier, _metadata(metadata.model_dump(mode="python", by_alias=True))),
        )

    def read_node(self, source: IcebergType) -> TypedNode:
        """Read a type without fabricating a schema document."""
        return self._type(source, ())

    def read_field(self, source: NestedField) -> Field:
        """Read a standalone field with its name as the provenance path."""
        return self._field(source, (source.name,))

    def _field(self, source: NestedField, path: tuple[str | int, ...]) -> Field:
        """Keep field defaults and descriptions outside the nullable type node."""
        metadata = {"field_id": source.field_id, "required": source.required}
        for name in ("initial_default", "write_default"):
            if name in source.model_fields_set:
                metadata[name] = iceberg_value(getattr(source, name))
        return Field(
            source.name,
            replace(self._type(source.field_type, path), nullable=not source.required),
            required=source.required,
            annotations={"description": source.doc} if "doc" in source.model_fields_set else {},
            extensions=Extensions({"x-iceberg": metadata}, self.extension_registry),
            provenance=Provenance("iceberg", path, _metadata(source.model_dump(mode="python", by_alias=True))),
        )

    def _type(self, source: IcebergType, path: tuple[str | int, ...]) -> TypedNode:
        """Translate supported Iceberg types into genuine typed child trees."""
        provenance = Provenance("iceberg", path)
        if type(source) in PRIMITIVES:
            kind, logical, width, timezone = PRIMITIVES[type(source)]
            return TypedNode(
                type=kind,
                nullable=False,
                provenance=provenance,
                hints=NodeHints(logical_type=logical, bit_width=width, timezone=timezone),
            )
        if isinstance(source, (DecimalType, FixedType)):
            hints = (
                NodeHints(logical_type="decimal", precision=source.precision, scale=source.scale)
                if isinstance(
                    source,
                    DecimalType,
                )
                else NodeHints(logical_type="fixed", length=source.root)
            )
            return TypedNode(
                type="number" if isinstance(source, DecimalType) else "string",
                nullable=False,
                hints=hints,
                provenance=provenance,
            )
        if isinstance(source, StructType):
            return TypedNode(
                type="object",
                nullable=False,
                provenance=provenance,
                properties=tuple(self._field(child, (*path, child.name)) for child in source.fields),
                additional_properties=TypedNode(type="never", nullable=False),
            )
        if isinstance(source, ListType):
            return TypedNode(
                type="array",
                nullable=False,
                provenance=provenance,
                items=replace(
                    self._type(source.element_type, (*path, "element")), nullable=not source.element_required
                ),
                extensions=Extensions(
                    {
                        "x-iceberg": {
                            "element_id": source.element_id,
                            "element_required": source.element_required,
                        }
                    },
                    self.extension_registry,
                ),
            )
        if isinstance(source, MapType):
            return TypedNode(
                type="map",
                nullable=False,
                provenance=provenance,
                key_type=self._type(source.key_type, (*path, "key")),
                value_type=replace(self._type(source.value_type, (*path, "value")), nullable=not source.value_required),
                extensions=Extensions(
                    {
                        "x-iceberg": {
                            "key_id": source.key_id,
                            "value_id": source.value_id,
                            "value_required": source.value_required,
                        }
                    },
                    self.extension_registry,
                ),
            )
        msg = f"No {self.mapping_target} mapping for Iceberg type: {source}"
        raise ConversionError(msg)
