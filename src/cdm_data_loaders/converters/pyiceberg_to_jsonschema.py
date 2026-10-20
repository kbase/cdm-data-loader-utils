"""Convert Iceberg types and loaded tables through typed readers and emitters."""

from typing import Any, cast

from pyiceberg.table import Table
from pyiceberg.types import IcebergType, NestedField, StructType

from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry
from cdm_data_loaders.converters.emitters.json_schema import JsonSchemaEmitter
from cdm_data_loaders.converters.readers.iceberg import IcebergReader


def convert_type(
    field_type: IcebergType, *, extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS
) -> dict[str, Any]:
    """Read and emit a single Iceberg type."""
    return cast(
        "dict[str, Any]",
        JsonSchemaEmitter().emit_node(
            IcebergReader(extension_registry, mapping_target="JSON Schema").read_node(field_type)
        ),
    )


def convert_field(field: NestedField, *, extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS) -> dict[str, Any]:
    """Read and emit a field with its nullable wrapper and scoped metadata."""
    return cast(
        "dict[str, Any]",
        JsonSchemaEmitter().emit_field(
            IcebergReader(extension_registry, mapping_target="JSON Schema").read_field(field)
        ),
    )


def convert_struct(
    fields: tuple[NestedField, ...], *, extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS
) -> dict[str, Any]:
    """Read and emit a struct without a document envelope."""
    return convert_type(StructType(*fields), extension_registry=extension_registry)


def table_to_json_schema(
    table: Table,
    identifier: tuple[str, ...],
    *,
    extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS,
    schema: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Read loaded table metadata and emit a JSON Schema document without I/O.

    :param table: Loaded Iceberg table.
    :param identifier: Namespace and table name.
    :param extension_registry: Registry for Iceberg extensions.
    :param schema: Optional namespace document to add the table's schema into, preserving existing entries.
    :returns: JSON Schema document, or the destination namespace document when ``schema`` is given.
    """
    document = IcebergReader(extension_registry, mapping_target="JSON Schema").read_table(table, identifier)
    emitter = JsonSchemaEmitter()
    return emitter.emit(document) if schema is None else emitter.emit_into(document, schema)
