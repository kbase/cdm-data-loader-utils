"""Convert Iceberg types and loaded tables through typed readers and emitters."""

from linkml_runtime.linkml_model.meta import SchemaDefinition, SlotDefinition
from pyiceberg.table import Table
from pyiceberg.types import IcebergType, NestedField, StructType

from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry
from cdm_data_loaders.converters.emitters.linkml import LinkMLEmitter
from cdm_data_loaders.converters.readers.iceberg import IcebergReader


def convert_type(
    field_type: IcebergType, *, extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS
) -> SlotDefinition:
    """Read and emit a single Iceberg type.

    :param field_type: Iceberg type to convert.
    :param extension_registry: Registry for Iceberg extensions.
    :returns: LinkML slot definition.
    """
    return LinkMLEmitter().emit_node(IcebergReader(extension_registry, mapping_target="LinkML").read_node(field_type))


def convert_field(field: NestedField, *, extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS) -> SlotDefinition:
    """Read and emit a field with its cardinality and description.

    :param field: Iceberg field to convert.
    :param extension_registry: Registry for Iceberg extensions.
    :returns: Named LinkML slot definition.
    """
    return LinkMLEmitter().emit_field(IcebergReader(extension_registry, mapping_target="LinkML").read_field(field))


def convert_struct(
    fields: tuple[NestedField, ...], *, extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS
) -> SlotDefinition:
    """Read and emit a struct without a document envelope.

    :param fields: Iceberg struct fields to convert.
    :param extension_registry: Registry for Iceberg extensions.
    :returns: LinkML slot referencing the struct class.
    """
    return convert_type(StructType(*fields), extension_registry=extension_registry)


def table_to_linkml(
    table: Table,
    identifier: tuple[str, ...],
    *,
    extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS,
    schema: SchemaDefinition | None = None,
) -> SchemaDefinition:
    """Read loaded table metadata and emit a LinkML schema model without I/O.

    :param table: Loaded Iceberg table.
    :param identifier: Namespace and table name.
    :param extension_registry: Registry for Iceberg extensions.
    :param schema: Optional schema to append table classes to, preserving existing definitions.
    :returns: Schema model suitable for LinkML's YAML dumper.
    """
    document = IcebergReader(extension_registry, mapping_target="LinkML").read_table(table, identifier)
    emitter = LinkMLEmitter()
    return emitter.emit(document) if schema is None else emitter.emit_into(document, schema)
