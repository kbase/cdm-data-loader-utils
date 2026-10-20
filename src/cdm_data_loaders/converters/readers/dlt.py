"""Build typed schema trees from dlt's normalized table forest."""

import logging
from collections.abc import Mapping
from dataclasses import dataclass, replace
from typing import Any, Final, Literal

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry, Extensions
from cdm_data_loaders.converters.core.ir import Field, NodeHints, NodeType, Provenance, SchemaDocument, TypedNode
from cdm_data_loaders.utils.dlt.dlt_normalization import DltTableNode, unflatten_tables

logger = logging.getLogger(__name__)
DLT_INTERNAL_PREFIX: Final = "_dlt_"
TYPES: Final[dict[str, NodeType]] = {
    "text": "string",
    "bigint": "integer",
    "double": "number",
    "decimal": "number",
    "bool": "boolean",
    "timestamp": "string",
    "date": "string",
    "time": "string",
    "binary": "string",
    "json": "any",
    "wei": "integer",
}


@dataclass(frozen=True, slots=True)
class DltReader:
    """Resolve filtering and ambiguous child-table shape at the reader boundary."""

    include_dlt_columns: bool = False
    include_variant_columns: bool = False
    child_table_mode: Literal["object", "array"] = "object"
    extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS

    def __post_init__(self) -> None:
        """Reject invalid policy configuration."""
        if self.child_table_mode not in {"object", "array"}:
            msg = "child_table_mode must be 'object' or 'array'"
            raise ValueError(msg)
        if type(self.include_dlt_columns) is not bool or type(self.include_variant_columns) is not bool:
            msg = "Column filters must be booleans"
            raise TypeError(msg)

    def read(self, source: Mapping[str, Any]) -> dict[str, SchemaDocument]:
        """Return ordered root documents from stored-schema or plain-table input."""
        if not isinstance(source, Mapping):
            msg = "Expected a dlt stored schema or tables mapping"
            raise ConversionError(msg)
        if "tables" in source and isinstance(source["tables"], Mapping):
            tables = source["tables"]
        elif all(isinstance(value, Mapping) and ("columns" in value or "parent" in value) for value in source.values()):
            tables = source
        else:
            msg = (
                "Input must be a dlt stored schema (with a 'tables' dict) or a "
                "{table_name: TTableSchema} mapping; no tables found."
            )
            raise ConversionError(msg)
        if not tables:
            msg = "The schema contains no tables."
            raise ConversionError(msg)
        roots = unflatten_tables(tables)
        envelope = {key: value for key, value in source.items() if key != "tables"} if "tables" in source else {}
        documents = {}
        for name, table in roots.items():
            root = self._table(table)
            documents[name] = SchemaDocument(
                root=root,
                name=name,
                provenance=Provenance(
                    "dlt",
                    (name,),
                    {
                        "envelope": envelope,
                        "child_table_mode": self.child_table_mode,
                        "include_dlt_columns": self.include_dlt_columns,
                        "include_variant_columns": self.include_variant_columns,
                    },
                ),
            )
        return documents

    def _columns(self, table: DltTableNode) -> dict[str, Any]:
        """Filter before scalar-value detection and required propagation."""
        return {
            name: column
            for name, column in table.table.get("columns", {}).items()
            if (self.include_dlt_columns or not name.startswith(DLT_INTERNAL_PREFIX))
            and (self.include_variant_columns or "__v_" not in name)
        }

    def _table(self, table: DltTableNode) -> TypedNode:
        """Build real child nodes while retaining isolated source table metadata.

        A child table's name always wins over a same-named scalar column
        (dlt's own naming can produce this when a column is later normalized
        into a nested table); the scalar column's `required` flag still
        contributes to the resulting field's `required` flag.
        """
        columns = self._columns(table)
        properties = {name: self._column(name, column, table.name) for name, column in columns.items()}
        required_names = [name for name, column in columns.items() if column.get("nullable", True) is False]
        for name, child in table.children.items():
            child_node = self._table(child)
            required = any(column.get("nullable", True) is False for column in self._columns(child).values())
            properties[name] = Field(name, child_node, required=required or name in required_names)
            if required:
                required_names.append(name)
        annotations = {"description": table.table["description"]} if "description" in table.table else {}
        provenance = Provenance("dlt", ("tables", table.name), table.table)
        object_node = TypedNode(
            type="object",
            nullable=False,
            properties=tuple(properties.values()),
            required_names=tuple(required_names),
            additional_properties=TypedNode(type="never", nullable=False),
            annotations=annotations,
            provenance=provenance,
        )
        if table.parent is None:
            return object_node
        scalar = not table.children and list(columns) == ["value"]
        if scalar or self.child_table_mode == "array":
            return TypedNode(
                type="array",
                nullable=False,
                items=properties["value"].node if scalar else replace(object_node, annotations={}),
                annotations=annotations,
                provenance=provenance,
            )
        return object_node

    def _column(self, name: str, column: Mapping[str, Any], table_name: str) -> Field:
        """Preserve source type and hints without adopting JSON target encodings."""
        data_type = column.get("data_type")
        kind = TYPES.get(data_type, "any" if data_type is None else "unknown")
        if kind == "unknown":
            logger.warning("Unknown dlt data_type %r at %s.%s retained in IR", data_type, table_name, name)
        provenance = Provenance("dlt", ("tables", table_name, "columns", name), column)
        nullable = column.get("nullable", True)
        hints = NodeHints(
            logical_type=data_type,
            precision=column.get("precision"),
            scale=column.get("scale"),
            bit_width=column.get("precision") if data_type == "bigint" else None,
            timezone=column.get("timezone"),
        )
        extensions = Extensions({"x-dlt": dict(column)}, self.extension_registry)
        node = TypedNode(
            type=kind,
            nullable=nullable,
            hints=hints,
            annotations={"description": column["description"]} if "description" in column else {},
            extensions=extensions,
            provenance=provenance,
        )
        return Field(name, node, required=nullable is False, provenance=provenance)
