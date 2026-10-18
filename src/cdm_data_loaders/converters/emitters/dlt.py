"""Emit dlt table graphs directly from typed schema documents."""

import logging
from collections.abc import Mapping
from dataclasses import replace
from decimal import Decimal
from typing import Any, Final, cast

import yaml
from dlt.common.data_types.typing import TDataType
from dlt.common.schema import Schema
from dlt.common.schema.typing import TColumnSchema, TStoredSchema, TTableSchema
from dlt.common.schema.utils import new_column, new_table
from pydantic import BaseModel, ConfigDict, Field, field_validator

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.inference import json_type_from_enum
from cdm_data_loaders.converters.core.ir import Field as SchemaField
from cdm_data_loaders.converters.core.ir import SchemaDocument, TypedNode
from cdm_data_loaders.converters.core.ir_values import Value, mutable_value
from cdm_data_loaders.utils.dlt.dlt_normalization import (
    DltTableNode,
    child_key,
    child_table_name,
    flatten_nodes,
    tables_parents_first,
)

logger = logging.getLogger(__name__)

SCHEMA_ENGINE_VERSION: Final = 11
TABLE_PATH_LENGTH: Final = 2
SCALARS: Final[dict[str, TDataType]] = {
    "null": "text",
    "boolean": "bool",
    "integer": "bigint",
    "object": "json",
    "array": "json",
    "map": "json",
    "any": "json",
}
FORMATS: Final[dict[str, TDataType]] = {
    "date-time": "timestamp",
    "datetime": "timestamp",
    "timestamp": "timestamp",
    "date": "date",
    "time": "time",
    "byte": "binary",
    "binary": "binary",
    "base64": "binary",
}
BINARY_FORMATS: Final = frozenset(name for name, kind in FORMATS.items() if kind == "binary")
LOGICAL: Final[dict[str, TDataType]] = {
    "decimal": "decimal",
    "wei": "wei",
    "timestamp": "timestamp",
    "timestamp-with-tz": "timestamp",
    "timestamp-without-tz": "timestamp",
    "date": "date",
    "time": "time",
    "binary": "binary",
    "fixed": "binary",
}
COLUMN_FACTS: Final = frozenset({"name", "data_type", "nullable", "description", "precision", "scale", "timezone"})


def _dlt_table_metadata(node: TypedNode) -> Mapping[str, Value] | None:
    """Return scoped table hints, excluding column provenance."""
    provenance = node.provenance
    if provenance is not None and provenance.source == "dlt" and len(provenance.path) == TABLE_PATH_LENGTH:
        return provenance.metadata
    return None


def _restore_table(tree: DltTableNode, node: TypedNode) -> DltTableNode:
    """Preserve source table hints around freshly derived structure."""
    metadata = _dlt_table_metadata(node)
    if metadata is None:
        return tree
    table = {key: mutable_value(value) for key, value in metadata.items() if key not in {"columns", "name", "parent"}}
    if "name" in metadata:
        table["name"] = tree.name
    if tree.parent is not None or "parent" in metadata:
        table["parent"] = tree.parent
    if tree.table["columns"] or "columns" in metadata:
        table["columns"] = tree.table["columns"]
    table.pop("description", None)
    if "description" in node.annotations:
        table["description"] = mutable_value(node.annotations["description"])
    return replace(tree, table=cast("TTableSchema", table))


def _column_hints(node: TypedNode) -> dict[str, Any]:
    """Translate logical precision, scale, widths and timezone without rounding."""
    precision = node.hints.precision
    if precision is None:
        precision = node.hints.bit_width or node.hints.length
    if precision is None and node.hints.logical_type in {"timestamp-with-tz", "timestamp-without-tz"}:
        precision = 6
    return {
        key: value
        for key, value in {"precision": precision, "scale": node.hints.scale, "timezone": node.hints.timezone}.items()
        if value is not None
    }


def _restore_column(column: TColumnSchema, node: TypedNode) -> TColumnSchema:
    """Copy compatible hints without resurrecting source types or removed fields."""
    metadata = node.extensions.get("x-dlt")
    if not isinstance(metadata, Mapping):
        return column
    restored = {key: mutable_value(value) for key, value in metadata.items() if key not in COLUMN_FACTS}
    restored.update(column)
    if "name" not in metadata:
        restored.pop("name", None)
    if "nullable" not in metadata and node.nullable is True:
        restored.pop("nullable", None)
    elif "nullable" in metadata and node.nullable is None:
        restored["nullable"] = None
    if node.type == "any" and node.hints.logical_type is None:
        restored.pop("data_type", None)
        if "data_type" in metadata:
            restored["data_type"] = None
    for key in ("precision", "scale", "timezone"):
        if key not in restored and key in metadata and metadata[key] is None:
            restored[key] = None
    return cast("TColumnSchema", restored)


def _is_json(node: TypedNode) -> bool:
    """Identify JSON-specific compatibility facts."""
    return node.provenance is not None and node.provenance.source == "json-schema"


def _boolean_schema(node: TypedNode) -> bool:
    """Distinguish JSON boolean schemas from empty schema objects."""
    return _is_json(node) and node.type in {"any", "never"} and not node.source_keywords


def _resolve_type(node: TypedNode) -> str:
    """Select a type without substituting branch structure or constraints."""
    if not _is_json(node):
        return node.type
    declaration = node.declared_type
    if isinstance(declaration, tuple):
        non_null = tuple(kind for kind in declaration if kind != "null")
        if len(non_null) != 1:
            if non_null:
                logger.warning("Collapsing multi-type union %r to 'text' (dlt has no union types).", declaration)
            return "text" if non_null else "null"
        declaration = non_null[0]
    if declaration is None:
        enum = node.constraints.get("enum")
        if isinstance(enum, tuple):
            return json_type_from_enum(list(enum))
        if node.inferred_type is not None:
            return node.inferred_type
    for name, branches in (("oneOf", node.one_of), ("anyOf", node.any_of)):
        if branches:
            logger.warning("Approximating '%s' by using only its first branch for type inference.", name)
            return "text" if _boolean_schema(branches[0]) else _resolve_type(branches[0])
    return declaration or "text"


def _single_item(node: TypedNode) -> TypedNode | None:
    """Select the first tuple member while preserving empty-prefix precedence."""
    items = node.prefix_items if node.prefix_items is not None else node.items
    item = (items[0] if items else None) if isinstance(items, tuple) else items
    if item is not None and _is_json(item) and not item.source_keywords and item.type != "any":
        return None
    return item


def _scalar_type(node: TypedNode) -> TDataType:
    """Map the selected type using the outer node's scalar constraints."""
    kind = _resolve_type(node)
    logical = node.hints.logical_type
    if not _is_json(node) and logical in LOGICAL:
        return LOGICAL[logical]
    if not _is_json(node) and kind in {"unknown", "never"}:
        msg = f"No dlt mapping for node type {kind!r} with logical type {logical!r}"
        raise ConversionError(msg)
    if kind == "number":
        return "decimal" if isinstance(node.constraints.get("multipleOf"), (int, float, Decimal)) else "double"
    if kind == "string":
        fmt = node.constraints.get("format")
        return FORMATS.get(fmt, "text") if isinstance(fmt, str) else "text"
    return SCALARS.get(kind, "text")


class DltEmitter(BaseModel):
    """Configure typed emission with JSON-to-dlt compatibility defaults."""

    model_config = ConfigDict(frozen=True)

    schema_name: str
    skip_nested_types: bool = False
    max_nesting: int = Field(default=10, ge=1)
    write_disposition: str | None = None
    flatten_scalars: bool = True

    @field_validator("schema_name")
    @classmethod
    def _validate_schema_name(cls, value: str) -> str:
        """Reject blank schema names."""
        if not value.strip():
            msg = "schema_name must be a non-empty string"
            raise ValueError(msg)
        return value

    def emit(self, document: SchemaDocument) -> TStoredSchema:
        """Emit one document as a stored schema with one root table."""
        root = document.root
        if not _is_json(root) and root.type != "object":
            msg = "A dlt table requires an object-compatible root"
            raise ConversionError(msg)
        native = root.provenance is not None and root.provenance.source == "dlt"
        name = document.name if native and document.name is not None else self.schema_name
        tree = self._table(root, name, None, 0, native=native)
        stored: TStoredSchema = {
            "name": self.schema_name,
            "version": 1,
            "previous_hashes": [],
            "engine_version": SCHEMA_ENGINE_VERSION,
            "tables": flatten_nodes({tree.name: tree}),
        }
        if native and document.provenance is not None and document.provenance.source == "dlt":
            envelope = document.provenance.metadata.get("envelope")
            if isinstance(envelope, Mapping):
                stored.update(
                    {key: mutable_value(value) for key, value in envelope.items() if key not in {"name", "tables"}}
                )
        return stored

    def to_schema(self, document: SchemaDocument) -> Schema:
        """Build a live schema, installing every parent before its children."""
        stored = self.emit(document)
        schema = Schema(stored["name"])
        for name, table in stored["tables"].items():
            table["name"] = name
            for column_name, column in table.setdefault("columns", {}).items():
                column["name"] = column_name
        for table in tables_parents_first(stored["tables"]):
            schema.update_table(table, normalize_identifiers=False)
        return schema

    def to_yaml(self, document: SchemaDocument) -> str:
        """Serialize the emitted stored schema without sorting table keys."""
        return yaml.safe_dump(self.emit(document), allow_unicode=True, default_flow_style=False, sort_keys=False)

    def _table(
        self, node: TypedNode, name: str, parent: str | None, depth: int, *, native: bool = False
    ) -> DltTableNode:
        """Recursively derive columns and child tables from object fields."""
        columns: dict[str, TColumnSchema] = {}
        children: dict[str, DltTableNode] = {}
        for prop in node.properties:
            child = self._child(prop, name, depth, native=native)
            if child is None:
                columns[prop.name] = self._column(prop, native=native)
            else:
                children[child.key] = child
        table = new_table(
            table_name=name,
            parent_table_name=parent,
            write_disposition=self.write_disposition if parent is None else None,
        )
        table["columns"] = columns
        if "description" in node.annotations:
            table["description"] = cast("str", mutable_value(node.annotations["description"]))
        tree = DltTableNode(name, parent, child_key(parent, name) if parent else name, table, children)
        if native:
            tree = _restore_table(tree, node)
        if parent is None and self.write_disposition is not None:
            tree.table["write_disposition"] = self.write_disposition
        return tree

    def _child(self, prop: SchemaField, parent: str, depth: int, *, native: bool = False) -> DltTableNode | None:
        """Flatten structures while retaining the legacy object-only depth limit."""
        node = prop.node
        if self.skip_nested_types:
            return None
        kind = _resolve_type(node)
        name = child_table_name(parent, prop.name)
        metadata = _dlt_table_metadata(node) if native else None
        if metadata is not None and node.provenance is not None:
            original = node.provenance.path[1]
            if (
                isinstance(original, str)
                and metadata.get("parent") == parent
                and child_key(parent, original) == prop.name
            ):
                name = original
        child = None
        if kind == "object":
            if depth < self.max_nesting:
                child = self._table(node, name, parent, depth + 1, native=native)
            else:
                logger.warning(
                    "Schema nesting limit (%d) reached at %r; ignoring object child.", self.max_nesting, prop.name
                )
        elif kind == "array" and (item := _single_item(node)) is not None:
            if _resolve_type(item) == "object":
                child = self._table(item, name, parent, depth + 1, native=native)
            elif self.flatten_scalars:
                value = (
                    new_column("value", _scalar_type(item))
                    if _is_json(node)
                    else self._column(SchemaField("value", item), native=native)
                )
                table = new_table(name, parent_table_name=parent)
                table["columns"] = {"value": value}
                child = DltTableNode(name, parent, child_key(parent, name), table)
        if child is not None:
            if native:
                child = _restore_table(child, node)
            if not _is_json(node) and "description" in prop.annotations:
                child.table["description"] = cast("str", mutable_value(prop.annotations["description"]))
        return child

    def _column(self, prop: SchemaField, *, native: bool = False) -> TColumnSchema:
        """Derive column types and JSON presence-based nullability."""
        node = prop.node
        if _boolean_schema(node):
            logger.warning("Boolean schema for property %r carries no type; using data_type 'text'.", prop.name)
        nullable = not prop.required if _is_json(node) or node.nullable is None else node.nullable
        column = new_column(prop.name, data_type=_scalar_type(node), nullable=nullable)
        if not _is_json(node):
            column.update(_column_hints(node))
        annotations = node.annotations if _is_json(node) else {**node.annotations, **prop.annotations}
        if "description" in annotations:
            column["description"] = cast("str", mutable_value(annotations["description"]))
        return _restore_column(column, node) if native else column
