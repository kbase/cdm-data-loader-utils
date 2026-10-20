"""Lossless, parent-linked normalization of dlt table definitions."""

from collections import deque
from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from typing import Final, cast

from dlt.common.schema.typing import TTableSchema

from cdm_data_loaders.converters.core.errors import ConversionError

NESTED_TABLE_SEPARATOR: Final[str] = "__"


class DltNormalizationError(ConversionError):
    """Raised for malformed tables, invalid parent graphs, or inconsistent trees."""


@dataclass(frozen=True, slots=True)
class DltTableNode:
    """A table with its mapping name, declared parent, property key, and raw definition.

    Roots use their mapping name as their key. Children are keyed by property name
    in declaration order. Table metadata and column definitions remain unfiltered.
    Identity fields are frozen; table and children mappings are independently owned.
    """

    name: str
    parent: str | None
    key: str
    table: TTableSchema
    children: dict[str, "DltTableNode"] = field(default_factory=dict)


def child_key(parent_name: str, child_name: str) -> str:
    """Remove exactly one declared-parent prefix, or return the full child name."""
    return child_name.removeprefix(f"{parent_name}{NESTED_TABLE_SEPARATOR}")


def child_table_name(parent_name: str, key: str) -> str:
    """Join a parent table name and property key without interpreting either path."""
    return f"{parent_name}{NESTED_TABLE_SEPARATOR}{key}"


def _copy_table(name: str, table: TTableSchema) -> TTableSchema:
    """Validate structural fields and copy all raw metadata without interpretation."""
    if not isinstance(name, str) or not name:
        message = "Table names must be nonempty strings."
        raise DltNormalizationError(message)
    if not isinstance(table, Mapping):
        message = f"Table {name!r} must be a mapping."
        raise DltNormalizationError(message)
    parent = table.get("parent")
    if parent is not None and (not isinstance(parent, str) or not parent):
        message = f"Table {name!r} parent must be a nonempty string or None."
        raise DltNormalizationError(message)
    columns = table.get("columns", {})
    if not isinstance(columns, Mapping):
        message = f"Table {name!r} columns must be a mapping."
        raise DltNormalizationError(message)
    for column_name, column in columns.items():
        if not isinstance(column_name, str) or not column_name:
            message = f"Table {name!r} column names must be nonempty strings."
            raise DltNormalizationError(message)
        if not isinstance(column, Mapping):
            message = f"Column {column_name!r} in table {name!r} must be a mapping."
            raise DltNormalizationError(message)
    return cast("TTableSchema", deepcopy(dict(table)))


def unflatten_tables(tables: Mapping[str, TTableSchema]) -> dict[str, DltTableNode]:
    """Build a forest using explicit parent links, preserving root and sibling order.

    Mapping keys identify tables, even when raw ``name`` metadata differs. Missing
    columns stay missing. Empty input yields an empty forest. Dangling parents,
    cycles (including disconnected cycles), and ambiguous child keys are errors.
    """
    if not isinstance(tables, Mapping):
        message = "Tables must be a mapping."
        raise DltNormalizationError(message)
    nodes: dict[str, DltTableNode] = {}
    roots: dict[str, DltTableNode] = {}
    for name, raw_table in tables.items():
        table = _copy_table(name, raw_table)
        parent = table.get("parent")
        key = child_key(parent, name) if parent is not None else name
        nodes[name] = DltTableNode(name=name, parent=parent, key=key, table=table)

    for name, node in nodes.items():
        if node.parent is None:
            roots[name] = node
            continue
        if node.parent not in nodes:
            message = f"Table {name!r} references parent {node.parent!r}, which does not exist."
            raise DltNormalizationError(message)
        siblings = nodes[node.parent].children
        if node.key in siblings:
            message = f"Tables under parent {node.parent!r} share child key {node.key!r}."
            raise DltNormalizationError(message)
        siblings[node.key] = node

    pending = deque(roots.values())
    visited: set[str] = set()
    while pending:
        node = pending.popleft()
        visited.add(node.name)
        pending.extend(node.children.values())
    if len(visited) != len(nodes):
        remaining = [name for name in nodes if name not in visited]
        prefix = "No root tables found: " if not roots else ""
        message = f"{prefix}Cycle in table parent links; unreachable tables: {remaining!r}."
        raise DltNormalizationError(message)
    return roots


def flatten_nodes(roots: Mapping[str, DltTableNode]) -> dict[str, TTableSchema]:
    """Copy a forest into tables in depth-first, parents-first declaration order.

    Raw definitions are preserved, including absent columns or parent fields.
    Node identity, mapping keys, and raw parent links must agree with the tree.
    Duplicate table names and repeated nodes (including cycles) are rejected.
    """
    tables: dict[str, TTableSchema] = {}
    pending: list[tuple[str, DltTableNode, str | None]] = [
        (key, node, None) for key, node in reversed(list(roots.items()))
    ]
    while pending:
        key, node, parent = pending.pop()
        if node.name in tables:
            message = f"Cycle or duplicate table name {node.name!r} in node tree."
            raise DltNormalizationError(message)
        table = _copy_table(node.name, node.table)
        expected_key = child_key(parent, node.name) if parent is not None else node.name
        if node.parent != parent or table.get("parent") != parent:
            message = f"Table {node.name!r} parent does not match its tree position."
            raise DltNormalizationError(message)
        if key != node.key or key != expected_key:
            message = f"Table {node.name!r} key does not match its tree position."
            raise DltNormalizationError(message)
        tables[node.name] = table
        pending.extend((child_name, child, node.name) for child_name, child in reversed(node.children.items()))
    return tables


def tables_parents_first(tables: Mapping[str, TTableSchema]) -> list[TTableSchema]:
    """Return isolated table definitions with each parent before its descendants."""
    return list(flatten_nodes(unflatten_tables(tables)).values())
