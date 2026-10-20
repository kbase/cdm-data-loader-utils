"""Pure table graph normalization tests."""

import json
from collections import UserDict
from copy import deepcopy
from dataclasses import FrozenInstanceError
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import pytest

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.utils.dlt.dlt_normalization import (
    DltNormalizationError,
    DltTableNode,
    child_key,
    child_table_name,
    flatten_nodes,
    tables_parents_first,
    unflatten_tables,
)

if TYPE_CHECKING:
    from dlt.common.schema.typing import TTableSchema


@pytest.mark.parametrize(
    ("parent", "name", "key"),
    [
        ("entry", "entry__property", "property"),
        ("entry", "entry__entry__property", "entry__property"),
        ("entry__member", "entry__member__db_reference", "db_reference"),
        ("entry", "other__property", "other__property"),
        ("entry", "entry", "entry"),
        ("entry", "entry__", ""),
        ("", "__value", "value"),
    ],
    ids=[
        "direct-child",
        "one-prefix-only",
        "nested-parent",
        "unrelated-name",
        "same-name",
        "empty-key",
        "empty-parent",
    ],
)
def test_child_key_pass_declared_parent(parent: str, name: str, key: str) -> None:
    """Strip one exact declared-parent prefix without XML decoding."""
    assert child_key(parent, name) == key
    assert child_key(parent, child_table_name(parent, key)) == key


@pytest.mark.parametrize(
    ("parent", "key", "expected"),
    [("entry", "db_reference", "entry__db_reference"), ("root", "root__tags", "root__root__tags"), ("", "", "__")],
    ids=["simple-key", "literal-path-key", "empty-fragments"],
)
def test_child_table_name_pass_literal_join(parent: str, key: str, expected: str) -> None:
    """Join names without decoding, stripping, or normalizing fragments."""
    assert child_table_name(parent, key) == expected


def test_unflatten_tables_pass_precise_forest_roundtrip() -> None:
    """Preserve every raw field while following explicit parents in declaration order."""
    tables = json.loads(Path("tests/data/dlt_normalization/tables.json").read_text(encoding="utf-8"))
    original = deepcopy(tables)
    roots = unflatten_tables(tables)
    leaf = DltTableNode("unrelated__leaf", "unrelated", "leaf", tables["unrelated__leaf"])
    unrelated = DltTableNode("unrelated", "root", "unrelated", tables["unrelated"], {"leaf": leaf})
    tags = DltTableNode("root__root__tags", "root", "root__tags", tables["root__root__tags"])
    assert roots == {
        "catalog__literal": DltTableNode("catalog__literal", None, "catalog__literal", tables["catalog__literal"]),
        "root": DltTableNode("root", None, "root", tables["root"], {"root__tags": tags, "unrelated": unrelated}),
        "root__looks_nested": DltTableNode(
            "root__looks_nested", None, "root__looks_nested", tables["root__looks_nested"]
        ),
    }
    assert list(roots) == ["catalog__literal", "root", "root__looks_nested"]
    assert list(roots["root"].children) == ["root__tags", "unrelated"]
    flattened = flatten_nodes(roots)
    assert flattened == original
    expected_order = [
        "catalog__literal",
        "root",
        "root__root__tags",
        "unrelated",
        "unrelated__leaf",
        "root__looks_nested",
    ]
    assert list(flattened) == expected_order
    assert tables_parents_first(tables) == [original[name] for name in expected_order]
    assert unflatten_tables(flattened) == roots
    assert flatten_nodes(unflatten_tables(tables)) == flattened
    assert tables == original


def test_unflatten_tables_pass_empty() -> None:
    """Empty forests and table mappings round-trip without introducing tables."""
    assert unflatten_tables({}) == {}
    assert flatten_nodes({}) == {}
    assert tables_parents_first({}) == []


def test_flatten_nodes_pass_mapping_input() -> None:
    """Accept ordered Mapping implementations without requiring reversible item views."""
    tables: dict[str, TTableSchema] = {"root": {"columns": {}}}
    roots = UserDict(unflatten_tables(UserDict(tables)))
    assert flatten_nodes(roots) == tables


@pytest.mark.parametrize(
    ("tables", "message"),
    [
        (None, "Tables must be a mapping"),
        ([], "Tables must be a mapping"),
        ({"": {}}, "Table names must be nonempty strings"),
        ({1: {}}, "Table names must be nonempty strings"),
        ({"root": None}, "must be a mapping"),
        ({"root": {"parent": ""}}, "parent must be"),
        ({"root": {"parent": []}}, "parent must be"),
        ({"root": {"parent": False}}, "parent must be"),
        ({"root": {"columns": None}}, "columns must be a mapping"),
        ({"root": {"columns": []}}, "columns must be a mapping"),
        ({"root": {"columns": {"": {}}}}, "column names must be nonempty strings"),
        ({"root": {"columns": {1: {}}}}, "column names must be nonempty strings"),
        ({"root": {"columns": {"bad": []}}}, "must be a mapping"),
        ({"root": {"parent": "missing"}}, "references parent 'missing'"),
        ({"root": {"parent": "root"}}, "No root tables found.*Cycle"),
        ({"first": {"parent": "second"}, "second": {"parent": "first"}}, "No root tables found.*Cycle"),
        ({"root": {}, "first": {"parent": "second"}, "second": {"parent": "first"}}, "Cycle"),
        ({"root": {}, "first": {"parent": "second"}, "second": {"parent": "second"}}, "Cycle"),
        ({"root": {}, "root__child": {"parent": "root"}, "child": {"parent": "root"}}, "share child key"),
    ],
    ids=[
        "null-tables",
        "list-tables",
        "empty-name",
        "numeric-name",
        "null-table",
        "empty-parent",
        "list-parent",
        "boolean-parent",
        "null-columns",
        "list-columns",
        "empty-column-name",
        "numeric-column-name",
        "list-column",
        "dangling-parent",
        "self-cycle",
        "rootless-cycle",
        "disconnected-cycle",
        "descendant-of-disconnected-cycle",
        "child-key-collision",
    ],
)
def test_unflatten_tables_fail_invalid_graph(tables: Any, message: str) -> None:  # noqa: ANN401
    """Reject structural errors and all cyclic components without mutating input."""
    original = deepcopy(tables)
    with pytest.raises(DltNormalizationError, match=message) as raised:
        unflatten_tables(tables)
    assert isinstance(raised.value, ConversionError)
    assert tables == original
    with pytest.raises(DltNormalizationError, match=message):
        tables_parents_first(tables)


def test_unflatten_tables_pass_isolated_copies() -> None:
    """Input, trees, flattened output, and siblings do not share mutable metadata."""
    metadata = {"values": ["original"]}
    tables: dict[str, Any] = {
        "root": {"columns": {"value": {"x-hint": metadata}}, "x-hint": metadata},
        "child": {"parent": "root", "columns": {"value": {"x-hint": metadata}}, "x-hint": metadata},
    }
    original = deepcopy(tables)
    roots = unflatten_tables(tables)
    second_tree = unflatten_tables(tables)
    flattened = cast("dict[str, Any]", flatten_nodes(roots))
    second_flattened = flatten_nodes(roots)
    flattened["child"]["columns"]["value"]["x-hint"]["values"].append("flat")
    assert tables == original
    assert roots == second_tree
    assert second_flattened == original
    root_table = cast("dict[str, Any]", roots["root"].table)
    root_table["x-hint"]["values"].append("tree")
    assert roots["root"].children["child"].table == original["child"]
    assert tables == original
    assert flatten_nodes(second_tree) == original
    metadata["values"].append("input")
    assert flatten_nodes(second_tree) == original
    with pytest.raises(FrozenInstanceError):
        roots["root"].name = "changed"


@pytest.mark.parametrize("depth", [1, 1200], ids=["single-root", "beyond-python-recursion-limit"])
def test_unflatten_tables_pass_deep_graph(depth: int) -> None:
    """Normalize deep parent graphs iteratively and flatten them in topological order."""
    tables: dict[str, TTableSchema] = {"table_0": {"columns": {}}}
    for index in range(1, depth):
        tables[f"table_{index}"] = {"parent": f"table_{index - 1}", "columns": {}}
    reversed_tables = dict(reversed(tables.items()))
    flattened = flatten_nodes(unflatten_tables(reversed_tables))
    assert flattened == tables
    assert list(flattened) == list(tables)


@pytest.mark.parametrize(
    ("node", "mapping_key", "message"),
    [
        (DltTableNode("root", None, "root", {"parent": "other"}), "root", "parent does not match"),
        (DltTableNode("root", "other", "root", {}), "root", "parent does not match"),
        (DltTableNode("root", None, "wrong", {}), "root", "key does not match"),
        (DltTableNode("root", None, "root", {}), "wrong", "key does not match"),
    ],
    ids=["raw-parent-mismatch", "node-parent-mismatch", "node-key-mismatch", "mapping-key-mismatch"],
)
def test_flatten_nodes_fail_inconsistent_tree(node: DltTableNode, mapping_key: str, message: str) -> None:
    """Reject conflicting identity fields instead of silently rewriting raw metadata."""
    with pytest.raises(DltNormalizationError, match=message):
        flatten_nodes({mapping_key: node})


@pytest.mark.parametrize("reuse_node", [True, False], ids=["cyclic-node", "duplicate-name"])
def test_flatten_nodes_fail_repeated_table(reuse_node: bool) -> None:
    """Reject cycles and duplicate names in manually constructed trees."""
    root = DltTableNode("root", None, "root", {})
    root.children["root"] = root if reuse_node else DltTableNode("root", "root", "root", {"parent": "root"})
    with pytest.raises(DltNormalizationError, match="Cycle or duplicate table name"):
        flatten_nodes({"root": root})
