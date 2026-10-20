"""dlt table normalization and source-preserving reader policy."""

from typing import Any

import pytest

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionError
from cdm_data_loaders.converters.readers.dlt import DltReader


def test_read_pass_forest_filtering_and_scalar_children() -> None:
    """Roots retain order and child classification uses filtered columns."""
    tables = {
        "other": {"columns": {}},
        "root__tags": {
            "parent": "root",
            "description": None,
            "columns": {
                "_dlt_id": {"data_type": "text", "nullable": False},
                "value": {"data_type": "text", "nullable": False},
                "value__v_int": {"data_type": "bigint"},
            },
        },
        "root": {"columns": {"id": {"data_type": "bigint", "precision": 32, "nullable": False}}},
        "root__obj": {"parent": "root", "columns": {"data": {"data_type": "json", "nullable": True}}},
    }
    documents = DltReader().read({"name": "stored", "tables": tables})
    assert tuple(documents) == ("other", "root")
    root = documents["root"].root
    assert tuple(prop.name for prop in root.properties) == ("id", "tags", "obj")
    assert tuple(prop.required for prop in root.properties) == (True, True, False)
    assert (root.properties[0].node.hints.bit_width,) == (32,)
    tags = root.properties[1].node
    assert tags.type == "array"
    assert tags.items.type == "string"
    assert tags.annotations == {"description": None}
    assert root.properties[2].node.properties[0].node.type == "any"
    assert root.properties[2].node.properties[0].node.nullable is True
    assert documents["root"].provenance.metadata["envelope"] == {"name": "stored"}
    tables["root"]["columns"].clear()
    assert "id" in root.provenance.metadata["columns"]
    included = DltReader(include_dlt_columns=True, include_variant_columns=True).read(tables)
    assert included["root"].root.properties[0].node.type == "object"


@pytest.mark.parametrize("mode", ["object", "array"], ids=["object-policy", "array-policy"])
def test_read_pass_child_required_propagation(mode: str) -> None:
    """Only direct visible required columns propagate child presence."""
    tables = {
        "root": {"columns": {}},
        "child": {"parent": "root", "columns": {"_dlt_id": {"nullable": False}}},
        "grandchild": {"parent": "child", "columns": {"id": {"data_type": "bigint", "nullable": False}}},
    }
    document = DltReader(child_table_mode=mode).read(tables)["root"]
    child = document.root.properties[0]
    assert child.required is False
    assert child.node.type == mode
    object_node = child.node if mode == "object" else child.node.items
    assert object_node.properties[0].required is True
    assert document.provenance.metadata["child_table_mode"] == mode


@pytest.mark.parametrize(
    ("data_type", "expected"),
    [
        ("text", "string"),
        ("bigint", "integer"),
        ("double", "number"),
        ("decimal", "number"),
        ("bool", "boolean"),
        ("timestamp", "string"),
        ("date", "string"),
        ("time", "string"),
        ("binary", "string"),
        ("json", "any"),
        ("wei", "integer"),
        (None, "any"),
    ],
    ids=[
        "text",
        "bigint",
        "double",
        "decimal",
        "bool",
        "timestamp",
        "date",
        "time",
        "binary",
        "json",
        "wei",
        "incomplete",
    ],
)
def test_read_pass_scalar_facts(data_type: str | None, expected: str) -> None:
    """Source logical types remain facts rather than target string fallbacks."""
    node = DltReader().read({"root": {"columns": {"value": {"data_type": data_type}}}})["root"].root.properties[0].node
    assert node.type == expected
    assert node.hints.logical_type == data_type
    assert node.hints.timezone is None
    assert node.nullable is True


def test_read_pass_custom_hints_and_unknown_type(caplog: pytest.LogCaptureFixture) -> None:
    """Unknown readable types retain their source name and log without coercion."""
    registry = DEFAULT_EXTENSIONS.extend_payload("x-dlt", {"x-custom": {"type": ["array", "null"]}})
    column = {"data_type": "future", "timezone": True, "x-custom": [None], "description": None}
    node = (
        DltReader(extension_registry=registry)
        .read({"root": {"columns": {"a": column}}})["root"]
        .root.properties[0]
        .node
    )
    assert node.type == "unknown"
    assert node.hints.timezone is True
    assert node.annotations == {"description": None}
    assert node.extensions["x-dlt"]["x-custom"] == (None,)
    assert "future" in caplog.text
    with pytest.raises(ExtensionError):
        DltReader().read({"root": {"columns": {"a": column}}})


@pytest.mark.parametrize(
    "source",
    [
        {},
        {"tables": {}},
        {"tables": []},
        {"orphan": {"parent": "missing", "columns": {}}},
        {"first": {"parent": "second"}, "second": {"parent": "first"}},
        {"root": {"columns": {}}, "cycle": {"parent": "cycle"}},
    ],
    ids=["empty", "empty-tables", "bad-tables", "orphan", "rootless-cycle", "disconnected-cycle"],
)
def test_read_fail_invalid_graph(source: dict[str, Any]) -> None:
    """Normalizer graph errors surface before schema reading."""
    with pytest.raises(ConversionError):
        DltReader().read(source)


def test_dlt_reader_fail_invalid_options() -> None:
    """Invalid reader policies fail immediately."""
    with pytest.raises(ValueError, match="child_table_mode"):
        DltReader(child_table_mode="guess")
    with pytest.raises(TypeError):
        DltReader(include_dlt_columns="yes")


def test_read_fail_non_mapping() -> None:
    """The public reader rejects non-mapping input."""
    with pytest.raises(ConversionError, match="mapping"):
        DltReader().read([])


def test_read_pass_child_collision_and_array_annotations() -> None:
    """Child replacements retain required ordering and array-only descriptions."""
    tables = {
        "root": {"columns": {"child": {"data_type": "text", "nullable": False}}},
        "root__child": {"parent": "root", "description": "child rows", "columns": {"optional": {"data_type": "text"}}},
    }
    root = DltReader(child_table_mode="array").read(tables)["root"].root
    assert root.required_names == ("child",)
    assert root.properties[0].required is True
    assert root.properties[0].node.annotations == {"description": "child rows"}
    assert root.properties[0].node.items.annotations == {}
