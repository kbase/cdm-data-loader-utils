"""Tests for the xsd reader module.

Test XSDs live under tests/data/xsd and are loaded through the load_schema
fixture defined in conftest.py.
"""

from collections.abc import Callable
from pathlib import Path

import pytest
import xmlschema

from cdm_data_loaders.readers.xsd import find_list_and_single_child_paths


def test_find_list_and_single_child_paths_pass_no_children(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """A root element with a simple type has no children, so both result sets are empty."""
    schema = load_schema("no_children.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == set()
    assert single_paths == set()


def test_find_list_and_single_child_paths_pass_empty_schema(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """A schema with no global elements yields no paths at all."""
    schema = load_schema("empty_schema.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == set()
    assert single_paths == set()


def test_find_list_and_single_child_paths_pass_attribute_only_type_no_children(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """A complex type with only attributes and no child elements yields no paths."""
    schema = load_schema("attribute_only.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == set()
    assert single_paths == set()


@pytest.mark.parametrize(
    ("child_name", "expected_set_name"),
    [
        pytest.param("requiredSingle", "single", id="required-single-occurs-once"),
        pytest.param("optionalSingle", "single", id="optional-single-occurs-zero-or-one"),
        pytest.param("optionalRepeated", "list", id="optional-repeated-occurs-zero-or-more"),
        pytest.param("requiredRepeated", "list", id="required-repeated-occurs-one-or-more"),
        pytest.param("boundedRepeated", "list", id="bounded-repeated-max-occurs-two"),
    ],
)
def test_find_list_and_single_child_paths_pass_occurs_variety(
    load_schema: Callable[[str], xmlschema.XMLSchema],
    child_name: str,
    expected_set_name: str,
) -> None:
    """Children are sorted into list_paths or single_paths based on their maxOccurs."""
    schema = load_schema("occurs_variety.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    path = ("root", child_name)
    if expected_set_name == "list":
        assert path in list_paths
        assert path not in single_paths
    else:
        assert path in single_paths
        assert path not in list_paths


def test_find_list_and_single_child_paths_pass_occurs_variety_excludes_impossible_child(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """A child with maxOccurs=0 can never appear, so it is excluded from both result sets."""
    schema = load_schema("occurs_variety.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    path = ("root", "neverOccurs")
    assert path not in list_paths
    assert path not in single_paths


def test_find_list_and_single_child_paths_pass_choice_group_repeated(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """Children of a repeatable choice group are list paths even with their own maxOccurs=1."""
    schema = load_schema("nested_choice.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == {("root", "optionA"), ("root", "optionB")}
    assert single_paths == set()


def test_find_list_and_single_child_paths_pass_nested_sequence_flattened(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """A child nested two sequence levels deep is still reported against its direct parent element."""
    schema = load_schema("nested_sequence_deep.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == {("root", "deepChild")}
    assert single_paths == set()


def test_find_list_and_single_child_paths_pass_ref_element_resolved(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """A child declared with ref= resolves to the referenced global element's local name."""
    schema = load_schema("ref_element.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert ("root", "shared") in list_paths
    assert ("root", "shared") not in single_paths


def test_find_list_and_single_child_paths_pass_self_referential_type_no_infinite_recursion(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """A type that nests itself as its own child is recorded once and does not recurse forever."""
    schema = load_schema("self_referential.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == {("node", "child")}
    assert single_paths == set()


def test_find_list_and_single_child_paths_pass_mutual_recursion_no_infinite_recursion(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """Two types that reference each other are both traversed once without infinite recursion."""
    schema = load_schema("mutual_recursion.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == {("a", "b")}
    assert single_paths == {("b", "a")}


def test_find_list_and_single_child_paths_pass_uniref_like_schema(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """Test the function on a more realistic schema."""
    schema = load_schema("uniref_like.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths == {
        ("UniRef50", "entry"),
        ("UniRef90", "entry"),
        ("UniRef100", "entry"),
        ("entry", "property"),
        ("entry", "member"),
        ("dbReference", "property"),
    }
    assert single_paths == {
        ("entry", "representativeMember"),
        ("member", "dbReference"),
        ("member", "sequence"),
        ("representativeMember", "dbReference"),
        ("representativeMember", "sequence"),
    }


def test_find_list_and_single_child_paths_pass_returns_disjoint_sets(
    load_schema: Callable[[str], xmlschema.XMLSchema],
) -> None:
    """No path can be classified as both a list path and a single path at the same time."""
    schema = load_schema("uniref_like.xsd")

    list_paths, single_paths = find_list_and_single_child_paths(schema)

    assert list_paths.isdisjoint(single_paths)


def test_find_list_and_single_child_paths_fail_invalid_schema_type() -> None:
    """Passing an object without a schema-like elements mapping raises AttributeError."""
    not_a_schema = "this is a string, not an XMLSchema"

    with pytest.raises(AttributeError):
        find_list_and_single_child_paths(not_a_schema)  # type: ignore[arg-type]


def test_find_list_and_single_child_paths_fail_nonexistent_file_path(test_data_dir: Path) -> None:
    """Loading a schema from a path that does not exist raises an OSError from xmlschema."""
    missing_path = test_data_dir / "does_not_exist.xsd"

    with pytest.raises(OSError):
        xmlschema.XMLSchema(str(missing_path))
