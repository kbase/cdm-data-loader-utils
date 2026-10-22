"""Tests for the xml_utils module."""

import pytest
from lxml.etree import Element

from cdm_data_loaders.utils.xml_utils import get_text, local_name


def make_element_with_text(raw_text: str | None) -> Element:
    """Create an element whose text is set to raw_text."""
    elem = Element("tag")
    elem.text = raw_text
    return elem


@pytest.mark.parametrize("default", [None, "n/a"], ids=["default_none", "default_custom"])
@pytest.mark.parametrize("raw_text", ["hello", "  hello  "], ids=["no_padding", "padded"])
def test_get_text_pass_returns_stripped_text(raw_text: str, default: str | None) -> None:
    """Non-empty element text is returned stripped regardless of the default."""
    assert get_text(make_element_with_text(raw_text), default) == "hello"


@pytest.mark.parametrize("default", [None, "n/a"], ids=["default_none", "default_custom"])
@pytest.mark.parametrize(
    "missing_text",
    [None, "", "   "],
    ids=["elem_none", "text_none", "text_whitespace"],
)
def test_get_text_pass_returns_default_when_text_missing(missing_text: str | None, default: str | None) -> None:
    """A None element or missing/whitespace-only text yields the default value."""
    elem = None if missing_text is None else make_element_with_text(missing_text)
    assert get_text(elem, default) == default


def test_get_text_fail_rejects_non_element() -> None:
    """A non-Element input raises AttributeError."""
    with pytest.raises(AttributeError):
        get_text(123)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("tag", "expected"),
    [
        ("{http://uniprot.org/uniprot}entry", "entry"),
        ("{}entry", "entry"),
        ("{http://example.org/{nested}}item", "item"),
        ("entry", "entry"),
        ("", ""),
    ],
    ids=[
        "clark_notation",
        "empty_namespace_braces",
        "nested_braces_in_namespace",
        "unprefixed_tag",
        "empty_string",
    ],
)
def test_local_name_pass_extracts_local_part(tag: str, expected: str) -> None:
    """local_name strips Clark-notation namespaces and preserves unprefixed names."""
    assert local_name(tag) == expected
