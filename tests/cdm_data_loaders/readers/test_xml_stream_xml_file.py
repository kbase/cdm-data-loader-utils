"""Tests for `stream_xml_file` from cdm_data_loaders.readers.xml."""

import gzip
import logging
from collections.abc import Callable
from pathlib import Path

import pytest
from lxml.etree import XMLSyntaxError

from cdm_data_loaders.readers.xml import stream_xml_file

SIMPLE_LIBRARY_XML = """<?xml version="1.0"?>
<library>
    <book id="1"><title>The Shining</title></book>
    <book id="2"><title>The Stand</title></book>
    <book id="3"><title>The Tommyknockers</title></book>
</library>
"""

MIXED_TAGS_XML = """<?xml version="1.0"?>
<catalog>
    <book id="1"><title>Cell</title></book>
    <magazine id="1"><title>National Geographic</title></magazine>
    <book id="2"><title>Cujo</title></book>
</catalog>
"""

NAMESPACE_URI = "http://example.com/lib"
NAMESPACED_SIMPLE_XML = SIMPLE_LIBRARY_XML.replace("<library>", f'<library xmlns="{NAMESPACE_URI}">')

NAMESPACED_MIXED_TAG_XML = MIXED_TAGS_XML.replace("<catalog>", f'<catalog xmlns="{NAMESPACE_URI}">')

MISMATCHED_TAG_XML = '<library><book id="1"><title>Broken</title></library>'

EXPECTED_SIMPLE_TAGS_IDS = [("book", "1"), ("book", "2"), ("book", "3")]
EXPECTED_MIXED_TAGS_IDS = [("book", "1"), ("book", "2")]


@pytest.fixture
def xml_path_factory(tmp_path: Path) -> Callable[..., Path]:
    """Return a factory that writes XML content to a real file, optionally gzip-compressed."""

    def _make(content: str, *, gzip_compress: bool = False, filename: str = "data.xml") -> Path:
        if gzip_compress:
            path = tmp_path / f"{filename}.gz"
            with gzip.open(path, "wb") as f:
                f.write(content.encode("utf-8"))
        else:
            path = tmp_path / filename
            path.write_text(content, encoding="utf-8")
        return path

    return _make


@pytest.mark.parametrize("as_path", [True, False])
def test_stream_xml_file_pass_yields_matching_elements_in_order(
    xml_path_factory: Callable[..., Path], as_path: bool
) -> None:
    """Verify elements matching the requested tag are yielded in document order, for both Path and str inputs."""
    path = xml_path_factory(SIMPLE_LIBRARY_XML)
    file_path: Path | str = path if as_path else str(path)

    tags_ids = [(elem.tag, elem.get("id")) for elem in stream_xml_file(file_path, "book")]

    assert tags_ids == EXPECTED_SIMPLE_TAGS_IDS


def test_stream_xml_file_pass_gzip_file_yields_same_elements_as_plain_file(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify a gzip-compressed (.gz) file is transparently decompressed and parsed identically to a plain file."""
    plain_path = xml_path_factory(SIMPLE_LIBRARY_XML, filename="plain.xml")
    gz_path = xml_path_factory(SIMPLE_LIBRARY_XML, gzip_compress=True, filename="compressed.xml")

    plain_ids = [(elem.tag, elem.get("id")) for elem in stream_xml_file(plain_path, "book")]
    gz_ids = [(elem.tag, elem.get("id")) for elem in stream_xml_file(gz_path, "book")]

    assert gz_ids == plain_ids == EXPECTED_SIMPLE_TAGS_IDS


def test_stream_xml_file_pass_filters_out_non_matching_tags(xml_path_factory: Callable[..., Path]) -> None:
    """Verify only elements matching the requested tag are yielded when multiple tag types are interleaved."""
    path = xml_path_factory(MIXED_TAGS_XML)

    tags_seen = [(elem.tag, elem.get("id")) for elem in stream_xml_file(path, "book")]

    assert tags_seen == EXPECTED_MIXED_TAGS_IDS


def test_stream_xml_file_pass_no_matching_tag_yields_nothing(xml_path_factory: Callable[..., Path]) -> None:
    """Verify requesting a tag absent from a well-formed file yields no elements, rather than erroring."""
    path = xml_path_factory(SIMPLE_LIBRARY_XML)

    results = list(stream_xml_file(path, "chapter"))

    assert results == []


@pytest.mark.parametrize(
    ("tag", "expected"),
    [
        ("library", [("library", None)]),
        ("title", [("title", "The Shining"), ("title", "The Stand"), ("title", "The Tommyknockers")]),
    ],
)
def test_stream_xml_file_pass_other_tags(xml_path_factory: Callable[..., Path], tag: str, expected: list) -> None:
    """Verify that other elements in the XML structure can also be found by stream_xml_file."""
    path = xml_path_factory(SIMPLE_LIBRARY_XML)
    results = [(elem.tag, elem.text) for elem in stream_xml_file(path, tag)]
    assert results == expected


@pytest.mark.parametrize(
    ("input_xml", "expected"),
    [(NAMESPACED_SIMPLE_XML, EXPECTED_SIMPLE_TAGS_IDS), (NAMESPACED_MIXED_TAG_XML, EXPECTED_MIXED_TAGS_IDS)],
)
def test_stream_xml_file_pass_namespaced_tag_matches_default_namespace_elements(
    xml_path_factory: Callable[..., Path], input_xml: str, expected: list[tuple[str, str]]
) -> None:
    """Verify a brace-qualified tag (e.g. '{ns}book') correctly matches elements under a default XML namespace."""
    file_path = xml_path_factory(input_xml)
    namespaced_tag = f"{{{NAMESPACE_URI}}}book"

    tags_ids = [(elem.tag, elem.get("id")) for elem in stream_xml_file(file_path, namespaced_tag)]
    assert tags_ids == [(namespaced_tag, el[1]) for el in expected]


def test_stream_xml_file_pass_bare_tag_does_not_match_namespaced_elements(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify a plain (non-namespace-qualified) tag does NOT match elements declared under an XML namespace."""
    path = xml_path_factory(NAMESPACED_SIMPLE_XML)

    results = list(stream_xml_file(path, "book"))

    assert results == []


def test_stream_xml_file_pass_blank_text_between_child_elements_is_removed(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify indentation whitespace between child elements is stripped due to `remove_blank_text=True`."""
    path = xml_path_factory(SIMPLE_LIBRARY_XML)

    book = next(stream_xml_file(path, "book"))

    assert book.text is None


def test_stream_xml_file_pass_element_content_intact_before_clear(xml_path_factory: Callable[..., Path]) -> None:
    """Verify a yielded element still has its title/attributes intact at the moment it's received by the caller."""
    path = xml_path_factory(SIMPLE_LIBRARY_XML)
    gen = stream_xml_file(path, "book")

    first = next(gen)

    assert first.get("id") == "1"
    assert first.find("title").text == "The Shining"


def test_stream_xml_file_pass_previously_yielded_elements_are_cleared_after_advancing(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify each element is cleared (text/children removed) once the generator advances past it, confirming streaming memory behavior."""
    path = xml_path_factory(SIMPLE_LIBRARY_XML)

    all_elements = list(stream_xml_file(path, "book"))

    assert all(elem.text is None for elem in all_elements)
    assert all(len(elem) == 0 for elem in all_elements)
    assert all(elem.find("title") is None for elem in all_elements)


def test_stream_xml_file_pass_logs_debug_message_per_yielded_element(
    xml_path_factory: Callable[..., Path], caplog: pytest.LogCaptureFixture
) -> None:
    """Verify the debug logs emitted by stream_xml_file."""
    path = xml_path_factory(SIMPLE_LIBRARY_XML)
    caplog.set_level(logging.DEBUG)
    n_expected_records = 4

    list(stream_xml_file(path, "book"))

    assert len(caplog.records) == n_expected_records
    assert caplog.records[0].message.startswith("Streaming XML from")
    log_levels = {r.levelno for r in caplog.records}
    assert log_levels == {logging.DEBUG}
    for rec in caplog.records[1:]:
        assert rec.message.startswith("<Element book at")


def test_stream_xml_file_fail_missing_file_raises_file_not_found_error(tmp_path: Path) -> None:
    """Verify a path to a non-existent file raises FileNotFoundError."""
    missing_path = tmp_path / "does_not_exist.xml"

    with pytest.raises(FileNotFoundError):
        list(stream_xml_file(missing_path, "book"))


def test_stream_xml_file_fail_missing_gzip_file_raises_file_not_found_error(tmp_path: Path) -> None:
    """Verify a path ending in .gz that doesn't exist raises FileNotFoundError rather than a gzip-specific error."""
    missing_path = tmp_path / "does_not_exist.xml.gz"

    with pytest.raises(FileNotFoundError):
        list(stream_xml_file(missing_path, "book"))


def test_stream_xml_file_fail_malformed_xml_raises_syntax_error(xml_path_factory: Callable[..., Path]) -> None:
    """Verify malformed XML (mismatched closing tags) raises lxml's XMLSyntaxError."""
    path = xml_path_factory(MISMATCHED_TAG_XML)

    with pytest.raises(XMLSyntaxError):
        list(stream_xml_file(path, "book"))


def test_stream_xml_file_fail_empty_file_raises_syntax_error(xml_path_factory: Callable[..., Path]) -> None:
    """Verify a zero-byte XML file raises lxml's XMLSyntaxError rather than silently yielding nothing."""
    path = xml_path_factory("")

    with pytest.raises(XMLSyntaxError):
        list(stream_xml_file(path, "book"))


def test_stream_xml_file_fail_non_gzip_content_with_gz_extension_raises_error(tmp_path: Path) -> None:
    """Verify a file named `*.gz` but containing plain (non-gzip-compressed) text raises an error when read."""
    fake_gz_path = tmp_path / "not_really_gzipped.xml.gz"
    fake_gz_path.write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    with pytest.raises(OSError, match="Not a gzipped file"):
        list(stream_xml_file(fake_gz_path, "book"))
