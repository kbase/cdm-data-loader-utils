"""Tests for `parse_head_matter` from cdm_data_loaders.readers.xml."""

import gzip
import logging
from collections.abc import Callable
from pathlib import Path

import pytest
from lxml.etree import XMLSyntaxError

from cdm_data_loaders.readers.xml import parse_head_matter

# XML with no namespace declarations at all.
NO_NAMESPACE_XML = """<?xml version="1.0"?>
<library>
    <book id="1"><title>The Shining</title></book>
    <book id="2"><title>The Stand</title></book>
</library>
"""

# XML with a single default namespace.
DEFAULT_NS_URI = "http://example.com/lib"
DEFAULT_NAMESPACE_XML = f"""<?xml version="1.0"?>
<library xmlns="{DEFAULT_NS_URI}">
    <book id="1"><title>The Shining</title></book>
</library>
"""

# XML with a single prefixed namespace.
LIB_NS_URI = "http://example.com/lib"
PREFIXED_NAMESPACE_XML = f"""<?xml version="1.0"?>
<lib:library xmlns:lib="{LIB_NS_URI}">
    <lib:book id="1"><lib:title>The Shining</lib:title></lib:book>
</lib:library>
"""

# XML with a default namespace plus multiple prefixed namespaces.
XSI_NS_URI = "http://www.w3.org/2001/XMLSchema-instance"
DC_NS_URI = "http://purl.org/dc/elements/1.1/"
MULTI_NAMESPACE_XML = f"""<?xml version="1.0"?>
<library xmlns="{DEFAULT_NS_URI}" xmlns:xsi="{XSI_NS_URI}" xmlns:dc="{DC_NS_URI}">
    <book id="1"><dc:title>The Shining</dc:title></book>
</library>
"""

# A real-world-style UniProt document header, modelled on the actual UniProtKB XML.
UNIPROT_NS_URI = "https://uniprot.org/uniprot"
UNIPROT_XML = f"""<?xml version="1.0" encoding="UTF-8"?>
<uniprot xmlns="{UNIPROT_NS_URI}" xmlns:xsi="{XSI_NS_URI}">
    <entry dataset="Swiss-Prot"><accession>P12345</accession></entry>
</uniprot>
"""

# Namespaces declared on a child element rather than the root are NOT in scope on the
# root, so parsing (which stops at the root's start event) must not report them.
CHILD_ONLY_NAMESPACE_XML = f"""<?xml version="1.0"?>
<library>
    <book xmlns:dc="{DC_NS_URI}" id="1"><dc:title>The Shining</dc:title></book>
</library>
"""

MISMATCHED_TAG_XML = '<library><book id="1"><title>Broken</title></library>'


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
def test_parse_head_matter_pass_no_namespaces_returns_empty_dict(
    xml_path_factory: Callable[..., Path], as_path: bool
) -> None:
    """Verify a well-formed document with no namespace declarations returns an empty dict, for Path and str inputs."""
    path = xml_path_factory(NO_NAMESPACE_XML)
    file_path: Path | str = path if as_path else str(path)

    assert parse_head_matter(file_path) == {}


def test_parse_head_matter_pass_default_namespace_stored_under_empty_key(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify a default (unprefixed) namespace is stored under the empty-string key."""
    path = xml_path_factory(DEFAULT_NAMESPACE_XML)

    assert parse_head_matter(path) == {"": DEFAULT_NS_URI}


def test_parse_head_matter_pass_prefixed_namespace_stored_under_prefix(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify a single prefixed namespace is keyed by its prefix."""
    path = xml_path_factory(PREFIXED_NAMESPACE_XML)

    assert parse_head_matter(path) == {"lib": LIB_NS_URI}


def test_parse_head_matter_pass_multiple_namespaces_all_captured(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify a document with a default namespace plus multiple prefixed namespaces captures all of them."""
    path = xml_path_factory(MULTI_NAMESPACE_XML)

    assert parse_head_matter(path) == {
        "": DEFAULT_NS_URI,
        "xsi": XSI_NS_URI,
        "dc": DC_NS_URI,
    }


def test_parse_head_matter_pass_uniprot_style_header(xml_path_factory: Callable[..., Path]) -> None:
    """Verify a realistic UniProt-style header yields the expected default + xsi namespaces."""
    path = xml_path_factory(UNIPROT_XML)

    assert parse_head_matter(path) == {"": UNIPROT_NS_URI, "xsi": XSI_NS_URI}


def test_parse_head_matter_pass_child_only_namespace_not_reported(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify namespaces declared only on a child element are not reported (parsing stops at the root)."""
    path = xml_path_factory(CHILD_ONLY_NAMESPACE_XML)

    assert parse_head_matter(path) == {}


def test_parse_head_matter_pass_gzip_file_matches_plain_file(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify a gzip-compressed (.gz) file is transparently decompressed and parsed identically to a plain file."""
    plain_path = xml_path_factory(MULTI_NAMESPACE_XML, filename="plain.xml")
    gz_path = xml_path_factory(MULTI_NAMESPACE_XML, gzip_compress=True, filename="compressed.xml")

    plain_result = parse_head_matter(plain_path)
    gz_result = parse_head_matter(gz_path)

    assert gz_result == plain_result == {"": DEFAULT_NS_URI, "xsi": XSI_NS_URI, "dc": DC_NS_URI}


def test_parse_head_matter_pass_does_not_read_whole_document(
    xml_path_factory: Callable[..., Path], caplog: pytest.LogCaptureFixture
) -> None:
    """Verify parsing stops early: only namespaces in scope on the root are returned, even for large bodies."""
    body = "\n".join(f'    <book id="{i}"><title>Book {i}</title></book>' for i in range(10_000))
    big_xml = f'<?xml version="1.0"?>\n<library xmlns="{DEFAULT_NS_URI}">\n{body}\n</library>\n'
    path = xml_path_factory(big_xml)
    caplog.set_level(logging.DEBUG)

    result = parse_head_matter(path)

    assert result == {"": DEFAULT_NS_URI}
    assert any("Found 1 namespace declaration(s)" in message for message in caplog.messages)


def test_parse_head_matter_pass_logs_debug_messages(
    xml_path_factory: Callable[..., Path], caplog: pytest.LogCaptureFixture
) -> None:
    """Verify the debug logs emitted by parse_head_matter."""
    path = xml_path_factory(MULTI_NAMESPACE_XML)
    caplog.set_level(logging.DEBUG)

    parse_head_matter(path)

    assert any(m.startswith("Parsing head matter from") for m in caplog.messages)
    assert any(m.startswith("Found 3 namespace declaration(s)") for m in caplog.messages)
    assert {r.levelno for r in caplog.records} == {logging.DEBUG}


def test_parse_head_matter_fail_missing_file_raises_file_not_found_error(tmp_path: Path) -> None:
    """Verify a path to a non-existent file raises FileNotFoundError."""
    missing_path = tmp_path / "does_not_exist.xml"

    with pytest.raises(FileNotFoundError):
        parse_head_matter(missing_path)


def test_parse_head_matter_fail_missing_gzip_file_raises_file_not_found_error(tmp_path: Path) -> None:
    """Verify a path ending in .gz that doesn't exist raises FileNotFoundError rather than a gzip-specific error."""
    missing_path = tmp_path / "does_not_exist.xml.gz"

    with pytest.raises(FileNotFoundError):
        parse_head_matter(missing_path)


def test_parse_head_matter_fail_empty_file_raises_syntax_error(xml_path_factory: Callable[..., Path]) -> None:
    """Verify a zero-byte XML file raises lxml's XMLSyntaxError rather than silently returning an empty dict."""
    path = xml_path_factory("")

    with pytest.raises(XMLSyntaxError):
        parse_head_matter(path)


def test_parse_head_matter_fail_malformed_root_raises_syntax_error(xml_path_factory: Callable[..., Path]) -> None:
    """Verify content that is malformed before/at the root element raises lxml's XMLSyntaxError."""
    path = xml_path_factory("not xml at all <<<")

    with pytest.raises(XMLSyntaxError):
        parse_head_matter(path)


def test_parse_head_matter_pass_body_malformation_after_root_is_tolerated(
    xml_path_factory: Callable[..., Path],
) -> None:
    """Verify malformation deep in the body is NOT surfaced, since parsing stops at the root's start event.

    This documents the intended streaming behaviour: only the header (root element and its
    namespace declarations) is read, so a document whose sole malformation lies after the
    root has opened still yields the root's namespaces without error.
    """
    path = xml_path_factory(f'<library xmlns="{DEFAULT_NS_URI}"><book></wrong></library>')

    assert parse_head_matter(path) == {"": DEFAULT_NS_URI}


def test_parse_head_matter_fail_non_gzip_content_with_gz_extension_raises_error(tmp_path: Path) -> None:
    """Verify a file named `*.gz` but containing plain (non-gzip-compressed) text raises an error when read."""
    fake_gz_path = tmp_path / "not_really_gzipped.xml.gz"
    fake_gz_path.write_text(DEFAULT_NAMESPACE_XML, encoding="utf-8")

    with pytest.raises(OSError, match="Not a gzipped file"):
        parse_head_matter(fake_gz_path)
