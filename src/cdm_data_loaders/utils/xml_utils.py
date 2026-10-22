"""Shared XML helper utilities used by XML readers and UniProt/UniRef parsers.

This module centralizes common operations:
- Safe text extraction
- Local (namespace-stripped) XML tag/attribute name extraction
"""

from lxml.etree import Element


def get_text(elem: Element | None, default: str | None = None) -> str | None:
    """Return elem.text if exists and non-empty."""
    if elem is None:
        return default
    if elem.text is None:
        return default
    text = elem.text.strip()
    return text or default


def local_name(tag: str) -> str:
    """Return the local (namespace-stripped) part of an XML tag or attribute name.

    Handles Clark notation (``"{http://example.com}entry"`` -> ``"entry"``) and returns
    already-unprefixed names unchanged.
    """
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) and "}" in tag else tag
