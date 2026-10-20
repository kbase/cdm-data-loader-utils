"""Shared XML helper utilities used by UniProt and UniRef parsers.

This module centralizes common operations:
- Safe text extraction
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
