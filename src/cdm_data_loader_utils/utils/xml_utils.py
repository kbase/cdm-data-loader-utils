"""
Shared XML helper utilities used by UniProt and UniRef parsers.

This module centralizes common operations:
- Safe text extraction
- Safe attribute extraction
- Property parsing
- Evidence / dbReference parsing
- Cleaning dictionaries
- Deduplicating lists
"""

import gzip
from collections.abc import Generator
from pathlib import Path
from typing import Any

from lxml.etree import Element, iterparse

from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger


def get_text(elem: Element | None, default: str | None = None) -> str | None:
    """Return elem.text if exists and non-empty."""
    if elem is None:
        return default
    if elem.text is None:
        return default
    text = elem.text.strip()
    return text or default


def get_attr(elem: Element | None, name: str, default: str | None = None) -> str | None:
    """Return elem.get(name) safely."""
    if elem is None:
        return default
    val = elem.get(name)
    return val.strip() if isinstance(val, str) else default


# ============================================================
# List / Node Finders
# ============================================================


def find_one(elem: Element, xpath: str, ns: dict[str, str]) -> Element | None:
    """Return first element matching xpath or None."""
    results = elem.findall(xpath, ns)
    return results[0] if results else None


def find_all_text(elem: Element, xpath: str, ns: dict[str, str]) -> list[str]:
    """Return list of text values from xpath matches (deduped)."""
    texts = []
    for node in elem.findall(xpath, ns):
        txt = get_text(node)
        if txt:
            texts.append(txt)
    return list(dict.fromkeys(texts))  # preserve order, dedupe


def safe_list(x) -> list[Any]:
    """Convert None → []."""
    if x is None:
        return []
    if isinstance(x, list):
        return x
    return [x]


# ============================================================
# dbReference / property parsing (shared by UniProt + UniRef)
# ============================================================


def parse_properties(dbref: Element | None, ns: dict[str, str]) -> dict[str, list[str]]:
    """
    Extract key/value pairs from <property type="..." value="..."> blocks.
    """
    if dbref is None:
        return {}
    props = {}
    for prop in dbref.findall("ns:property", ns):
        ptype = prop.attrib.get("type")
        pval = prop.attrib.get("value")
        if ptype and pval:
            if ptype not in props:
                props[ptype] = []
            props[ptype].append(pval)
    return props


def parse_db_references(elem: Element, ns: dict[str, str], pub_types=("PubMed", "DOI")):
    """Generic dbReference parser.

    - Identify publication IDs (PubMed, DOI)
    - Identify other cross-references (dbType:dbId)
    """
    publications = []
    others = []

    for dbref in elem.findall("ns:dbReference", ns):
        db_type = dbref.get("type")
        db_id = dbref.get("id")

        if not db_type or not db_id:
            continue

        if db_type in pub_types:
            publications.append(f"{db_type.upper()}:{db_id}")
        else:
            others.append(f"{db_type}:{db_id}")

    return publications, others


# ============================================================
# Dict Cleaning
# ============================================================


def clean_dict(d: dict[str, Any]) -> dict[str, Any]:
    """
    Remove keys whose value is None or empty list.
    """
    return {k: v for k, v in d.items() if v not in (None, [], {})}


def stream_xml_file(file_path: str | Path, element_with_ns: str) -> Generator[Element, Any]:
    """Stream XML elements from a file.

    :param file_path: path to the XML file; file can be gzipped or not.
    :type file_path: str | Path
    :param element_with_ns: name of the element (including namespace, in braces) to return; e.g. f"{{{UNIPROT_NS}}}entry"
    :type element_with_ns: str
    :yield: elements from the file
    :rtype: Generator[Element, Any]
    """
    logger = get_cdm_logger()
    if isinstance(file_path, Path):
        file_path = str(file_path)
    logger.info("Streaming XML from %s", file_path)
    open_fn = open
    if file_path.endswith(".gz"):
        open_fn = gzip.open

    with open_fn(file_path, "rb") as f:
        # with open(filepath, "rb") as f:
        for _, elem in iterparse(f, tag=(element_with_ns), remove_blank_text=True):
            logger.debug(elem)
            yield elem
            elem.clear()


def parse_head_matter(fh, ns_dict) -> dict[str, str]:
    """Parse the head matter of an XML file and return the data source information.

    :param fh: open file handle / data source
    :type fh:
    :raises ValueError: if the expected UniRef namespace is not found
    :raises ValueError:
    :return: data source information
    :rtype: dict[str, str]
    """
    ns = {}
    data_src = {}
    ctx = iterparse(fh, events=("start-ns", "end-ns", "start", "end"))
    for event, elem in ctx:
        if event == "start-ns":
            # output will be namespaces
            ns[elem[0]] = elem[1]
            continue

        if event == "end-ns":
            continue

        if event == "start" and elem.tag.endswith("entry"):
            # we're done here
            break
        elem.clear()

    if "" not in ns:
        msg = "No default namespace found!"
        raise ValueError(msg)

    if ns[""] != ns_dict["ns"]:
        if ns[""] in ns_dict.values():
            get_cdm_logger().warning("xmlns set to '%s'", ns[""])
            ns_dict["ns"] = ns[""]
        else:
            msg = f"Unexpected default namespace: got '{ns['']}', expected '{ns_dict['ns']}'"
            raise ValueError(msg)

    return data_src
