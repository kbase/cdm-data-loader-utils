"""Parser for UniProt metalink XML files.

These metadata files provide information and links for UniProt and related downloads.


"""

import datetime
from pathlib import Path
from typing import Any
from xml.etree.ElementTree import Element

from defusedxml.ElementTree import parse

from cdm_data_loader_utils.parsers.uniprot.xml_utils import get_text
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger, log_and_die

NS = {"ns": "http://www.metalinker.org/"}
NOW = datetime.datetime.now(tz=datetime.UTC)
COLUMNS = ["id", "db", "xref"]

logger = get_cdm_logger()


def parse_metalink(metalink_xml_path: Path | str) -> Element | None:
    """Parse the metalink file and return the root node."""
    document = parse(str(metalink_xml_path))
    root = document.getroot()
    if root is not None:
        return root

    return log_and_die("Could not find root for metalink file", RuntimeError)


def generate_data_source_table(metalink_xml_path: Path | str) -> dict[str, Any]:
    """Generate the data source information for the ID Mapping data."""
    root = parse_metalink(metalink_xml_path)
    if root is None:
        return {}

    data_source = {
        "license": get_text(root.find("./ns:license/ns:name", NS)),
        "publisher": get_text(root.find("./ns:publisher/ns:name", NS)),
        "resource_type": "dataset",
        "version": get_text(root.find("./ns:version", NS)),
    }
    missing = [k for k in data_source if not data_source[k]]
    if missing:
        log_and_die(
            f"Missing required elements from metalink file: {', '.join(missing)}",
            RuntimeError,
        )

    return data_source


def get_files(metalink_xml_path: Path | str, files_to_find: list[str] | None = None) -> dict[str, Any]:
    """Generate the data source information for the ID Mapping data."""
    root = parse_metalink(metalink_xml_path)
    if root is None:
        return {}

    if files_to_find is not None and files_to_find == []:
        logger.warning("Empty file list supplied to get_files: aborting.")
        return {}

    files = {}
    for f in root.findall("./ns:files/ns:file", NS):
        # get the name, size, any verification info
        name = f.get("name")
        # skip now if the file is not of interest
        if files_to_find and name not in files_to_find:
            continue

        size = get_text(f.find("./ns:size", NS))
        checksum = f.find("./ns:verification/ns:hash", NS)
        if checksum is not None:
            checksum_fn = checksum.get("type")
            checksum_value = get_text(checksum)
        else:
            checksum_fn = checksum_value = None
        dl_url = get_text(f.find("./ns:resources/ns:url[@location='us']", NS))
        files[name] = {
            "name": name,
            "size": size,
            "checksum": checksum_value,
            "checksum_fn": checksum_fn,
            "url": dl_url,
        }

    # report on unfound files
    if files_to_find:
        not_found = {f for f in files_to_find if f not in files}
        if not_found:
            msg = "The following files were not found: " + ", ".join(not_found)
            logger.warning(msg)

    return files
