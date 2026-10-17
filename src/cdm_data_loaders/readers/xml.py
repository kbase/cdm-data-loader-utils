"""Common reusable XML pipeline elements for the xmltodict-based pipelines."""

from collections.abc import Callable, Generator, Iterable
from logging import Logger, getLogger
from pathlib import Path
from typing import Any

import dlt
import xmltodict
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.extract import DltResource
from dlt.extract.items import DataItemWithMeta
from dlt.sources.filesystem import filesystem
from lxml.etree import Element, iterparse, tostring

from cdm_data_loaders.core.fields import DEFAULT_XML_FILE_GLOB
from cdm_data_loaders.core.settings import CtsSettings
from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings
from cdm_data_loaders.utils.batcher import get_file_batches
from cdm_data_loaders.utils.buffer import DictBuffer, ListBuffer
from cdm_data_loaders.utils.gz import open_maybe_gzip

logger: Logger = getLogger(__name__)

DEFAULT_XMLTODICT_ARGS: dict[str, str | bool | int | None] = {"attr_prefix": "_"}


def parse_head_matter(file_path: str | Path) -> dict[str, str]:
    """Parse the namespace declarations from the head matter of an XML file.

    Reads the beginning of an XML document (gzipped or plain) and collects all XML
    namespace declarations (``xmlns`` / ``xmlns:prefix`` attributes) that are in scope
    on the root element. Parsing stops as soon as the root element is fully opened, so
    the whole document is never loaded into memory.

    The returned mapping uses namespace prefixes as keys and namespace URIs as values.
    The default (unprefixed) namespace is stored under the empty-string key ``""``.

    :param file_path: path to the XML file; the file can be gzipped or not.
    :type file_path: str | Path
    :raises FileNotFoundError: if the file does not exist.
    :raises lxml.etree.XMLSyntaxError: if the file is empty or malformed.
    :return: mapping of namespace prefix to namespace URI.
    :rtype: dict[str, str]
    """
    logger.debug("Parsing head matter from %s", file_path)

    namespaces: dict[str, str] = {}
    with open_maybe_gzip(file_path) as fh:
        ctx = iterparse(fh, events=("start-ns", "start"))
        for event, elem in ctx:
            if event == "start-ns":
                # elem is a (prefix, uri) tuple; the default namespace has an empty prefix.
                prefix, uri = elem
                namespaces[prefix] = uri
                continue

            # The first "start" event corresponds to the root element; by this point every
            # namespace declared on the root has already been emitted as a "start-ns" event.
            break

    logger.debug("Found %d namespace declaration(s) in %s", len(namespaces), file_path)
    return namespaces


def stream_xml_file(file_path: str | Path, element_with_ns: str) -> Generator[Element, Any]:
    """Stream XML elements from a file.

    :param file_path: path to the XML file; file can be gzipped or not.
    :type  file_path: str | Path
    :param element_with_ns: name of the element (including namespace, in braces) to return; e.g. f"{{{UNIPROT_NS}}}entry"
    :type  element_with_ns: str
    :yield: elements from the file
    :rtype: Generator[Element, Any]
    """
    logger.debug("Streaming XML from %s", file_path)

    with open_maybe_gzip(file_path) as f:
        for _, elem in iterparse(f, tag=(element_with_ns), remove_blank_text=True):
            logger.debug(elem)
            yield elem
            elem.clear()


def _log_interval_progress(n_entries: int, log_interval: int) -> None:
    """Log the processed-entry count every log_interval elements."""
    if (n_entries + 1) % log_interval == 0:
        logger.debug("Processed %d entries", n_entries + 1)


def _log_interval_final(n_entries: int, file_path: Path, log_interval: int) -> None:
    """Log the trailing entry count if the final checkpoint was not an interval boundary."""
    if n_entries >= 0 and (n_entries + 1) % log_interval != 0:
        logger.debug("Processed %d entries from %s", n_entries + 1, file_path.name)


def process_xml_file_to_dict(settings: XmlToDictSettings, file_path: Path) -> Generator[DataItemWithMeta]:
    """Convert XML to dictionary form.

    Note: each incoming element produces a single dictionary as output.

    Use process_xml_file for parsing functions that return data spread across several tables.

    :param settings: pipeline config, including xml_tag, table_name, and buffer_size
    :type settings: XmlToDictSettings
    :param file_path: file path, as a string
    :type file_path: Path
    :yield: pages of dictionary items
    :rtype: Generator[DataItemWithMeta]
    """
    logger.info("Reading from %s", str(file_path))
    n_entries = -1

    buffer = ListBuffer(table_name=settings.table_name, max_items=settings.buffer_size)
    for n_entries, element in enumerate(stream_xml_file(file_path, settings.xml_tag)):
        parsed_element = xmltodict.parse(tostring(element), **DEFAULT_XMLTODICT_ARGS, **settings.xmltodict_args)
        if parsed_element:
            for k, v in parsed_element.items():
                # remove the xmlns declarations
                parsed_element[k] = {kv: val for kv, val in v.items() if not kv.startswith("_xmlns")}
            yield from buffer.add_item(parsed_element)
        _log_interval_progress(n_entries, settings.log_interval)

    _log_interval_final(n_entries, file_path, settings.log_interval)
    yield from buffer.flush()


def process_xml_file(
    settings: CtsSettings,
    xml_tag: str,
    parse_fn: Callable,
    file_path: Path,
) -> Generator[DataItemWithMeta, Any]:
    """Core generator shared by XML-based dlt pipeline resources.

    This processor is expected to return a dictionary of table names and lists of rows, unlike the
    xmltodict parser, which creates a single dictionary for each element.

    :param settings: pipeline config with input_dir and buffer_size
    :type  settings: CtsSettings
    :param xml_tag: XML element tag to stream
    :type  xml_tag: str
    :param parse_fn: callable(element, timestamp, file_path) -> dict[str, rows]
    :type  parse_fn: Callable
    :param file_path: path to the XML file to be processed
    :type  file_path: Path
    :yield: table-tagged rows
    :rtype: Generator[DataItemWithMeta, Any]
    """
    logger.info("Reading from %s", str(file_path))
    n_entries = -1
    buffer = DictBuffer(max_items=settings.buffer_size)
    for n_entries, element in enumerate(stream_xml_file(file_path, xml_tag)):
        parsed_element = parse_fn(entry=element, file_path=file_path)
        yield from buffer.add_items(parsed_element)

        _log_interval_progress(n_entries, settings.log_interval)

    _log_interval_final(n_entries, file_path, settings.log_interval)
    yield from buffer.flush()


def process_xml_file_batches(
    settings: CtsSettings,
    xml_tag: str,
    parse_fn: Callable,
) -> Generator[DataItemWithMeta, Any]:
    """Generate a list of XML files to process using the NumericFileSequenceBatcher.

    :param settings: pipeline config with input_dir
    :type settings: CtsSettings
    :param xml_tag: XML element tag to stream
    :type xml_tag: str
    :param parse_fn: function for parsing the XML
    :type parse_fn: Callable
    """
    for files in get_file_batches(settings):
        for file_path in files:
            yield from process_xml_file(
                settings=settings,
                xml_tag=xml_tag,
                parse_fn=parse_fn,
                file_path=file_path,
            )


def process_xml_file_items(
    items: Iterable[FileItemDict],
    settings: CtsSettings,
    xml_tag: str,
    parse_fn: Callable,
) -> Generator[DataItemWithMeta, Any]:
    """Process each file item emitted by the dlt filesystem source.

    Each file is passed to the parser as `settings.input_dir` joined to the item's relative path,
    so recorded source file names do not depend on the absolute location of the input directory.

    :param items: file items to read
    :type  items: Iterable[FileItemDict]
    :param settings: pipeline config with input_dir, buffer_size and log_interval
    :type  settings: CtsSettings
    :param xml_tag: XML element tag to stream
    :type  xml_tag: str
    :param parse_fn: function for parsing the XML
    :type  parse_fn: Callable
    :yield: table-tagged rows
    :rtype: Generator[DataItemWithMeta, Any]
    """
    for file_item in items:
        yield from process_xml_file(
            settings=settings,
            xml_tag=xml_tag,
            parse_fn=parse_fn,
            file_path=settings.input_dir / file_item["relative_path"],
        )


def build_xml_file_resource(
    settings: CtsSettings,
    xml_tag: str,
    parse_fn: Callable,
    resource_name: str,
) -> DltResource:
    """Build a dlt resource that parses the XML files in `settings.input_dir`.

    Files are read using the dlt filesystem source.

    :param settings: pipeline config with input_dir, buffer_size and log_interval
    :type  settings: CtsSettings
    :param xml_tag: XML element tag to stream
    :type  xml_tag: str
    :param parse_fn: function for parsing the XML
    :type  parse_fn: Callable
    :param resource_name: name of the resulting dlt resource
    :type  resource_name: str
    :return: resource yielding table-tagged rows
    :rtype: DltResource
    """
    files = filesystem(bucket_url=str(settings.input_dir), file_glob=DEFAULT_XML_FILE_GLOB)
    resource = dlt.transformer(
        process_xml_file_items,
        data_from=files,
        name=resource_name,
        file_format="parquet",
        parallelized=True,
    )
    resource.bind(settings, xml_tag, parse_fn)
    return resource
