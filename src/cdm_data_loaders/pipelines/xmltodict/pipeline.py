"""XML to dictionary pipeline for the KBase CTS."""

from collections.abc import Generator, Iterator
from pathlib import Path
from typing import Any

from dlt.common.pipeline import LoadInfo
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.common.typing import TDataItems

from cdm_data_loaders.pipelines.core import (
    PipelineCli,
    build_file_ingest_kwargs,
    build_filesystem_resource,
    run_pipeline,
)
from cdm_data_loaders.pipelines.xmltodict.settings import PIPELINE_NAME, XmlToDictSettings
from cdm_data_loaders.readers.xml import process_xml_file_to_dict


def _read_xmltodict_items(
    items: Iterator[FileItemDict], settings: XmlToDictSettings
) -> Generator[TDataItems, Any, Any]:
    """Read each file in items. Yields one dict per `xml_tag` element.

    :param items: file items to read
    :type  items: Iterator[FileItemDict]
    :param settings: pipeline config, including xml_tag
    :type  settings: XmlToDictSettings
    :yield: one dict per matched element
    :rtype: Generator[TDataItems, Any, Any]
    """
    for file_item in items:
        yield from process_xml_file_to_dict(settings, file_path=Path(file_item.local_file_path))


def run_xmltodict_pipeline(settings: XmlToDictSettings) -> LoadInfo | None:
    """Run the XML to dict pipeline on matching files in settings.input_dir.

    :param settings: pipeline configuration
    :type  settings: XmlToDictSettings
    """
    xmltodict_resource = build_filesystem_resource(settings, _read_xmltodict_items, settings, name="xmltodict_reader")

    return run_pipeline(
        settings=settings,
        resource=xmltodict_resource,
        **build_file_ingest_kwargs(settings, PIPELINE_NAME),
    )


cli = PipelineCli(XmlToDictSettings, run_xmltodict_pipeline)


if __name__ == "__main__":
    cli()
