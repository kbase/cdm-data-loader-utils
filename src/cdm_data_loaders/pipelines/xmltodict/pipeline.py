"""XML to dictionary pipeline for the KBase CTS."""

from collections.abc import Generator, Iterator
from logging import Logger, getLogger
from pathlib import Path
from typing import TYPE_CHECKING, Any

import dlt
from dlt.common.pipeline import LoadInfo
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.common.typing import TDataItems

from cdm_data_loaders.pipelines.core import PipelineCli, filesystem_resource, run_pipeline
from cdm_data_loaders.pipelines.xmltodict.settings import PIPELINE_NAME, XmlToDictSettings
from cdm_data_loaders.readers.xml import process_xml_file_to_dict

if TYPE_CHECKING:
    from dlt.extract import DltResource

logger: Logger = getLogger(__name__)


def _read_items(items: Iterator[FileItemDict], settings: XmlToDictSettings) -> Generator[TDataItems, Any, Any]:
    """Read each file in items. Yields one dict per `xml_tag` element.

    :param items: file items to read
    :type  items: Iterator[FileItemDict]
    :param settings: pipeline config, including xml_tag and xsd-derived paths
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
    reader = dlt.transformer(_read_items, name="xmltodict_reader", parallelized=True)

    files = filesystem_resource(bucket_url=settings.input_dir, file_glob=settings.file_glob)

    xmltodict_resource: DltResource = files | reader(settings)

    return run_pipeline(
        settings=settings,
        resource=xmltodict_resource,
        pipeline_kwargs={
            "pipeline_name": PIPELINE_NAME,
            "dataset_name": settings.dataset_name,
        },
        pipeline_run_kwargs={
            "loader_file_format": str(settings.loader_file_format),
        },
    )


cli = PipelineCli(XmlToDictSettings, run_xmltodict_pipeline)


if __name__ == "__main__":
    cli()
