"""XML to pyarrow pipeline for the KBase CTS."""

from collections.abc import Generator, Iterator
from pathlib import Path
from typing import TYPE_CHECKING, Any

import dlt
from dlt.common.pipeline import LoadInfo
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.common.typing import TDataItems
from dlt.sources.filesystem import filesystem

from cdm_data_loaders.pipelines.core import PipelineCli, run_pipeline
from cdm_data_loaders.pipelines.xmltopyarrow.settings import PIPELINE_NAME, XmlToPyArrowSettings
from cdm_data_loaders.readers.xml import process_xml_file_to_pyarrow

if TYPE_CHECKING:
    from dlt.extract import DltResource


def _read_items(items: Iterator[FileItemDict], settings: XmlToPyArrowSettings) -> Generator[TDataItems, Any, Any]:
    """Read each file in items. Yields PyArrow tables for matched `xml_tag` elements.

    :param items: file items to read
    :type  items: Iterator[FileItemDict]
    :param settings: pipeline config, including xml_tag and XSD schema
    :type  settings: XmlToPyArrowSettings
    :yield: PyArrow tables of matched elements
    :rtype: Generator[TDataItems, Any, Any]
    """
    for file_item in items:
        yield from process_xml_file_to_pyarrow(settings, file_path=Path(file_item.local_file_path))


def run_xmltopyarrow_pipeline(settings: XmlToPyArrowSettings) -> LoadInfo | None:
    """Run the XML to PyArrow pipeline on matching files in settings.input_dir.

    :param settings: pipeline configuration
    :type  settings: XmlToPyArrowSettings
    """
    reader = dlt.transformer(_read_items, name="xmltopyarrow_reader", parallelized=True)
    files = filesystem(bucket_url=str(settings.input_dir), file_glob=settings.file_glob)

    xmltopyarrow_resource: DltResource = files | reader(settings)
    pipeline_name = str(settings.model_config.get("cli_prog_name") or PIPELINE_NAME)

    return run_pipeline(
        settings=settings,
        resource=xmltopyarrow_resource,
        pipeline_kwargs={
            "pipeline_name": pipeline_name,
            "dataset_name": settings.dataset_name,
        },
        pipeline_run_kwargs={
            "loader_file_format": str(settings.loader_file_format),
        },
    )


cli = PipelineCli(XmlToPyArrowSettings, run_xmltopyarrow_pipeline)


if __name__ == "__main__":
    cli()
