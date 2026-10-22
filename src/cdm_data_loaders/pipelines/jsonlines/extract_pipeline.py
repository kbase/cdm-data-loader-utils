"""JSONL load pipeline for the KBase CTS.

Reads JSONL files. No validation occurs other than that records are valid JSON.

`table_name` names the destination table for all valid records; files
may live anywhere under `input_dir`. Records that fail JSON parsing go
to a `<table_name>_invalid` table. Each invalid record includes source
file, line number, raw text, and parse error detail.
"""

from collections.abc import Generator, Iterator
from typing import Any

import dlt
from dlt.common.pipeline import LoadInfo
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.common.typing import TDataItems

from cdm_data_loaders.pipelines.core import (
    PipelineCli,
    build_file_ingest_kwargs,
    build_filesystem_resource,
    run_pipeline,
)
from cdm_data_loaders.pipelines.jsonlines.settings import (
    EXTRACT_PIPELINE_NAME,
    JsonlExtractSettings,
)
from cdm_data_loaders.readers.jsonlines import route_jsonl_lines


@dlt.transformer(name="jsonl_reader", parallelized=True)
def jsonl_reader(items: Iterator[FileItemDict], settings: JsonlExtractSettings) -> Generator[TDataItems, Any, Any]:
    """Read each file in items and route each line by the configured table name.

    :param items: file items to read
    :type  items: Iterator[FileItemDict]
    :yield: one dict per line
    :rtype: Generator[TDataItems, Any, Any]
    """
    yield from route_jsonl_lines(
        items, buffer_size=settings.buffer_size, table_name_of=lambda _file_item: settings.table_name
    )


def run_jsonlines_pipeline(settings: JsonlExtractSettings) -> LoadInfo | None:
    """Run the JSONL pipeline.

    :param settings: pipeline configuration
    :type settings: JsonlExtractSettings
    :return: load information for the pipeline
    :rtype: LoadInfo | None
    """
    jsonl_resource = build_filesystem_resource(settings, jsonl_reader, settings, name="jsonl_reader")

    return run_pipeline(
        settings=settings,
        resource=jsonl_resource,
        **build_file_ingest_kwargs(settings, EXTRACT_PIPELINE_NAME, fallback_dataset_name=True),
    )


cli = PipelineCli(JsonlExtractSettings, run_jsonlines_pipeline)


if __name__ == "__main__":
    cli()
