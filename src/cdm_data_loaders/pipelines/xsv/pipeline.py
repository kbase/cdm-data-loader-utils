"""XSV (CSV/TSV) ingestion pipeline for the KBase CTS.

Cleans and validates each XSV file against `settings.schema_file` using qsv (see
`readers.jsonschema_xsv.xsv_validator`), then loads each valid row, cast to the types declared
in the schema, into `table_name`. Files that fail cleaning/validation are routed to
`<table_name>_rejected`, one row per recorded error; qsv's own per-row validation diagnostics are
left in its scratch output directory and are not currently surfaced as rows.
"""

from collections.abc import Generator, Iterator
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Final

import dlt
from dlt.common.pipeline import LoadInfo
from dlt.common.schema.typing import TSchemaContractDict
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.common.typing import TDataItems

from cdm_data_loaders.pipelines.core import PipelineCli, filesystem_resource, run_pipeline
from cdm_data_loaders.pipelines.xsv.settings import PIPELINE_NAME, XsvIngestSettings
from cdm_data_loaders.readers.jsonschema_xsv.xsv_reader import read_validated_rows
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.helpers import CleanerValidatorArgs, ErrorRecord
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.qsv import clean_validate_file, qsv_check
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import generate_header
from cdm_data_loaders.utils.buffer import ListBuffer

SCHEMA_CONTRACT: Final[TSchemaContractDict] = {
    "tables": "evolve",
    "columns": "evolve",
    "data_type": "discard_row",
}
REJECTED_SUFFIX: Final[str] = "_rejected"


@dataclass
class XsvWorkPaths:
    """Scratch resources shared across all files processed in one pipeline run.

    :param qsv_cmd: path to the qsv binary
    :type qsv_cmd: str
    :param header_file_path: header file used when an input file is missing a header row
    :type header_file_path: Path
    :param tmp_dir: scratch directory for intermediate cleaning/validation files
    :type tmp_dir: Path
    :param qsv_output_dir: directory qsv writes failure diagnostics to
    :type qsv_output_dir: Path
    :param validated_dir: directory the final cleaned, validated file is written to
    :type validated_dir: Path
    """

    qsv_cmd: str
    header_file_path: Path
    tmp_dir: Path
    qsv_output_dir: Path
    validated_dir: Path


def _process_xsv_file(
    file_item: FileItemDict,
    settings: XsvIngestSettings,
    work_paths: XsvWorkPaths,
) -> Generator[TDataItems, Any, Any]:
    """Clean, validate, and load one XSV file. Yield rows for table_name and its rejected table.

    :param file_item: raw XSV file as a FileItemDict
    :type  file_item: FileItemDict
    :param settings: pipeline configuration
    :type  settings: XsvIngestSettings
    :param work_paths: scratch resources shared across the pipeline run
    :type  work_paths: XsvWorkPaths
    :yield: rows routed to table_name and/or ``{table_name}_rejected``
    :rtype: Generator[TDataItems, Any, Any]
    """
    file_path = Path(settings.input_dir) / file_item["relative_path"]
    args = CleanerValidatorArgs(
        qsv_cmd=work_paths.qsv_cmd,
        xsv_file_path=file_path,
        header_file_path=work_paths.header_file_path,
        first_pass_schema=settings.first_pass_schema,
        post_norm_schema=settings.validated_schema,
        tmp_dir_path=work_paths.tmp_dir,
        qsv_output_dir_path=work_paths.qsv_output_dir,
        validated_file_dir_path=work_paths.validated_dir,
        **settings.parsing_config,
    )

    validated_file_name = clean_validate_file(args)

    valid_buffer = ListBuffer(table_name=settings.table_name, max_items=settings.buffer_size)
    rejected_buffer = ListBuffer(table_name=f"{settings.table_name}{REJECTED_SUFFIX}", max_items=settings.buffer_size)

    if validated_file_name:
        properties = settings.validated_schema.jsonschema.get("properties", {})
        rows = read_validated_rows(
            work_paths.validated_dir / validated_file_name,
            properties,
            delimiter=args.delimiter,
            quote=args.quote,
            escape=args.escape,
        )
        for row in rows:
            yield from valid_buffer.add_item(row)
    elif not args.errors:
        args.errors.append(
            ErrorRecord(
                file=args.file_name,
                message="File failed cleaning/validation with no rows surviving; see qsv output directory for details.",
            )
        )

    for error in args.errors:
        yield from rejected_buffer.add_item(error.model_dump())

    yield from valid_buffer.flush()
    yield from rejected_buffer.flush()


def _read_items(
    items: Iterator[FileItemDict],
    settings: XsvIngestSettings,
    work_paths: XsvWorkPaths,
) -> Generator[TDataItems, Any, Any]:
    """Clean, validate, and load each file in file_paths.

    :param items: file items to read
    :type  items: Iterator[FileItemDict]
    :param settings: pipeline configuration
    :type  settings: XsvIngestSettings
    :param work_paths: scratch resources shared across the pipeline run
    :type  work_paths: XsvWorkPaths
    :yield: rows routed to table_name and/or ``{table_name}_rejected``
    :rtype: Generator[TDataItems, Any, Any]
    """
    for file_item in items:
        yield from _process_xsv_file(file_item, settings, work_paths)


def run_xsv_pipeline(settings: XsvIngestSettings) -> LoadInfo | None:
    """Run the XSV ingestion pipeline on matching files in settings.input_dir.

    :param settings: pipeline configuration
    :type  settings: XsvIngestSettings
    :return: load information for the pipeline
    :rtype: LoadInfo | None
    """
    qsv_cmd = qsv_check()

    xsv_reader = dlt.transformer(_read_items, name="xsv_reader", parallelized=True)

    files = filesystem_resource(bucket_url=settings.input_dir, file_glob=settings.file_glob)

    with (
        TemporaryDirectory(prefix="xsv_tmp_") as tmp_dir,
        TemporaryDirectory(prefix="xsv_qsv_output_") as qsv_output_dir,
        TemporaryDirectory(prefix="xsv_validated_") as validated_dir,
    ):
        work_paths = XsvWorkPaths(
            qsv_cmd=qsv_cmd,
            header_file_path=generate_header(settings.validated_schema, Path(tmp_dir)),
            tmp_dir=Path(tmp_dir),
            qsv_output_dir=Path(qsv_output_dir),
            validated_dir=Path(validated_dir),
        )

        xsv_resource = files | xsv_reader(settings, work_paths)
        xsv_resource.apply_hints(schema_contract=TSchemaContractDict(**SCHEMA_CONTRACT))  # pyright: ignore[reportCallIssue]

        return run_pipeline(
            settings=settings,
            resource=xsv_resource,
            pipeline_kwargs={
                "pipeline_name": PIPELINE_NAME,
                "dataset_name": settings.dataset_name,
            },
            pipeline_run_kwargs={
                "loader_file_format": str(settings.loader_file_format),
            },
        )


cli = PipelineCli(XsvIngestSettings, run_xsv_pipeline)


if __name__ == "__main__":
    cli()
