"""Shared machinery for the JSONL registry-ingest pipelines.

Each pipeline reads one JSONL input subdirectory per table, routes valid records
to a table named after the entity, and routes failures to ``<table>_rejected``
(or ``_invalid``) with provenance columns. This module holds the pieces all
three pipelines share; each pipeline module supplies only its validator and
registry loader.
"""

from collections.abc import Callable, Generator
from logging import Logger, getLogger
from typing import TYPE_CHECKING, Any, Final

import dlt
from dlt.common.pipeline import LoadInfo
from dlt.common.schema.typing import TAnySchemaColumns, TSchemaContractDict
from dlt.sources.filesystem import filesystem

from cdm_data_loaders.pipelines.core import run_pipeline
from cdm_data_loaders.readers.jsonlines import stream_jsonl_lines
from cdm_data_loaders.utils.buffer import ListBuffer

if TYPE_CHECKING:
    from dlt.extract.resource import DltResource

logger: Logger = getLogger(__name__)

SCHEMA_CONTRACT: Final[TSchemaContractDict] = {
    "tables": "evolve",
    "columns": "evolve",
    "data_type": "discard_row",
}

DEFAULT_REJECTED_SUFFIX: Final[str] = "_rejected"

REJECTED_COLUMNS: Final[tuple[str, ...]] = ("source_file", "line_no", "raw_record")

ParseErrorResolver = Callable[[dict[str, Any]], str | None]
ResourceFactory = Callable[[str, dict[str, Any], Any], "DltResource"]


def make_page_router(  # noqa: PLR0913
    table_name: str,
    buffer_size: int,
    resolve_error: ParseErrorResolver,
    *,
    rejected_suffix: str = DEFAULT_REJECTED_SUFFIX,
    include_record: bool = False,
    transform_valid: Callable[[Any], Any] | None = None,
) -> Callable[[list[dict[str, Any]]], Generator[Any, Any, Any]]:
    """Build a function that routes pages of line records to valid and rejected tables.

    For each item, `resolve_error` inspects the record and returns None for a
    valid record, or the error detail string otherwise. Parse errors bypass
    `resolve_error` entirely. Valid records yield `item["record"]` passed
    through `transform_valid` (identity when None) into `table_name`; failed
    records yield a provenance row into ``{table_name}{rejected_suffix}``. With
    `include_record`, the rejected row also carries the parsed `record` value.

    :param table_name: table name for valid records
    :type  table_name: str
    :param buffer_size: number of rows to buffer per table before yielding a page
    :type  buffer_size: int
    :param resolve_error: function returning None for valid records or error detail
    :type  resolve_error: ParseErrorResolver
    :param rejected_suffix: suffix for the rejected table name
    :type  rejected_suffix: str
    :param include_record: whether rejected rows carry the parsed record value
    :type  include_record: bool
    :param transform_valid: function applied to each valid record before yielding
    :type  transform_valid: Callable[[Any], Any] | None
    :return: function taking one page of line records and yielding routed pages
    :rtype : Callable[[list[dict[str, Any]]], Generator[Any, Any, Any]]
    """
    rejected_table = f"{table_name}{rejected_suffix}"

    def _route(page: list[dict[str, Any]]) -> Generator[Any, Any, Any]:
        valid_buffer = ListBuffer(table_name=table_name, max_items=buffer_size)
        rejected_buffer = ListBuffer(table_name=rejected_table, max_items=buffer_size)
        for item in page:
            error_detail = item["parse_error"]
            if error_detail is None:
                error_detail = resolve_error(item["record"])
            if error_detail is not None:
                rejected = {column: item[column] for column in REJECTED_COLUMNS}
                rejected["error_detail"] = error_detail
                if include_record:
                    rejected["record"] = item["record"]
                yield from rejected_buffer.add_item(rejected)
            else:
                record = item["record"]
                if transform_valid is not None:
                    record = transform_valid(record)
                yield from valid_buffer.add_item(record)
        yield from valid_buffer.flush()
        yield from rejected_buffer.flush()

    return _route


def build_entity_resource(  # noqa: PLR0913
    table_name: str,
    settings: Any,  # noqa: ANN401
    resolve_error: ParseErrorResolver,
    *,
    rejected_suffix: str = DEFAULT_REJECTED_SUFFIX,
    include_record: bool = False,
    transform_valid: Callable[[Any], Any] | None = None,
    table_columns: TAnySchemaColumns | None = None,
) -> "DltResource":
    """Read plain or gzip JSONL files for one entity and route validated pages.

    A missing entity directory produces no rows. Nested values are preserved
    rather than normalized into child tables.

    `table_columns` optionally sets the entity table's schema. It is applied as
    a dlt table variant, so the columns land on `table_name` only; the rejected
    table keeps dlt's column inference for its provenance rows. Accepted forms
    are a column dict, a list of column schemas, a pydantic model, or a
    pyarrow schema.

    :param table_name: table name; also the input_dir subdirectory name
    :type  table_name: str
    :param settings: pipeline configuration with input_dir, file_glob, buffer_size
    :type  settings: Any
    :param resolve_error: function returning None for valid records or error detail
    :type  resolve_error: ParseErrorResolver
    :param rejected_suffix: suffix for the rejected table name
    :type  rejected_suffix: str
    :param include_record: whether rejected rows carry the parsed record value
    :type  include_record: bool
    :param transform_valid: function applied to each valid record before yielding
    :type  transform_valid: Callable[[Any], Any] | None
    :param table_columns: optional column schema for the entity table
    :type  table_columns: TAnySchemaColumns | None
    :return: resource yielding validated and rejected records for this entity
    :rtype: DltResource
    """
    entity_dir = settings.input_dir / table_name
    if entity_dir.is_dir():
        files = filesystem(bucket_url=str(entity_dir), file_glob=settings.file_glob)
    else:
        logger.warning("Entity directory does not exist, skipping: %s", entity_dir)
        files = dlt.resource([], name=f"{table_name}_files")

    raw = dlt.transformer(
        stream_jsonl_lines,
        data_from=files,
        name=f"{table_name}_raw",
        max_table_nesting=0,
        parallelized=True,
    )
    raw.bind(settings.buffer_size)
    validated = dlt.transformer(
        make_page_router(
            table_name,
            settings.buffer_size,
            resolve_error,
            rejected_suffix=rejected_suffix,
            include_record=include_record,
            transform_valid=transform_valid,
        ),
        data_from=raw,
        name=f"{table_name}_validated",
        max_table_nesting=0,
        parallelized=True,
    )
    validated.apply_hints(schema_contract=TSchemaContractDict(**SCHEMA_CONTRACT))  # pyright: ignore[reportCallIssue]
    if table_columns is not None:
        # pipeline.run applies loader_file_format to the base hints, but a variant shadows them
        # for the entity table, so the variant must carry the file format itself
        validated.apply_hints(  # pyright: ignore[reportCallIssue]
            table_name=table_name,
            columns=table_columns,
            file_format=str(settings.loader_file_format),
            create_table_variant=True,
        )
    return validated


def run_registry_ingest(
    settings: Any,  # noqa: ANN401
    table_to_schema: dict[str, dict[str, Any]],
    build_resource: ResourceFactory,
    pipeline_name: str,
    pipeline_run_kwargs: dict[str, Any] | None = None,
) -> LoadInfo | None:
    """Load the requested tables using the shared pipeline runner.

    Defaults to every table in the registry when `settings.table_names` is None.
    Raises ValueError for unknown table names.

    :param settings: pipeline configuration
    :type  settings: Any
    :param table_to_schema: mapping of table name to its schema entry
    :type  table_to_schema: dict[str, dict[str, Any]]
    :param build_resource: function building one entity's resource
    :type  build_resource: ResourceFactory
    :param pipeline_name: registered pipeline name
    :type  pipeline_name: str
    :param pipeline_run_kwargs: extra kwargs for pipeline.run
    :type  pipeline_run_kwargs: dict[str, Any] | None
    :raises ValueError: if settings.table_names includes a name not in the registry
    :return: load information for the pipeline
    :rtype:  LoadInfo | None
    """
    table_names = settings.table_names or sorted(table_to_schema)
    unknown = sorted(set(table_names) - set(table_to_schema))
    if unknown:
        err_msg = f"Unknown table name(s) requested: {unknown}; available: {sorted(table_to_schema)}"
        raise ValueError(err_msg)
    resources = [build_resource(table_name, table_to_schema[table_name], settings) for table_name in table_names]
    return run_pipeline(
        settings=settings,
        resource=resources,
        pipeline_kwargs={
            "pipeline_name": pipeline_name,
            "dataset_name": settings.dataset_name or pipeline_name,
        },
        pipeline_run_kwargs=pipeline_run_kwargs or {},
    )
