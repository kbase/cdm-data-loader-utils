"""Post-load compaction for xml2db "reused" tables written to a filesystem destination.

`cdm_data_loaders.readers.xml2db_doc.process_xml_file_with_xml2db` parses each source file --
and, when chunked parsing is enabled (`Xml2DbSettings.xml_tag`), each chunk within a
file -- into its own independent xml2db `Document`. xml2db's content-hash deduplication of
"reused" tables only ever happens *within* a single `Document`, so identical content parsed by
two different `Document` instances -- two separate source files, two chunks of one file, or any
combination -- still produces two physically separate rows once loaded.

`cdm_data_loaders.readers.xml2db_doc.flatten_xml2db_document` already gives those rows the *same*
content-addressed primary key (see its `_content_addressed_key` helper), so any literal
duplicates are byte-for-byte identical, and every foreign key referencing them is already
correct -- nothing needs to be rewritten. This module simply collapses those literal duplicates
after a pipeline run, for the (jsonl or parquet) filesystem destination, by keeping one arbitrary
row per distinct primary key value.
"""

import gzip
import json
import logging
from pathlib import Path
from typing import Final

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq
from xml2db import DataModel

from cdm_data_loaders.core.fields import PARQUET
from cdm_data_loaders.readers.xml2db_doc import is_xml2db_content_addressed_table

logger: logging.Logger = logging.getLogger(__name__)

_PARQUET_GLOB: Final[str] = "*.parquet"
_JSON_GLOB: Final[str] = "*.json*"
_COMPACTED_FILE_STEM: Final[str] = "compacted"


def reused_table_names(model: DataModel) -> list[str]:
    """Return the destination table names of every table whose rows use a content-addressed key.

    This intentionally excludes the schema's root table even though it is itself "reused": see
    `cdm_data_loaders.readers.xml2db_doc.is_xml2db_content_addressed_table` for why the root's
    key isn't content-addressed, and therefore never produces the kind of duplicate this module
    cleans up.

    :param model: an xml2db `DataModel`.
    :type model: DataModel
    :return: destination table names (`DataModelTable.name`) of every content-addressed table.
    :rtype: list[str]
    """
    return [table.name for table in model.tables.values() if is_xml2db_content_addressed_table(model, table)]


def compact_reused_tables(dataset_dir: Path, model: DataModel, loader_file_format: str) -> dict[str, int]:
    """Deduplicate literal duplicate rows out of every "reused" table's output files, in place.

    Reads and dedups ALL of a table's output fragment files together, in one pass -- across
    every source file and chunk the pipeline run produced, not just within one of them. Safe to
    call on a dataset that has no duplicates at all: tables with nothing to remove are left
    completely untouched.

    :param dataset_dir: the pipeline's output dataset directory (`output_dir / dataset_name`),
        containing one subdirectory per table, as produced by dlt's `filesystem` destination.
    :type dataset_dir: Path
    :param model: the xml2db `DataModel` used to produce the dataset; used to find every
        "reused" table's destination name and primary key column (`pk_<table name>`).
    :type model: DataModel
    :param loader_file_format: `"parquet"` or `"jsonl"`, matching the pipeline's
        `loader_file_format` setting.
    :type loader_file_format: str
    :return: mapping of table name to number of literal duplicate rows removed. Tables that were
        never written, or that had no duplicates, are omitted.
    :rtype: dict[str, int]
    """
    connection = duckdb.connect()
    removed: dict[str, int] = {}
    for table_name in reused_table_names(model):
        table_dir = dataset_dir / table_name
        if not table_dir.is_dir():
            continue
        n_removed = _compact_table(connection, table_dir, table_name, loader_file_format)
        if n_removed:
            removed[table_name] = n_removed
    return removed


def _source_files(table_dir: Path, loader_file_format: str) -> list[Path]:
    """List a table's current output files, in a stable order."""
    glob_pattern = _PARQUET_GLOB if loader_file_format == PARQUET else _JSON_GLOB
    return sorted(f for f in table_dir.glob(glob_pattern) if f.is_file())


def _read_expr(source_files: list[Path], loader_file_format: str) -> str:
    """Build a DuckDB table-function expression that reads every one of a table's source files.

    Table/column identifiers used elsewhere in this module's SQL come from the xml2db model
    (never from XML document content), so building these queries by string formatting carries no
    injection risk; only the file paths interpolated here originate from the filesystem, hence
    the defensive quote-escaping.
    """
    read_fn = "read_parquet" if loader_file_format == PARQUET else "read_json_auto"
    quoted_paths = ", ".join("'" + f.as_posix().replace("'", "''") + "'" for f in source_files)
    return f"{read_fn}([{quoted_paths}])"


def _compact_table(
    connection: duckdb.DuckDBPyConnection,
    table_dir: Path,
    table_name: str,
    loader_file_format: str,
) -> int:
    """Deduplicate one table's output files in place; returns the number of rows removed."""
    source_files = _source_files(table_dir, loader_file_format)
    if not source_files:
        return 0

    read_expr = _read_expr(source_files, loader_file_format)
    pk_column = f"pk_{table_name}"

    total_rows: int = connection.execute(f"SELECT count(*) FROM {read_expr}").fetchone()[0]  # noqa: S608
    deduped: pa.Table = connection.execute(
        f"SELECT * FROM {read_expr} QUALIFY row_number() OVER (PARTITION BY {pk_column}) = 1"  # noqa: S608
    ).to_arrow_table()

    n_removed = total_rows - deduped.num_rows
    if n_removed == 0:
        return 0

    for source_file in source_files:
        source_file.unlink()
    _write_compacted(deduped, table_dir, loader_file_format)

    logger.info("Compacted %d duplicate row(s) out of table %s", n_removed, table_name)
    return n_removed


def _write_compacted(deduped: pa.Table, table_dir: Path, loader_file_format: str) -> None:
    """Write a deduplicated pyarrow table as a table's sole remaining output file."""
    if loader_file_format == PARQUET:
        pq.write_table(deduped, table_dir / f"{_COMPACTED_FILE_STEM}.parquet")
        return

    with gzip.open(table_dir / f"{_COMPACTED_FILE_STEM}.jsonl.gz", "wt", encoding="utf-8") as fh:
        for row in deduped.to_pylist():
            fh.write(json.dumps(row, default=str))
            fh.write("\n")
