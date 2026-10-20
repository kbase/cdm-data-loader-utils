"""Misc helpers for integration tests."""

import gzip
import json
from collections import defaultdict
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from datetime import datetime
from logging import Logger, getLogger
from pathlib import Path
from typing import Any, Final

import dlt
import duckdb
import pendulum
import pyarrow.parquet as pq
from dlt.common.pipeline import LoadInfo
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.destinations import filesystem

from cdm_data_loaders.core.fields import PARQUET, LoaderFileFormatEnum
from cdm_data_loaders.core.settings import CtsSettings
from cdm_data_loaders.pipelines.core import LOAD_INFO_TABLE_NAME
from tests.conftest import LOCAL_FS

logger: Logger = getLogger(__name__)


WORKER_CONFIGS = ((1, 1, 1), (3, 2, 2), (5, 3, 4))
BUFFER_SIZES = (1, 5, 100, 500)
LOADER_FILE_FORMATS = list(LoaderFileFormatEnum)

DLT_METADATA_PREFIX: Final[str] = "_dlt"

XMLNS_COLUMNS: Final[frozenset[str]] = frozenset({"entry__axmlns", "entry__axmlns_xsi"})

DLT_TO_XML_NAMES: Final[dict[str, str]] = {
    "aid": "@id",
    "aupdated": "@updated",
    "atype": "@type",
    "avalue": "@value",
    "alength": "@length",
    "achecksum": "@checksum",
    "_text": "#text",
    "representative_member": "representativeMember",
    "db_reference": "dbReference",
}


def is_metadata_table(table_name: str) -> bool:
    """Return whether a table contains dlt pipeline metadata rather than loaded records."""
    return table_name.startswith("_dlt") or table_name == LOAD_INFO_TABLE_NAME


def read_jsonl_file(file_path: Path) -> list[dict[str, Any]]:
    """Read a plain or gzipped JSONLines file, detecting gzip by magic bytes."""
    with file_path.open("rb") as f:
        is_gzipped = f.read(2) == b"\x1f\x8b"

    opener = gzip.open if is_gzipped else open
    rows: list[dict[str, Any]] = []
    with opener(file_path, "rt", encoding="utf-8") as f:
        for raw_line in f:
            stripped = raw_line.strip()
            if not stripped:
                continue
            try:
                rows.append(json.loads(stripped))
            except json.JSONDecodeError:
                logger.warning("Skipping unparseable line in %s: %r", file_path, raw_line[:200])
    return rows


def read_parquet_tables(output_dir: Path) -> dict[str, list[dict[str, Any]]]:
    """Read all parquet files under a pipeline output dataset dir, grouped by table name."""
    tables: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for parquet_file in sorted(output_dir.rglob("*.parquet")):
        tables[parquet_file.parent.name].extend(pq.read_table(parquet_file).to_pylist())
    return dict(tables)


def read_pipeline_tables(output_dir: Path, output_format: str | None = None) -> dict[str, list[dict[str, Any]]]:
    """Read all json files under a pipeline output dataset dir, grouped by table name.

    Defaults to reading jsonlines files; supply "parquet" as file format to read parquet files.
    """
    if output_format and output_format == PARQUET:
        return read_parquet_tables(output_dir)

    tables: dict[str, list[dict[str, Any]]] = defaultdict(list)

    for json_file in sorted(output_dir.rglob("*.json*")):
        if not json_file.is_file():
            continue
        # ignore metadata dirs
        if is_metadata_table(json_file.parent.name):
            continue
        tables[json_file.parent.name].extend(read_jsonl_file(json_file))

    return dict(tables)


def read_destination_rows(base_dir: Path, table_name: str) -> list[dict[str, Any]]:
    """Read the rows the filesystem destination wrote for table_name under base_dir, without dlt columns.

    Handles parquet, jsonl and gzipped jsonl files. Returns an empty list if the table has no files.
    """
    rows: list[dict[str, Any]] = []
    for data_file in sorted(base_dir.glob(f"**/{table_name}/*")):
        if data_file.name.endswith(".parquet"):
            rows.extend(pq.read_table(data_file).to_pylist())
        elif ".jsonl" in data_file.name:
            rows.extend(read_jsonl_file(data_file))
    return [strip_metadata_columns(row) for row in rows]


def make_cts_settings[S: CtsSettings](
    settings_cls: type[S],
    base_dir: Path,
    **overrides: Any,  # noqa: ANN401
) -> S:
    """Create standard input/output directories and log config under base_dir and instantiate settings_cls."""
    input_dir = base_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir = base_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    log_config_file = base_dir / "logging.json"
    if not log_config_file.exists():
        log_config_file.write_text('{"version": 1}', encoding="utf-8")

    kwargs: dict[str, Any] = {
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "log_config_file": str(log_config_file),
        "use_destination": LOCAL_FS,
        "dlt_dev_mode": False,
        "use_output_dir_for_pipeline_metadata": False,
    }
    kwargs.update(overrides)
    return settings_cls(**kwargs)


def make_file_item(file_path: Path, mime_type: str = "application/octet-stream") -> FileItemDict:
    """Build a dlt FileItemDict for a local file."""
    return FileItemDict(
        {
            "file_url": file_path.as_uri(),
            "file_name": file_path.name,
            "relative_path": file_path.name,
            "mime_type": mime_type,
            "modification_date": pendulum.now(),
            "size_in_bytes": file_path.stat().st_size,
        }
    )


def fetch_all_duckdb(pipeline: dlt.Pipeline, query: str) -> list[tuple[Any, ...]]:
    """Run a query against a pipeline's DuckDB destination and return all rows."""
    with (
        pipeline.sql_client() as client,
        client.execute_query(query) as cur,
    ):
        return list(cur.fetchall())


@contextmanager
def duckdb_reference_from_jsonl(jsonl_path: Path) -> Iterator[duckdb.DuckDBPyConnection]:
    """Load a reference JSONL file into an in-memory DuckDB connection with a `reference` table."""
    connection = duckdb.connect()
    connection.execute(
        "CREATE TABLE reference AS SELECT * FROM read_json_auto(?, sample_size=200)",
        [str(jsonl_path)],
    )
    try:
        yield connection
    finally:
        connection.close()


def sorted_json(entries: list[dict[str, Any]]) -> list[str]:
    """Entries as sorted JSON strings for order-insensitive comparison."""
    return sorted(json.dumps(entry, sort_keys=True, default=str) for entry in entries)


def settings_to_argv(settings: Any) -> list[str]:
    """Derive CLI arguments from a settings instance.

    Every serialized field becomes a ``--kebab-case`` flag, with bools rendered
    as their ``true``/``false`` CLI spelling, sequence fields expanded after the
    flag, and None-valued or computed fields skipped.

    :param settings: a CtsSettings (or derivative) instance
    :type settings: Any
    :return: command-line argument list reproducing the settings
    :rtype: list[str]
    """
    argv: list[str] = []
    computed_fields: set[str] = set(getattr(type(settings), "model_computed_fields", {}).keys())
    for name, value in settings.model_dump(exclude=computed_fields).items():
        flag = f"--{name.replace('_', '-')}"
        if isinstance(value, bool):
            argv.extend([flag, "true" if value else "false"])
        elif value is None:
            continue
        elif isinstance(value, list | tuple):
            if value:
                argv.extend([flag, *(str(item) for item in value)])
        else:
            argv.extend([flag, str(value)])
    return argv


def assert_run_ok(load_info: LoadInfo | None) -> LoadInfo:
    """Assert a pipeline run succeeded, and return the load info for chained use."""
    assert load_info is not None
    assert not load_info.has_failed_jobs
    return load_info


def read_data_rows(dataset_dir: Path, file_format: str | None = None) -> dict[str, list[dict[str, Any]]]:
    """Read a pipeline output dataset, dropping dlt metadata columns and metadata tables."""
    tables = read_pipeline_tables(dataset_dir, file_format)
    return {
        name: [strip_metadata_columns(row) for row in rows]
        for name, rows in tables.items()
        if not is_metadata_table(name)
    }


def assert_pipeline_report_written(output_dir: Path, load_info: LoadInfo, *, metadata_moved: bool = False) -> None:
    """Assert that save_pipeline_metadata wrote the expected report file (and moved metadata if configured)."""
    report_files = list(output_dir.glob(f"report_{load_info.loads_ids[0]}__*.txt"))
    assert len(report_files) == 1
    pipeline_e2e_id = report_files[0].stem.removeprefix("report_")
    if metadata_moved:
        assert (output_dir / pipeline_e2e_id).is_dir()
    report_text = report_files[0].read_text(encoding="utf-8")
    report_lines = report_text.splitlines()
    started_line = report_lines[2]
    finished_line = report_lines[3]
    assert started_line.startswith("Started: ")
    assert finished_line.startswith("Finished: ")
    started_at = datetime.fromisoformat(started_line.removeprefix("Started: "))
    finished_at = datetime.fromisoformat(finished_line.removeprefix("Finished: "))
    assert pipeline_e2e_id == f"{load_info.loads_ids[0]}__{finished_at}"
    assert started_at <= load_info.started_at <= finished_at
    if load_info.finished_at is not None:
        assert load_info.started_at <= load_info.finished_at <= finished_at
    assert report_text == (
        f"Pipeline: {load_info.pipeline.pipeline_name}\n"
        f"Dataset: {load_info.dataset_name}\n"
        f"Started: {started_at}\n"
        f"Finished: {finished_at}\n"
        "Status: completed\n"
        f"Completed loads: {', '.join(load_info.loads_ids)}\n"
        "Pipeline load info saved\n"
    )


def run_and_read(
    run_pipeline: Any,
    chunk_dir: str,
    **overrides: Any,
) -> tuple[LoadInfo, Path]:
    """Run the pipeline on a chunk dir and return the run result with its tables."""
    (load_info, output_dir) = run_pipeline(chunk_dir, **overrides)
    assert load_info is not None
    assert not load_info.has_failed_jobs
    assert_pipeline_report_written(output_dir, load_info)
    return (load_info, output_dir)


def extract_table_data(dataset: dlt.Dataset, table_name: str) -> list[str]:
    """Convert table data into a JSON string, removing the _dlt fields first."""
    return sorted_json(
        [
            {k: v for k, v in datum.items() if not is_metadata_table(k)}
            for datum in dataset.table(table_name).df().to_dict("records")
        ]  # pyright: ignore[reportArgumentType]
    )


def get_pipeline(load_info: LoadInfo, output_dir: Path) -> dlt.Pipeline:
    """Return a Pipeline object for the given dataset and destination."""
    return dlt.pipeline(destination=filesystem(bucket_url=str(output_dir)), dataset_name=load_info.dataset_name)


def strip_metadata_columns(row: dict[str, Any]) -> dict[str, Any]:
    """Drop dlt metadata columns from a row dict."""
    return {key: value for key, value in row.items() if not key.startswith(DLT_METADATA_PREFIX)}


def denormalize_identifier(name: str) -> str:
    """Restore an xmltodict key from a dlt-normalized fragment."""
    return DLT_TO_XML_NAMES.get(name, name)


def denormalize_path(path: str) -> list[str]:
    """Split a dlt name into xml key fragments, unmangling each one."""
    return [denormalize_identifier(fragment) for fragment in path.split("__")]


def set_nested(target: dict[str, Any], path: list[str], value: Any) -> None:
    """Set a value at a nested path, creating intermediate dicts as needed."""
    for fragment in path[:-1]:
        target = target.setdefault(fragment, {})
    target[path[-1]] = value


def rows_as_multiset(tables: dict[str, list[dict[str, Any]]], table_name: str) -> list[str]:
    """Rows of a table as sorted JSON strings for order-insensitive comparison."""
    return sorted(
        json.dumps(strip_metadata_columns(row), sort_keys=True, default=str) for row in tables.get(table_name, [])
    )


def parse_json_string(value: Any) -> Any:
    """Parse a string holding a JSON object or array; pass any other value through.

    Parquet output with max_table_nesting=0 stores nested content as JSON strings,
    while jsonl output keeps it as parsed structures.
    """
    if isinstance(value, str) and value[:1] in ("{", "["):
        try:
            return json.loads(value)
        except json.JSONDecodeError:
            return value
    return value


def reconstruct_entries(tables: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    """Rebuild nested dicts from the pipeline's flattened output tables.

    The main table holds scalar columns; child tables hold list contents,
    linked by _dlt_parent_id and ordered by _dlt_list_idx. A child table's parent
    may be a main entry row or another nested row, so the walk follows parentage
    rather than table-name prefixes.
    """
    top_tables = [t for t in tables if "__" not in t]
    if not top_tables or len(top_tables) > 1:
        err_msg = "Could not find an appropriate top-level table for dataset"
        raise ValueError(err_msg)
    top_table = top_tables[0]

    main_rows = tables.get(top_table, [])
    if not main_rows:
        return []

    parent_index: dict[str, list[tuple[str, dict[str, Any]]]] = defaultdict(list)
    for table_name, rows in tables.items():
        if table_name == top_table or is_metadata_table(table_name):
            continue
        for row in rows:
            parent_index[row["_dlt_parent_id"]].append((table_name, row))

    entries: list[dict[str, Any]] = []
    for row in main_rows:
        entry: dict[str, Any] = {}
        for column, value in strip_metadata_columns(row).items():
            if column in XMLNS_COLUMNS or value is None:
                continue
            set_nested(entry, denormalize_path(column.removeprefix(f"{top_table}__")), parse_json_string(value))
        _attach_nested_rows(entry, row["_dlt_id"], parent_index, [], depth_limit=100, top_table=top_table)
        entries.append(entry)
    return entries


def _attach_nested_rows(  # noqa: PLR0917
    target: dict[str, Any],
    parent_id: str,
    parent_index: dict[str, list[tuple[str, dict[str, Any]]]],
    parent_root_path: list[str],
    depth_limit: int,
    top_table: str,
) -> None:
    """Attach child-table rows (and their own children, recursively) to a target dict.

    parent_root_path is the target's own path from the root entry; child paths are
    made relative to it, since a child table name encodes the path from the root.
    """
    if depth_limit <= 0:
        return
    grouped: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for table_name, row in parent_index.get(parent_id, []):
        grouped[table_name].append(row)

    for table_name, rows in grouped.items():
        rows.sort(key=lambda r: r["_dlt_list_idx"])
        full_path = _child_table_path(table_name, top_table)
        relative_path = full_path[len(parent_root_path) :] if parent_root_path else full_path
        children = []
        for row in rows:
            child: dict[str, Any] = {}
            for column, raw in strip_metadata_columns(row).items():
                if raw is None:
                    continue
                set_nested(child, denormalize_path(column), parse_json_string(raw))
            _attach_nested_rows(child, row["_dlt_id"], parent_index, full_path, depth_limit - 1, top_table)
            children.append(child)
        set_nested(target, relative_path, children)


def _child_table_path(table_name: str, top_table: str) -> list[str]:
    """Derive the xml key path for a child table from its dlt name.

    Table names are prefixed with the parent chain (entry__entry__property); the
    leading fragments repeat the root table name and are dropped.
    """
    fragments = denormalize_path(table_name)
    while len(fragments) > 1 and fragments[0] == top_table:
        fragments = fragments[1:]
    return fragments


def reference_entries_from_duckdb(connection: duckdb.DuckDBPyConnection) -> list[dict[str, Any]]:
    """Extract reference entries from DuckDB as canonical nested Python dicts."""
    rows = connection.execute("SELECT entry FROM reference ORDER BY entry['_id']").fetchall()
    return [_normalize_reference_entry(entry) for (entry,) in rows]


def _normalize_reference_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Normalize a DuckDB STRUCT entry dict to the shape of reconstructed pipeline entries."""
    normalized: dict[str, Any] = {}
    for key, raw_value in entry.items():
        if raw_value is None:
            continue
        value = raw_value
        if isinstance(value, str) and value.startswith(("{", "[")):
            value = json.loads(value)
        normalized[key] = _normalize_reference_value(value)
    return normalized


def _normalize_reference_value(value: Any) -> Any:
    """Recursively normalize DuckDB STRUCT/LIST/JSON values to plain Python structures."""
    if isinstance(value, list):
        return [_normalize_reference_value(item) for item in value]
    if isinstance(value, dict):
        return {key: _normalize_reference_value(item) for key, item in value.items() if item is not None}
    return value


def assert_datasets_equal(
    load_info_and_dir: dict[str, tuple[LoadInfo, Path]],
) -> None:
    """Compare datasets to check whether they are identical."""
    if len(load_info_and_dir) < 2:
        err_msg = "assert_datasets_equal requires at least two datasets to compare"
        raise ValueError(err_msg)

    datasets = {key: get_pipeline(*v) for key, v in load_info_and_dir.items()}

    all_keys = list(load_info_and_dir.keys())

    sorted_row_data = {k: {} for k in all_keys}
    dataset_zero = all_keys[0]
    table_names = {
        table_name for table_name in datasets[dataset_zero].dataset().tables if not is_metadata_table(table_name)
    }
    rows_per_table = {
        table_name: datasets[dataset_zero].dataset().table(table_name).df().shape[0] for table_name in table_names
    }
    col_names = {t: set(datasets[dataset_zero].dataset().table(t).df().columns) for t in table_names}
    sorted_row_data[dataset_zero] = {t: extract_table_data(datasets[dataset_zero].dataset(), t) for t in table_names}

    for key in all_keys[1:]:
        # check the table names are identical
        assert {
            table_name for table_name in datasets[key].dataset().tables if not is_metadata_table(table_name)
        } == table_names
        # check row counts
        assert {
            table_name: datasets[key].dataset().table(table_name).df().shape[0] for table_name in table_names
        } == rows_per_table
        # column names
        for t in table_names:
            assert set(datasets[key].dataset().table(t).df().columns) == col_names[t]
            sorted_row_data[key][t] = extract_table_data(datasets[key].dataset(), t)
            assert sorted_row_data[key][t] == sorted_row_data[dataset_zero][t]


def assert_dataset_matches_reference(
    load_info_and_dir: dict[str, tuple[LoadInfo, Path]],
    sorted_reference_json: list[str],
    output_format: str | None = None,
    reconstruct: Callable[[dict[str, list[dict[str, Any]]]], list[dict[str, Any]]] = reconstruct_entries,
    canonicalize: Callable[[dict[str, Any]], dict[str, Any]] | None = None,
) -> None:
    """Ensure that a dataset matches the reference dataset."""
    # parsed JSON output files for each table, indexed by key
    tables: dict[str, dict[str, Any]] = {}
    # original data structure, reconstructed from the tables, indexed by key
    reconstructed: dict[str, list[Any]] = {}

    for key, (load_info, output_dir) in load_info_and_dir.items():
        tables[key] = read_pipeline_tables(output_dir / load_info.dataset_name, output_format)
        entries = reconstruct(tables[key])
        if canonicalize:
            entries = [canonicalize(entry) for entry in entries]
        reconstructed[key] = sorted_json(entries)

    for value in reconstructed.values():
        assert value == sorted_reference_json
