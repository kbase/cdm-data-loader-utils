"""JSONLines file reading utilities."""

import json
from collections.abc import Callable, Generator, Iterator
from typing import Any, Literal

from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.common.typing import TDataItems

from cdm_data_loaders.core.fields import GZIP_SUFFIX
from cdm_data_loaders.utils.buffer import DictBuffer

LineRecord = dict[str, Any]


def compression_mode(file_name: str) -> Literal["enable", "disable"]:
    """Return the dlt compression mode for a file name.

    Returns "enable" for names ending in .gz and "disable" otherwise.
    Auto-detection is not used, since it depends on remote file metadata
    that local files do not carry.
    """
    return "enable" if file_name.endswith(GZIP_SUFFIX) else "disable"


def parse_jsonl_line(stripped: str) -> LineRecord:
    """Parse one stripped JSONL line into the shared line-record shape.

    The returned dict has: record (parsed JSON, or None on parse
    failure), raw_record, parse_error (None on success).
    """
    parse_error: str | None = None
    record: dict[str, Any] | None = None
    try:
        record = json.loads(stripped)
    except json.JSONDecodeError as e:
        parse_error = str(e)
    return {"record": record, "raw_record": stripped, "parse_error": parse_error}


def stream_jsonl_lines(items: Iterator[FileItemDict], buffer_size: int) -> Generator[TDataItems, Any, Any]:
    """Read each file in items. Yield pages of line records.

    Each record has: record (parsed JSON, or None on parse failure),
    raw_record, parse_error (None on success), source_file, line_no.
    Records are accumulated into pages of `buffer_size` before being
    yielded, so downstream transformers receive pages rather than
    individual rows.

    :param items: file items to read
    :type items: Iterator[FileItemDict]
    :param buffer_size: number of line records per yielded page
    :type buffer_size: int
    :yield: pages of line records
    :rtype: Generator[TDataItems, Any, Any]
    """
    page: list[LineRecord] = []
    for file_item in items:
        with file_item.open(mode="rt", compression=compression_mode(file_item["file_name"]), encoding="utf-8") as fh:
            for line_no, line in enumerate(fh, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                page.append(
                    {
                        **parse_jsonl_line(stripped),
                        "source_file": file_item["relative_path"],
                        "line_no": line_no,
                    }
                )
                if len(page) >= buffer_size:
                    yield page
                    page = []
    if page:
        yield page


def route_jsonl_lines(
    items: Iterator[FileItemDict],
    buffer_size: int,
    table_name_of: Callable[[FileItemDict], str],
) -> Generator[TDataItems, Any, Any]:
    """Read each file in items and route each line by its file's table name.

    Valid records yield into the file's table with the given name.
    Records that fail JSON parsing yield into ``{table}_invalid`` with
    raw_record, record (None), source_file, line_no, and parse_error.

    :param items: file items to read
    :type items: Iterator[FileItemDict]
    :param buffer_size: number of rows to buffer per table before yielding
    :type buffer_size: int
    :param table_name_of: function returning the table name for a file item
    :type table_name_of: Callable[[FileItemDict], str]
    :yield: table-tagged rows
    :rtype: Generator[TDataItems, Any, Any]
    """
    for file_item in items:
        buf = DictBuffer(max_items=buffer_size)
        table_name = table_name_of(file_item)
        with file_item.open(mode="rt", compression=compression_mode(file_item["file_name"]), encoding="utf-8") as fh:
            for line_no, line in enumerate(fh, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                line_record = parse_jsonl_line(stripped)
                if line_record["parse_error"] is None:
                    yield from buf.add_item(table_name, line_record["record"])
                else:
                    yield from buf.add_item(
                        f"{table_name}_invalid",
                        {
                            "raw_record": stripped,
                            "record": None,
                            "source_file": file_item["relative_path"],
                            "line_no": line_no,
                            "parse_error": line_record["parse_error"],
                        },
                    )
        yield from buf.flush()
