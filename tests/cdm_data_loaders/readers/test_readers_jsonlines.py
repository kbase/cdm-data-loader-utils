"""Tests for the JSONL line-reading helpers in readers.jsonlines."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from dlt.sources.filesystem import filesystem

from cdm_data_loaders.readers.jsonlines import (
    compression_mode,
    parse_jsonl_line,
    route_jsonl_lines,
    stream_jsonl_lines,
)


@pytest.mark.parametrize(
    ("file_name", "expected"),
    [pytest.param("data.jsonl.gz", "enable", id="gzip"), pytest.param("data.jsonl", "disable", id="plain")],
)
def test_compression_mode_pass_returns_mode_for_suffix(file_name: str, expected: str) -> None:
    """Gzip-suffixed names enable compression; others disable it."""
    assert compression_mode(file_name) == expected


def test_parse_jsonl_line_pass_valid_json_returns_record() -> None:
    """A valid JSON line parses into the record with no error."""
    result = parse_jsonl_line('{"widget_id": "a", "count": 1}')

    assert result == {
        "record": {"widget_id": "a", "count": 1},
        "raw_record": '{"widget_id": "a", "count": 1}',
        "parse_error": None,
    }


def test_parse_jsonl_line_fail_malformed_json_captures_error() -> None:
    """A line that is not valid JSON sets record to None and fills parse_error."""
    result = parse_jsonl_line("{not valid json")

    assert result["record"] is None
    assert result["raw_record"] == "{not valid json"
    assert result["parse_error"] is not None


def read_items(directory: Path, file_glob: str = "*.jsonl*") -> list[Any]:
    """List the filesystem items in a directory."""
    return list(filesystem(bucket_url=str(directory), file_glob=file_glob))


def read_records(directory: Path, buffer_size: int = 10) -> list[dict[str, Any]]:
    """Stream every line record in a directory, flattened across pages."""
    return [
        record for page in stream_jsonl_lines(iter(read_items(directory)), buffer_size=buffer_size) for record in page
    ]


def test_stream_jsonl_lines_pass_yields_pages_of_line_records(tmp_path: Path) -> None:
    """stream_jsonl_lines yields pages of at most buffer_size records with provenance."""
    directory = tmp_path / "widget"
    directory.mkdir()
    lines = [f'{{"widget_id": "w{i}", "count": {i}}}' for i in range(3)]
    (directory / "batch.jsonl").write_text("\n".join(lines) + "\n", encoding="utf-8")

    pages = list(stream_jsonl_lines(iter(read_items(directory)), buffer_size=2))

    assert [len(page) for page in pages] == [2, 1]
    first = pages[0][0]
    assert first["record"] == {"widget_id": "w0", "count": 0}
    assert first["parse_error"] is None
    assert first["source_file"] == "batch.jsonl"
    assert first["line_no"] == 1


def test_stream_jsonl_lines_pass_blank_lines_are_skipped(tmp_path: Path) -> None:
    """Blank and whitespace-only lines produce no output records."""
    directory = tmp_path / "widget"
    directory.mkdir()
    (directory / "blank.jsonl").write_text('{"widget_id": "a"}\n\n   \n{"widget_id": "b"}\n', encoding="utf-8")

    records = read_records(directory)

    assert [record["record"]["widget_id"] for record in records] == ["a", "b"]


def test_stream_jsonl_lines_fail_malformed_json_line_captured_not_raised(tmp_path: Path) -> None:
    """A line that is not valid JSON sets record to None and fills parse_error; it does not raise."""
    directory = tmp_path / "widget"
    directory.mkdir()
    (directory / "bad.jsonl").write_text("{not valid json\n", encoding="utf-8")

    records = read_records(directory)

    assert len(records) == 1
    assert records[0]["record"] is None
    assert records[0]["parse_error"] is not None
    assert records[0]["raw_record"] == "{not valid json"
    assert records[0]["line_no"] == 1


def test_stream_jsonl_lines_pass_empty_file_yields_nothing(tmp_path: Path) -> None:
    """An empty file produces no output records."""
    directory = tmp_path / "widget"
    directory.mkdir()
    (directory / "empty.jsonl").write_text("", encoding="utf-8")

    assert read_records(directory) == []


def test_stream_jsonl_lines_pass_line_numbers_restart_per_file(tmp_path: Path) -> None:
    """Line numbers restart at 1 in each file; source_file names the file."""
    directory = tmp_path / "widget"
    directory.mkdir()
    (directory / "first.jsonl").write_text('{"widget_id": "a"}\n{"widget_id": "b"}\n', encoding="utf-8")
    (directory / "second.jsonl").write_text('{"widget_id": "c"}\n', encoding="utf-8")

    records = read_records(directory)

    by_source: dict[str, list[int]] = {}
    for record in records:
        by_source.setdefault(Path(record["source_file"]).name, []).append(record["line_no"])

    assert by_source == {"first.jsonl": [1, 2], "second.jsonl": [1]}


def test_stream_jsonl_lines_pass_gzip_file_is_decompressed(
    tmp_path: Path, write_gzip_file: Callable[[Path, str, str], Path]
) -> None:
    """A .gz-suffixed file is decompressed before its lines are read."""
    directory = tmp_path / "widget"
    write_gzip_file(directory, "compressed.jsonl.gz", '{"widget_id": "a", "count": 1}\n')

    records = read_records(directory)

    assert len(records) == 1
    assert records[0]["record"] == {"widget_id": "a", "count": 1}
    assert records[0]["parse_error"] is None


def test_route_jsonl_lines_pass_routes_valid_and_invalid_rows(tmp_path: Path) -> None:
    """Valid lines route to the table; parse failures route to <table>_invalid."""
    directory = tmp_path / "widget"
    directory.mkdir()
    (directory / "data.jsonl").write_text('{"widget_id": "a"}\n{broken\n{"widget_id": "b"}\n', encoding="utf-8")

    pages = list(
        route_jsonl_lines(
            iter(read_items(directory)),
            buffer_size=10,
            table_name_of=lambda file_item: Path(file_item["file_url"]).parent.name,
        )
    )

    by_table: dict[str, list[dict[str, Any]]] = {}
    for page in pages:
        by_table.setdefault(page.meta.table_name, []).extend(page.data)

    assert list(by_table["widget"]) == [{"widget_id": "a"}, {"widget_id": "b"}]
    invalid = by_table["widget_invalid"]
    assert len(invalid) == 1
    assert invalid[0]["record"] is None
    assert invalid[0]["raw_record"] == "{broken"
    assert invalid[0]["line_no"] == 2
    assert invalid[0]["source_file"] == "data.jsonl"
    assert invalid[0]["parse_error"] is not None
