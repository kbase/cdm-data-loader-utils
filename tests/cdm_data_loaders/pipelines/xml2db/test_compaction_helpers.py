"""Unit tests for the private helpers and edge cases in `cdm_data_loaders.pipelines.xml2db.compaction`."""

import gzip
import json
from pathlib import Path
from typing import Any, Final

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from xml2db import DataModel

from cdm_data_loaders.core.fields import JSONL, PARQUET
from cdm_data_loaders.pipelines.xml2db.compaction import (
    _read_expr,
    _source_files,
    _write_compacted,
    compact_reused_tables,
)
from cdm_data_loaders.readers.xml2db_doc import build_xml2db_model

REFERENCE_XSD: Final[Path] = Path("tests") / "data" / "uniprot" / "uniref" / "uniref.xsd"
ROOT_SHORT_NAME: Final[str] = "uniref_test"


@pytest.fixture(scope="module")
def uniref_model() -> DataModel:
    """Build the xml2db DataModel from uniref.xsd once for the whole test module."""
    return build_xml2db_model(REFERENCE_XSD, short_name=ROOT_SHORT_NAME)


def _write_jsonl_gz(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        for row in rows:
            fh.write(json.dumps(row) + "\n")


def _write_parquet(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


# _source_files


def test_source_files_pass_selects_parquet_only(tmp_path: Path) -> None:
    """Only .parquet files are selected when the loader format is parquet."""
    table_dir = tmp_path / "property"
    _write_parquet(table_dir / "b.parquet", [])
    _write_parquet(table_dir / "a.parquet", [])
    (table_dir / "ignored.jsonl.gz").write_text("not a parquet file", encoding="utf-8")
    (table_dir / "ignored.txt").write_text("neither is this", encoding="utf-8")

    files = _source_files(table_dir, PARQUET)

    assert [f.name for f in files] == ["a.parquet", "b.parquet"]


def test_source_files_pass_selects_all_json_variants(tmp_path: Path) -> None:
    """Every .json* variant is selected when the loader format is jsonl, in stable order."""
    table_dir = tmp_path / "property"
    _write_jsonl_gz(table_dir / "b.jsonl.gz", [])
    _write_jsonl_gz(table_dir / "a.jsonl.gz", [])
    _write_jsonl_gz(table_dir / "c.json", [])
    (table_dir / "ignored.parquet").write_bytes(b"not a json file")

    files = _source_files(table_dir, JSONL)

    assert [f.name for f in files] == ["a.jsonl.gz", "b.jsonl.gz", "c.json"]


def test_source_files_pass_empty_directory(tmp_path: Path) -> None:
    """A table directory with no matching files produces an empty list."""
    table_dir = tmp_path / "property"
    table_dir.mkdir()

    assert _source_files(table_dir, PARQUET) == []
    assert _source_files(table_dir, JSONL) == []


def test_source_files_pass_subdirectories_are_ignored(tmp_path: Path) -> None:
    """Subdirectories matching the glob pattern are not selected."""
    table_dir = tmp_path / "property"
    (table_dir / "nested.parquet").mkdir(parents=True)

    assert _source_files(table_dir, PARQUET) == []


# _read_expr


def test_read_expr_pass_parquet_bracketed_file_list(tmp_path: Path) -> None:
    """The parquet read expression reads every file via a bracketed list."""
    source_files = [tmp_path / "a.parquet", tmp_path / "b.parquet"]

    expr = _read_expr(source_files, PARQUET)

    assert expr == f"read_parquet(['{source_files[0].as_posix()}', '{source_files[1].as_posix()}'])"


def test_read_expr_pass_jsonl_uses_read_json_auto(tmp_path: Path) -> None:
    """The jsonl read expression uses read_json_auto."""
    expr = _read_expr([tmp_path / "a.jsonl.gz"], JSONL)

    assert expr == f"read_json_auto(['{tmp_path / 'a.jsonl.gz'}'])" or expr.startswith("read_json_auto")


def test_read_expr_pass_quote_escaping(tmp_path: Path) -> None:
    """A single quote in a path is escaped, so the SQL expression stays valid."""
    tricky_dir = tmp_path / "it's a dir"
    tricky_dir.mkdir()
    source_files = [tricky_dir / "a.parquet"]

    expr = _read_expr(source_files, PARQUET)

    assert "it''s" in expr


# _write_compacted


def test_write_compacted_pass_parquet(tmp_path: Path) -> None:
    """A deduplicated pyarrow table is written as a single compacted.parquet file."""
    table_dir = tmp_path / "property"
    table_dir.mkdir()
    deduped = pa.Table.from_pylist([{"pk_property": "property:aaa", "value": "1"}])

    _write_compacted(deduped, table_dir, PARQUET)

    written = pq.read_table(table_dir / "compacted.parquet").to_pylist()
    assert written == [{"pk_property": "property:aaa", "value": "1"}]


def test_write_compacted_pass_jsonl_gz(tmp_path: Path) -> None:
    """A deduplicated pyarrow table is written as a single compacted.jsonl.gz file."""
    table_dir = tmp_path / "property"
    table_dir.mkdir()
    deduped = pa.Table.from_pylist([{"pk_property": "property:aaa", "value": "1"}])

    _write_compacted(deduped, table_dir, JSONL)

    with gzip.open(table_dir / "compacted.jsonl.gz", "rt", encoding="utf-8") as fh:
        rows = [json.loads(line) for line in fh if line.strip()]
    assert rows == [{"pk_property": "property:aaa", "value": "1"}]


# compact_reused_tables, remaining edge cases


def test_compact_reused_tables_pass_dedupes_across_mixed_file_counts(uniref_model: DataModel, tmp_path: Path) -> None:
    """Duplicates spread across more than two files collapse to the distinct-key count."""
    dataset_dir = tmp_path / "dataset"
    rows = [[{"pk_property": f"property:{n:03d}", "value": str(n)}] for n in range(4)]
    # three more copies of the first key, in separate files
    rows.extend([[{"pk_property": "property:000", "value": "0"}]] * 3)

    table_dir = dataset_dir / "property"
    for index, file_rows in enumerate(rows):
        _write_parquet(table_dir / f"part_{index}.parquet", file_rows)

    removed = compact_reused_tables(dataset_dir, uniref_model, PARQUET)

    assert removed == {"property": 3}

    remaining = pq.read_table(table_dir / "compacted.parquet").to_pylist()
    assert len(remaining) == 4  # noqa: PLR2004 -- one row per distinct key
    assert {row["pk_property"] for row in remaining} == {f"property:{n:03d}" for n in range(4)}


def test_compact_reused_tables_pass_jsonl_duplicates_across_files(uniref_model: DataModel, tmp_path: Path) -> None:
    """The jsonl loader format deduplicates across multiple .jsonl.gz files."""
    dataset_dir = tmp_path / "dataset"
    _write_jsonl_gz(dataset_dir / "property" / "a.jsonl.gz", [{"pk_property": "property:aaa", "value": "1"}])
    _write_jsonl_gz(dataset_dir / "property" / "b.jsonl.gz", [{"pk_property": "property:aaa", "value": "1"}])

    removed = compact_reused_tables(dataset_dir, uniref_model, JSONL)

    assert removed == {"property": 1}
    table_dir = dataset_dir / "property"
    remaining_files = sorted(table_dir.iterdir())
    assert [f.name for f in remaining_files] == ["compacted.jsonl.gz"]


def test_compact_reused_tables_pass_empty_table_directory_is_skipped(uniref_model: DataModel, tmp_path: Path) -> None:
    """A table directory that exists but holds no matching files is skipped, not an error."""
    dataset_dir = tmp_path / "dataset"
    (dataset_dir / "property").mkdir(parents=True)
    (dataset_dir / "property" / "notes.txt").write_text("not an output file", encoding="utf-8")

    assert compact_reused_tables(dataset_dir, uniref_model, PARQUET) == {}
