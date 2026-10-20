"""Unit tests for `cdm_data_loaders.pipelines.xml2db.compaction`.

These build synthetic output files matching the shape `flatten_xml2db_document` produces (see
`cdm_data_loaders.readers.xml`), rather than running a whole pipeline, so they can exercise
compaction's file-rewriting logic directly and cheaply for both loader file formats.
"""

import gzip
import json
from pathlib import Path
from typing import Any, Final

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from xml2db import DataModel

from cdm_data_loaders.core.fields import JSONL, PARQUET
from cdm_data_loaders.pipelines.xml2db.compaction import compact_reused_tables, reused_table_names
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


def _write_table_files(
    dataset_dir: Path, table_name: str, file_rows: list[list[dict[str, Any]]], loader_file_format: str
) -> None:
    """Write one output file per list of rows in `file_rows`, mimicking multiple chunks/files."""
    for index, rows in enumerate(file_rows):
        if loader_file_format == PARQUET:
            _write_parquet(dataset_dir / table_name / f"part_{index}.parquet", rows)
        else:
            _write_jsonl_gz(dataset_dir / table_name / f"part_{index}.jsonl.gz", rows)


def _read_table_rows(dataset_dir: Path, table_name: str, loader_file_format: str) -> list[dict[str, Any]]:
    table_dir = dataset_dir / table_name
    rows: list[dict[str, Any]] = []
    if loader_file_format == PARQUET:
        for file_path in table_dir.glob("*.parquet"):
            rows.extend(pq.read_table(file_path).to_pylist())
    else:
        for file_path in table_dir.glob("*.jsonl*"):
            with gzip.open(file_path, "rt", encoding="utf-8") as fh:
                rows.extend(json.loads(line) for line in fh if line.strip())
    return rows


def test_reused_table_names_pass_excludes_root_table(uniref_model: DataModel) -> None:
    """The virtual root table is never content-addressed, so it must not be a compaction target."""
    names = reused_table_names(uniref_model)
    assert ROOT_SHORT_NAME not in names
    assert {"entry", "property", "representativeMember"} <= set(names)


@pytest.mark.parametrize("loader_file_format", [JSONL, PARQUET])
def test_compact_reused_tables_pass_removes_literal_duplicate_rows(
    uniref_model: DataModel, tmp_path: Path, loader_file_format: str
) -> None:
    """Two files/chunks contributing an overlapping property collapse to the distinct-key count."""
    dataset_dir = tmp_path / "dataset"
    # simulate two chunks: the second repeats "property:aaa" (as flatten_xml2db_document would,
    # since identical content always gets the same content-addressed key).
    _write_table_files(
        dataset_dir,
        "property",
        [
            [
                {"pk_property": "property:aaa", "type": "isSeed", "value": "true"},
                {"pk_property": "property:bbb", "type": "member count", "value": "1"},
            ],
            [
                {"pk_property": "property:aaa", "type": "isSeed", "value": "true"},
                {"pk_property": "property:ccc", "type": "member count", "value": "2"},
            ],
        ],
        loader_file_format,
    )

    removed = compact_reused_tables(dataset_dir, uniref_model, loader_file_format)

    assert removed == {"property": 1}
    remaining = _read_table_rows(dataset_dir, "property", loader_file_format)
    assert len(remaining) == 3  # noqa: PLR2004
    assert {row["pk_property"] for row in remaining} == {"property:aaa", "property:bbb", "property:ccc"}
    # the surviving "aaa" row still carries its original (identical, either copy is fine) content.
    aaa_row = next(row for row in remaining if row["pk_property"] == "property:aaa")
    assert aaa_row == {"pk_property": "property:aaa", "type": "isSeed", "value": "true"}


def test_compact_reused_tables_pass_no_duplicates_is_a_no_op(uniref_model: DataModel, tmp_path: Path) -> None:
    """A table with no literal duplicates is left out of the result and its files untouched."""
    dataset_dir = tmp_path / "dataset"
    _write_table_files(
        dataset_dir,
        "property",
        [[{"pk_property": "property:aaa", "type": "isSeed", "value": "true"}]],
        JSONL,
    )
    original_files = sorted((dataset_dir / "property").iterdir())

    removed = compact_reused_tables(dataset_dir, uniref_model, JSONL)

    assert removed == {}
    assert sorted((dataset_dir / "property").iterdir()) == original_files


def test_compact_reused_tables_pass_missing_table_directory_is_skipped(uniref_model: DataModel, tmp_path: Path) -> None:
    """A reused table that was never written (no output directory at all) is silently skipped."""
    dataset_dir = tmp_path / "dataset"
    dataset_dir.mkdir()

    assert compact_reused_tables(dataset_dir, uniref_model, JSONL) == {}


def test_compact_reused_tables_pass_only_touches_tables_with_duplicates(
    uniref_model: DataModel, tmp_path: Path
) -> None:
    """One table's duplicates are compacted while a sibling table without duplicates is untouched."""
    dataset_dir = tmp_path / "dataset"
    _write_table_files(
        dataset_dir,
        "property",
        [
            [{"pk_property": "property:aaa", "type": "isSeed", "value": "true"}],
            [{"pk_property": "property:aaa", "type": "isSeed", "value": "true"}],
        ],
        JSONL,
    )
    _write_table_files(
        dataset_dir,
        "entry",
        [[{"pk_entry": "entry:xxx", "id": "UniRef50_1", "name": "Cluster: x"}]],
        JSONL,
    )

    removed = compact_reused_tables(dataset_dir, uniref_model, JSONL)

    assert removed == {"property": 1}
    assert len(_read_table_rows(dataset_dir, "entry", JSONL)) == 1
