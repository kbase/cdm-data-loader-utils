"""Shared fixtures for xmltopyarrow pipeline integration tests."""

import gzip
import json
from collections.abc import Callable, Generator
from pathlib import Path
from typing import Any

import duckdb
import pytest
import xmltodict
from dlt.common.pipeline import LoadInfo
from frozendict import frozendict

from cdm_data_loaders.pipelines.xmltopyarrow.pipeline import run_xmltopyarrow_pipeline
from cdm_data_loaders.pipelines.xmltopyarrow.settings import XmlToPyArrowSettings
from cdm_data_loaders.readers.xml import DEFAULT_XMLTODICT_ARGS
from tests.cdm_data_loaders.pipelines.helpers import (
    denormalize_identifier,
    duckdb_reference_from_jsonl,
    make_cts_settings,
    read_jsonl_file,
    strip_metadata_columns,
)
from tests.conftest import N_REFERENCE_XML_ENTRIES, REFERENCE_XML_TAG
from tests.xml_samples import LIBRARY_XSD

DEFAULT_XMLTOPYARROW_SETTINGS: frozendict = frozendict(
    {
        "buffer_size": 10,
        "dlt_dev_mode": False,
        "file_glob": "*.xml*",
        "log_interval": 1000,
        "use_output_dir_for_pipeline_metadata": False,
    }
)


@pytest.fixture
def library_xsd(tmp_path: Path) -> Path:
    """Write LIBRARY_XSD to tmp_path and return its path."""
    xsd_path = tmp_path / "library.xsd"
    xsd_path.write_text(LIBRARY_XSD, encoding="utf-8")
    return xsd_path


@pytest.fixture
def settings_factory(
    tmp_path: Path,
    dlt_destination_config: str,
    library_xsd: Path,
) -> Callable[..., XmlToPyArrowSettings]:
    """Return a factory that builds a valid XmlToPyArrowSettings, with fields open to override."""

    def _factory(**overrides: Any) -> XmlToPyArrowSettings:  # noqa: ANN401
        defaults: dict[str, Any] = {
            **DEFAULT_XMLTOPYARROW_SETTINGS,
            "dataset_name": "xml_test_dataset",
            "table_name": "book",
            "use_destination": dlt_destination_config,
            "xml_tag": "book",
            "xsd_file": str(library_xsd),
        }
        defaults.update(overrides)
        return make_cts_settings(XmlToPyArrowSettings, tmp_path, **defaults)

    return _factory


def import_jsonl_data(jsonl_path: Path) -> list[dict[str, Any]]:
    """Read the reference JSONL into a list of canonical nested entry dicts."""
    return read_jsonl_file(jsonl_path)


def write_reference_xml_as_jsonl(source_xml: Path, jsonl_path: Path) -> int:
    """Parse the UniRef XML with xmltodict and write the reference JSONL."""
    with gzip.open(source_xml, "rb") as fh:
        document = xmltodict.parse(fh.read(), **DEFAULT_XMLTODICT_ARGS)

    entries = document["UniRef50"]["entry"]
    with jsonl_path.open("w", encoding="utf-8") as out:
        for entry in entries:
            out.write(json.dumps({"entry": entry}, ensure_ascii=False) + "\n")
    return len(entries)


@pytest.fixture(scope="session")
def reference_xml_as_jsonl(reference_xml_data_dir: Path) -> Path:
    """Generate the reference JSONL from the chunk_100_el source, if not present."""
    jsonl_path = reference_xml_data_dir / "jsonl" / "uniref_100.jsonl"
    if not jsonl_path.exists():
        source_xml = reference_xml_data_dir / "chunk_100_el" / "uniref_100_part_01.xml.gz"
        n_entries = write_reference_xml_as_jsonl(source_xml, jsonl_path)
        assert n_entries == N_REFERENCE_XML_ENTRIES

    return jsonl_path


@pytest.fixture(scope="session")
def reference_xml_entries(reference_xml_as_jsonl: Path) -> list[dict[str, Any]]:
    """The canonical nested entry dicts read back from the reference JSONL."""
    ref_json = import_jsonl_data(reference_xml_as_jsonl)
    assert isinstance(ref_json, list)
    for line in ref_json:
        assert list(line.keys()) == ["entry"]
        entry = line.get("entry", {})
        for key in ("_id", "name", "property", "representativeMember"):
            assert key in entry
        assert entry.get("_id") is not None
        assert entry.get("_id").startswith("UniRef50_")
    assert len(ref_json) == N_REFERENCE_XML_ENTRIES
    return ref_json


@pytest.fixture(scope="session")
def sorted_reference_xml_entries(reference_xml_entries: list[dict[str, Any]]) -> list[str]:
    """Sorted XSD-normalized entries for the uniref reference dataset."""
    return sorted(
        json.dumps(normalize_uniref_entry_for_xsd(entry), sort_keys=True, default=str)
        for entry in reference_xml_entries
    )


@pytest.fixture(scope="session")
def reference_xml_dataset(reference_xml_as_jsonl: Path) -> Generator[duckdb.DuckDBPyConnection]:
    """Load the reference JSONL into an in-memory DuckDB connection with nested structures."""
    with duckdb_reference_from_jsonl(reference_xml_as_jsonl) as conn:
        yield conn


@pytest.fixture
def run_xmltopyarrow_pipeline_fxt(
    reference_xml_data_dir: Path,
    dlt_destination_config: str,
    isolated_run_factory: Callable[..., Callable[..., tuple[LoadInfo | None, Path]]],
) -> Callable[..., Any]:
    """Return a factory that runs the xmltopyarrow pipeline on a chunk dir with full isolation.

    Each run gets a unique output_dir and dataset_name under tmp_path, so
    parallel-safe repeated runs never collide or accumulate rows.
    """
    base_kwargs = {
        **DEFAULT_XMLTOPYARROW_SETTINGS,
        "settings_cls": XmlToPyArrowSettings,
        "table_name": "entry",
        "use_destination": dlt_destination_config,
        "xml_tag": REFERENCE_XML_TAG,
        "xsd_file": str(reference_xml_data_dir / "uniref.xsd"),
    }
    runner = isolated_run_factory(
        dataset_prefix="reference_run",
        run_fn=run_xmltopyarrow_pipeline,
        base_kwargs=base_kwargs,
    )

    def _factory(chunk_dir: str, **overrides: Any) -> tuple[LoadInfo | None, Path]:  # noqa: ANN401
        return runner(input_dir=str(reference_xml_data_dir / chunk_dir), **overrides)

    return _factory


def _strip_none_values(value: Any) -> Any:  # noqa: ANN401
    """Recursively remove keys with None values from nested dicts and lists."""
    if isinstance(value, dict):
        return {k: _strip_none_values(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [_strip_none_values(item) for item in value]
    return value


def strip_xsd_arrow_row(row: dict[str, Any]) -> dict[str, Any]:
    """Strip dlt metadata columns, unmangle top-level column names, and drop None struct fields."""
    unmangled = {denormalize_identifier(k): v for k, v in strip_metadata_columns(row).items()}
    return _strip_none_values(unmangled)


def _normalize_xsd_reference_member(member: dict[str, Any]) -> dict[str, Any]:
    """Normalize an xmltodict member dict to match XSD-typed cardinality and scalar types."""
    normalized = dict(member)
    db_ref = dict(normalized["dbReference"])
    if "property" in db_ref and isinstance(db_ref["property"], dict):
        db_ref["property"] = [db_ref["property"]]
    normalized["dbReference"] = db_ref
    if "sequence" in normalized and isinstance(normalized["sequence"], dict):
        seq = dict(normalized["sequence"])
        if "_length" in seq:
            seq["_length"] = int(seq["_length"])
        normalized["sequence"] = seq
    return normalized


def normalize_uniref_entry_for_xsd(raw_entry: dict[str, Any]) -> dict[str, Any]:
    """Normalize a raw xmltodict UniRef entry dict to match XSD-backed Arrow output."""
    entry = dict(raw_entry.get("entry", raw_entry))
    if "property" in entry and isinstance(entry["property"], dict):
        entry["property"] = [entry["property"]]
    entry["representativeMember"] = _normalize_xsd_reference_member(entry["representativeMember"])
    if "member" in entry:
        members = entry["member"] if isinstance(entry["member"], list) else [entry["member"]]
        entry["member"] = [_normalize_xsd_reference_member(m) for m in members]
    return entry


def reconstruct_xsd_entries(tables: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    """Extract normalized XSD entry dicts from the single Arrow-backed entry table."""
    return [strip_xsd_arrow_row(row) for row in tables.get("entry", [])]
