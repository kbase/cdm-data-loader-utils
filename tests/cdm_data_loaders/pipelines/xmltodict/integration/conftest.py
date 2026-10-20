"""Shared fixtures for xml pipeline integration tests."""

import gzip
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import duckdb
import pytest
import xmltodict
from dlt.common.pipeline import LoadInfo
from frozendict import frozendict

from cdm_data_loaders.pipelines.xmltodict.pipeline import run_xmltodict_pipeline
from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings
from cdm_data_loaders.readers.xml import DEFAULT_XMLTODICT_ARGS
from tests.conftest import N_REFERENCE_XML_ENTRIES, REFERENCE_XML_FIXTURE_DIR, REFERENCE_XML_TAG

DEFAULT_XMLTODICT_SETTINGS: frozendict = frozendict(
    {
        "buffer_size": 10,
        "dlt_dev_mode": False,
        "file_glob": "*.xml*",
        "log_interval": 1000,
        "use_output_dir_for_pipeline_metadata": False,
    }
)


def make_settings(tmp_path: Path, dlt_destination_config: str, **overrides: Any) -> dict[str, Any]:
    """Generate a dict of settings for the XmlToDict ingestion route."""
    log_config_file = overrides.get("log_config_file", tmp_path / "logging.conf")
    log_config_file.touch()
    output_dir = overrides.get("output_dir", tmp_path / "output")
    output_dir.mkdir(exist_ok=True)
    input_dir = overrides.get("input_dir", tmp_path / "input")

    kwargs: dict[str, Any] = {
        **DEFAULT_XMLTODICT_SETTINGS,
        "dataset_name": "xml_test_dataset",
        "input_dir": str(input_dir),
        "log_config_file": str(log_config_file),
        "output_dir": str(output_dir),
        "table_name": "book",
        "use_destination": dlt_destination_config,
        "xml_tag": "book",
    }
    kwargs.update(overrides)
    return kwargs


@pytest.fixture
def settings_factory(tmp_path: Path, dlt_destination_config: str) -> Callable[..., XmlToDictSettings]:
    """Return a factory that builds a valid XmlToDictSettings, with fields open to override."""

    def _factory(**overrides: Any) -> XmlToDictSettings:
        settings_dict = make_settings(tmp_path, dlt_destination_config, **overrides)
        # ensure that the input dir exists
        input_dir = Path(settings_dict["input_dir"])
        input_dir.mkdir(exist_ok=True)

        return XmlToDictSettings(**settings_dict)

    return _factory


@pytest.fixture
def scenario_input_dir(test_data_dir: Path) -> Callable[[str], str]:
    """Return a function that maps a scenario name to its XML data directory."""

    def _resolve(scenario: str) -> str:
        path = test_data_dir / "xml" / scenario
        assert path.is_dir(), f"missing test fixture directory: {path}"
        return str(path)

    return _resolve


@pytest.fixture(scope="session")
def reference_xml_data_dir() -> Path:
    """Path to the uniref fixture directory."""
    return REFERENCE_XML_FIXTURE_DIR


def import_jsonl_data(jsonl_path: Path) -> list[dict[str, Any]]:
    """Read the reference JSONL into a list of canonical nested entry dicts."""
    with jsonl_path.open(encoding="utf-8") as fh:
        return [json.loads(line) for line in fh if line.strip()]


def write_reference_xml_as_jsonl(source_xml: Path, jsonl_path: Path) -> int:
    """Parse the UniRef XML with xmltodict and write the reference JSONL.

    The top-level element is ignored; one JSON object per entry element is written,
    wrapped as {"entry": ...} so that nested content is preserved verbatim.
    """
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
    """Sorted entries for the uniref reference dataset."""
    return sorted(json.dumps(entry.get("entry"), sort_keys=True, default=str) for entry in reference_xml_entries)


@pytest.fixture(scope="session")
def reference_xml_dataset(reference_xml_as_jsonl: Path) -> duckdb.DuckDBPyConnection:
    """Load the reference JSONL into an in-memory DuckDB connection with nested structures.

    read_json_auto infers native STRUCT/STRUCT[] types for the nested content;
    `member` infers as JSON because its shape varies across entries.
    """
    connection = duckdb.connect()
    jsonl_path = str(reference_xml_as_jsonl)
    connection.execute(
        "CREATE TABLE reference AS SELECT * FROM read_json_auto(?, sample_size=200)",
        [jsonl_path],
    )
    return connection


@pytest.fixture
def run_xmltodict_pipeline_fxt(
    reference_xml_data_dir: Path,
    dlt_destination_config: str,
    isolated_run_factory: Callable[..., Callable[..., tuple[LoadInfo | None, Path]]],
) -> Callable[..., Any]:
    """Return a factory that runs the xmltodict pipeline on a chunk dir with full isolation.

    Each run gets a unique output_dir and dataset_name under tmp_path, so
    parallel-safe repeated runs never collide or accumulate rows.
    """
    base_kwargs = {
        **DEFAULT_XMLTODICT_SETTINGS,
        "settings_cls": XmlToDictSettings,
        "table_name": "entry",
        "use_destination": dlt_destination_config,
        "xml_tag": REFERENCE_XML_TAG,
    }
    runner = isolated_run_factory(
        dataset_prefix="reference_run",
        run_fn=run_xmltodict_pipeline,
        base_kwargs=base_kwargs,
    )

    def _factory(chunk_dir: str, **overrides: Any) -> tuple[LoadInfo | None, Path]:
        return runner(input_dir=str(reference_xml_data_dir / chunk_dir), **overrides)

    return _factory
