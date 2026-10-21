"""Tests for the UniProt DLT pipeline.

The fixture directory contains 15 entries across four chunk files; exact
per-table counts below are derived from those fixtures.
"""

from pathlib import Path
from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.pipelines.uniprot_kb_cdm import (
    PIPELINE_NAME,
    UniProtSettings,
    parse_uniprot,
)
from tests.cdm_data_loaders.pipelines.conftest import duckdb_pipeline, read_destination_rows
from tests.cdm_data_loaders.pipelines.framework import PipelineTestRunner
from tests.cdm_data_loaders.pipelines.helpers import (
    assert_pipeline_report_written,
    assert_run_ok,
    fetch_all_duckdb,
)
from tests.conftest import LOCAL_FS, TEST_DATA_DIR

UNIPROT_FIXTURE_DIR = TEST_DATA_DIR / "uniprot" / "uniprot_kb" / "chunk_4"
UNIPROT_FIXTURE_FILES = ("chunk_00001.xml", "chunk_00002.xml", "chunk_00003.xml", "chunk_00004.xml")
UNIPROT_ENTRY_COUNT = 15

# Row counts loaded from the fixtures by a single pipeline run, per table.
EXPECTED_ROW_COUNTS = frozendict(
    {
        "entity": UNIPROT_ENTRY_COUNT,
        "identifier": 696,
        "name": 62,
        "protein": UNIPROT_ENTRY_COUNT,
        "entity_x_publication": 21,
        "entity_x_source_file": UNIPROT_ENTRY_COUNT,
    }
)


# Integration tests for the UniProt pipeline
@pytest.fixture
def duckdb_uniprot_settings(tmp_path: Path) -> UniProtSettings:
    """Provide UniProtSettings pointing at the real UniProt XML fixtures."""
    output_dir = tmp_path / "output_dir"
    output_dir.mkdir()

    return UniProtSettings(
        input_dir=str(UNIPROT_FIXTURE_DIR),
        output_dir=str(output_dir),
        use_destination=LOCAL_FS,
        use_output_dir_for_pipeline_metadata=True,
    )  # pyright: ignore[reportCallIssue]


def _run_uniprot_duckdb_pipeline(settings: UniProtSettings, tmp_path: Path, name: str) -> tuple[Any, Any]:
    """Run ``parse_uniprot`` through a real DuckDB pipeline and return (pipeline, load_info).

    dlt writes the duckdb destination file into the working directory and appends
    on re-runs, so each test gets its own database file and pipeline name under
    tmp_path to stay isolated from other tests.
    """
    pipeline = duckdb_pipeline(tmp_path, f"test_uniprot_pipeline_{name}", "test_uniprot")
    load_info = pipeline.run(parse_uniprot(settings))
    return pipeline, load_info


def test_integration_uniprot_pipeline_loads_into_duckdb(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
) -> None:
    """Full integration test: parse the real fixture files into a DuckDB destination.

    Confirms that the pipeline runs end-to-end without failed jobs and that the
    expected CDM tables are populated with the exact expected row counts.
    """
    pipeline, load_info = _run_uniprot_duckdb_pipeline(duckdb_uniprot_settings, tmp_path, "loads")

    assert_run_ok(load_info)

    for table_name, expected_count in EXPECTED_ROW_COUNTS.items():
        ((count,),) = fetch_all_duckdb(pipeline, f"SELECT COUNT(*) FROM {table_name}")  # noqa: S608
        assert count == expected_count, f"expected {expected_count} rows in '{table_name}', got {count}"


def test_integration_uniprot_pipeline_entity_and_source_file_content(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
) -> None:
    """Entity rows are uniprot proteins and source-file rows reference every fixture file."""
    pipeline, load_info = _run_uniprot_duckdb_pipeline(duckdb_uniprot_settings, tmp_path, "content")

    assert_run_ok(load_info)

    # every entity must have a "uniprot:"-prefixed entity_id and be a protein.
    rows = fetch_all_duckdb(pipeline, "SELECT entity_id, entity_type FROM entity")
    assert len(rows) == UNIPROT_ENTRY_COUNT
    assert all(entity_id.startswith("uniprot:") for entity_id, _ in rows)
    assert all(entity_type == "protein" for _, entity_type in rows)

    # entity_x_source_file must reference exactly the fixture files we loaded from
    source_files = {
        row[0] for row in fetch_all_duckdb(pipeline, "SELECT DISTINCT source_file FROM entity_x_source_file")
    }
    assert source_files == {str(UNIPROT_FIXTURE_DIR.resolve() / name) for name in UNIPROT_FIXTURE_FILES}


def test_integration_cli_uniprot_pipeline_output_validated(
    duckdb_uniprot_settings: UniProtSettings,
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Exercise the UniProt pipeline end-to-end across core, cli, and docker execution modes."""
    load_info = assert_run_ok(pipeline_runner.run(PIPELINE_NAME, settings=duckdb_uniprot_settings))
    assert load_info.pipeline.pipeline_name == PIPELINE_NAME
    assert load_info.dataset_name == PIPELINE_NAME
    assert_pipeline_report_written(Path(duckdb_uniprot_settings.output_dir), load_info, metadata_moved=True)

    dataset_dir = Path(duckdb_uniprot_settings.output_dir) / load_info.dataset_name
    for table_name, expected_count in EXPECTED_ROW_COUNTS.items():
        rows = read_destination_rows(dataset_dir, table_name)
        assert len(rows) == expected_count, f"expected {expected_count} rows in '{table_name}', got {len(rows)}"
