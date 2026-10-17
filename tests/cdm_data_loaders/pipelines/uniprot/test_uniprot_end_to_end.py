"""Tests for the UniProt DLT pipeline.

The fixture directory contains 15 entries across four chunk files; exact
per-table counts below are derived from those fixtures.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.core.fields import LOCAL_FS
from cdm_data_loaders.pipelines import uniprot_kb as uniprot_module
from cdm_data_loaders.pipelines.uniprot_kb import (
    UniProtSettings,
    cli,
    parse_uniprot,
)
from tests.cdm_data_loaders.pipelines.conftest import duckdb_pipeline
from tests.conftest import TEST_DATA_DIR

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
        use_output_dir_for_pipeline_metadata=False,
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


def _fetch_all(pipeline: Any, query: str) -> list[tuple[Any, ...]]:
    """Run a query against the pipeline's DuckDB destination and return all rows."""
    with (
        pipeline.sql_client() as client,
        client.execute_query(query) as cur,
    ):
        return list(cur.fetchall())


def test_integration_uniprot_pipeline_loads_into_duckdb(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
) -> None:
    """Full integration test: parse the real fixture files into a DuckDB destination.

    Confirms that the pipeline runs end-to-end without failed jobs and that the
    expected CDM tables are populated with the exact expected row counts.
    """
    pipeline, load_info = _run_uniprot_duckdb_pipeline(duckdb_uniprot_settings, tmp_path, "loads")

    assert not load_info.has_failed_jobs

    for table_name, expected_count in EXPECTED_ROW_COUNTS.items():
        ((count,),) = _fetch_all(pipeline, f"SELECT COUNT(*) FROM {table_name}")  # noqa: S608
        assert count == expected_count, f"expected {expected_count} rows in '{table_name}', got {count}"


def test_integration_uniprot_pipeline_entity_and_source_file_content(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
) -> None:
    """Entity rows are uniprot proteins and source-file rows reference every fixture file."""
    pipeline, load_info = _run_uniprot_duckdb_pipeline(duckdb_uniprot_settings, tmp_path, "content")

    assert not load_info.has_failed_jobs

    # every entity must have a "uniprot:"-prefixed entity_id and be a protein.
    rows = _fetch_all(pipeline, "SELECT entity_id, entity_type FROM entity")
    assert len(rows) == UNIPROT_ENTRY_COUNT
    assert all(entity_id.startswith("uniprot:") for entity_id, _ in rows)
    assert all(entity_type == "protein" for _, entity_type in rows)

    # entity_x_source_file must reference exactly the fixture files we loaded from
    source_files = {row[0] for row in _fetch_all(pipeline, "SELECT DISTINCT source_file FROM entity_x_source_file")}
    assert source_files == {str(UNIPROT_FIXTURE_DIR / name) for name in UNIPROT_FIXTURE_FILES}


def test_integration_cli_uniprot_pipeline_output_validated(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the real ``cli.run()`` wiring end-to-end against a DuckDB destination."""
    captured: dict[str, Any] = {}

    def fake_run_cli(
        settings_cls: type[UniProtSettings],
        pipeline_fn: Callable[..., Any],
        cli_args: list[str] | None = None,
        settings_kwargs: dict[str, Any] | None = None,
    ) -> Any:
        """Replacement for core.run_cli that dispatches to pipeline_fn with the fixture settings."""
        assert settings_cls is UniProtSettings
        return pipeline_fn(duckdb_uniprot_settings)

    def fake_run_pipeline(
        *,
        settings: UniProtSettings,
        resource: Any,
        pipeline_kwargs: dict[str, Any],
    ) -> None:
        """Replacement for core.run_pipeline that runs the resource through DuckDB."""
        assert settings is duckdb_uniprot_settings
        assert pipeline_kwargs == {"pipeline_name": "uniprot_kb", "dataset_name": "uniprot_kb"}
        pipeline = duckdb_pipeline(tmp_path, "test_uniprot_cli_pipeline", "test_uniprot_cli")
        captured["load_info"] = pipeline.run(resource)
        captured["pipeline"] = pipeline

    monkeypatch.setattr("cdm_data_loaders.pipelines.core.run_cli", fake_run_cli)
    monkeypatch.setattr(uniprot_module, "run_pipeline", fake_run_pipeline)

    cli.run()

    load_info = captured["load_info"]
    pipeline = captured["pipeline"]
    assert not load_info.has_failed_jobs

    ((entity_count,),) = _fetch_all(pipeline, "SELECT COUNT(*) FROM entity")
    assert entity_count == UNIPROT_ENTRY_COUNT
