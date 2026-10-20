"""Tests for the UniRef DLT pipeline.

The fixture file contains three clusters; exact per-table counts below are
derived from that fixture. Every test runs once per valid UniRef variant.
"""

from pathlib import Path
from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.pipelines.uniref import (
    UNIREF_VARIANTS,
    UnirefSettings,
    parse_uniref,
)
from tests.cdm_data_loaders.pipelines.conftest import duckdb_pipeline, read_destination_rows
from tests.cdm_data_loaders.pipelines.framework import PipelineTestRunner
from tests.cdm_data_loaders.pipelines.helpers import assert_pipeline_report_written
from tests.conftest import LOCAL_FS, TEST_DATA_DIR

TEST_DEFAULT_UNIREF_VARIANT = "50"

UNIREF_FIXTURE_DIR = TEST_DATA_DIR / "uniprot" / "uniref" / "integration"
UNIREF_FIXTURE_FILE = "uniref_chunk_00001.xml"
UNIREF_CLUSTER_COUNT = 3

EXPECTED_ROW_COUNTS = frozendict(
    {
        "entity": UNIREF_CLUSTER_COUNT,
        "cluster": UNIREF_CLUSTER_COUNT,
        "clustermember": 4,
        "entity_x_source_file": UNIREF_CLUSTER_COUNT,
    }
)


# Integration test for the UniRef pipeline
# Data output goes to a local DuckDB destination.
@pytest.fixture(params=UNIREF_VARIANTS)
def duckdb_uniref_settings(tmp_path: Path, request: pytest.FixtureRequest) -> UnirefSettings:
    """Provide UnirefSettings pointing at the real UniRef XML fixtures, one per variant."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return UnirefSettings(
        variant=request.param,
        input_dir=str(UNIREF_FIXTURE_DIR),
        output_dir=str(output_dir),
        use_destination=LOCAL_FS,
        use_output_dir_for_pipeline_metadata=True,
    )  # pyright: ignore[reportCallIssue]


def _run_uniref_duckdb_pipeline(settings: UnirefSettings, tmp_path: Path, name: str) -> tuple[Any, Any]:
    """Run ``parse_uniref`` through a real DuckDB pipeline and return (pipeline, load_info).

    dlt writes the duckdb destination file into the working directory and appends
    on re-runs, so each test gets its own database file and pipeline name under
    tmp_path to stay isolated from other tests and parametrized cases.
    """
    pipeline = duckdb_pipeline(tmp_path, f"test_uniref_pipeline_{name}", "test_uniref")
    load_info = pipeline.run(parse_uniref(settings))
    return pipeline, load_info


def _fetch_all(pipeline: Any, query: str) -> list[tuple[Any, ...]]:
    """Run a query against the pipeline's DuckDB destination and return all rows."""
    with (
        pipeline.sql_client() as client,
        client.execute_query(query) as cur,
    ):
        return list(cur.fetchall())


def test_integration_uniref_pipeline_loads_into_duckdb(
    duckdb_uniref_settings: UnirefSettings,
    tmp_path: Path,
) -> None:
    """Full integration test: parse the real fixture files into a DuckDB destination.

    Confirms that the pipeline runs end-to-end without failed jobs and that the
    expected CDM tables are populated with the exact expected row counts.
    """
    pipeline, load_info = _run_uniref_duckdb_pipeline(duckdb_uniref_settings, tmp_path, "loads")

    assert not load_info.has_failed_jobs

    for table_name, expected_count in EXPECTED_ROW_COUNTS.items():
        ((count,),) = _fetch_all(pipeline, f"SELECT COUNT(*) FROM {table_name}")  # noqa: S608
        assert count == expected_count, f"expected {expected_count} rows in '{table_name}', got {count}"


def test_integration_uniref_pipeline_entity_and_cluster_content(
    duckdb_uniref_settings: UnirefSettings,
    tmp_path: Path,
) -> None:
    """Entity rows are UniRef clusters and every cluster carries the injected protocol label."""
    pipeline, load_info = _run_uniref_duckdb_pipeline(duckdb_uniref_settings, tmp_path, "content")

    assert not load_info.has_failed_jobs

    # every entity must have a "uniref:"-prefixed entity_id and be a Cluster.
    rows = _fetch_all(pipeline, "SELECT entity_id, entity_type FROM entity")
    assert len(rows) == UNIREF_CLUSTER_COUNT
    assert all(entity_id.startswith("uniref:") for entity_id, _ in rows)
    assert all(entity_type == "Cluster" for _, entity_type in rows)

    # every cluster must carry the injected "UniRef <variant>" protocol label
    protocols = {row[0] for row in _fetch_all(pipeline, "SELECT DISTINCT protocol FROM cluster")}
    assert protocols == {f"UniRef {duckdb_uniref_settings.variant}"}

    # entity_x_source_file must reference exactly the fixture file we loaded from
    source_files = {row[0] for row in _fetch_all(pipeline, "SELECT DISTINCT source_file FROM entity_x_source_file")}
    assert source_files == {str(UNIREF_FIXTURE_DIR.resolve() / UNIREF_FIXTURE_FILE)}


def test_integration_cli_uniref_pipeline_output_validated(
    duckdb_uniref_settings: UnirefSettings,
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Exercise the UniRef pipeline end-to-end across core, cli, and docker execution modes."""
    load_info = pipeline_runner.run("uniref", settings=duckdb_uniref_settings)

    assert load_info is not None
    assert not load_info.has_failed_jobs
    assert load_info.pipeline.pipeline_name == "uniref"
    assert load_info.dataset_name == f"uniref_{duckdb_uniref_settings.variant}"
    assert_pipeline_report_written(Path(duckdb_uniref_settings.output_dir), load_info, metadata_moved=True)

    dataset_dir = Path(duckdb_uniref_settings.output_dir) / load_info.dataset_name
    for table_name, expected_count in EXPECTED_ROW_COUNTS.items():
        rows = read_destination_rows(dataset_dir, table_name)
        assert len(rows) == expected_count, f"expected {expected_count} rows in '{table_name}', got {len(rows)}"
