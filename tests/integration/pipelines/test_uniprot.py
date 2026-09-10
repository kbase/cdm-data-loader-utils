"""Tests for the UniProt DLT pipeline."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import dlt
import pytest
from frozendict import frozendict

from cdm_data_loaders.pipelines import uniprot_kb as uniprot_module
from cdm_data_loaders.pipelines.uniprot_kb import (
    UNIPROT_LOG_INTERVAL,
    UniProtSettings,
    cli,
    parse_uniprot,
)
from tests.cdm_data_loaders.core.conftest import (
    TEST_BATCH_FILE_SETTINGS,
    TEST_BATCH_FILE_SETTINGS_RECONCILED,
    make_settings_autofill_config,
)
from tests.conftest import TEST_DATA_DIR

# Directory of real UniProt XML fixtures (chunk_00001.xml ... chunk_00004.xml),
# named so they match the NumericFileSequenceBatcher file-sequence regex.
UNIPROT_FIXTURE_DIR = Path("tests") / "data" / "uniprot" / "uniprot_kb" / "chunk_4"


@pytest.fixture
def test_settings() -> UniProtSettings:
    """Provide a minimal valid UniProtSettings object."""
    return make_settings_autofill_config(UniProtSettings)  # type: ignore[reportReturnType]


TEST_SETTINGS = frozendict(
    {**TEST_BATCH_FILE_SETTINGS, "log_interval": UNIPROT_LOG_INTERVAL},
)

TEST_SETTINGS_RECONCILED = frozendict({**TEST_BATCH_FILE_SETTINGS_RECONCILED, "log_interval": UNIPROT_LOG_INTERVAL})


UNIPROT_FIXTURE_DIR = TEST_DATA_DIR / "uniprot" / "uniprot_kb" / "chunk_4"


# Integration tests for the UniProt pipeline
@pytest.fixture
def duckdb_uniprot_settings_args(tmp_path: Path) -> frozendict:
    """Arguments for initialising a UniProtSettings object."""
    output_dir = tmp_path / "output_dir"
    output_dir.mkdir()

    return frozendict(
        {
            "input_dir": str(UNIPROT_FIXTURE_DIR),
            "output_dir": str(output_dir),
            "use_destination": "local_fs",
            "use_output_dir_for_pipeline_metadata": False,
        }
    )


# These use a local duckdb instance to exercise the full UniProt pipeline.
@pytest.fixture
def duckdb_uniprot_settings(duckdb_uniprot_settings_args: frozendict) -> UniProtSettings:
    """Provide UniProtSettings pointing at the real UniProt XML fixtures."""
    return make_settings_autofill_config(UniProtSettings, duckdb_uniprot_settings_args)  # pyright: ignore[reportReturnType]


def _run_uniprot_duckdb_pipeline(settings: UniProtSettings, tmp_path: Path) -> tuple[Any, Any]:
    """Run ``parse_uniprot`` through a real DuckDB pipeline and return (pipeline, load_info)."""
    pipeline = dlt.pipeline(
        pipeline_name="test_uniprot_pipeline",
        destination="duckdb",
        dataset_name="test_uniprot",
        pipelines_dir=str(tmp_path / "pipelines"),
    )
    load_info = pipeline.run(parse_uniprot(settings))
    return pipeline, load_info


def test_integration_uniprot_pipeline_loads_into_duckdb(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
) -> None:
    """Full integration test: parse the real fixture files into a DuckDB destination.

    Confirms that the pipeline runs end-to-end without failed jobs and that the
    expected CDM tables are populated with data queried back out of DuckDB.
    """
    pipeline, load_info = _run_uniprot_duckdb_pipeline(duckdb_uniprot_settings, tmp_path)

    assert not load_info.has_failed_jobs

    # the parser emits an "entity" row per UniProt entry; the fixtures contain
    # several entries across the chunk files, so this table must not be empty.
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT COUNT(*) FROM entity") as cur,
    ):
        (entity_count,) = cur.fetchone()

    assert entity_count > 0

    # every entity must have a "uniprot:"-prefixed entity_id and be a protein.
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT entity_id, entity_type FROM entity") as cur,
    ):
        rows = cur.fetchall()

    assert all(entity_id.startswith("uniprot:") for entity_id, _ in rows)
    assert all(entity_type == "protein" for _, entity_type in rows)


def test_integration_uniprot_pipeline_populates_related_tables(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
) -> None:
    """The pipeline should populate the related CDM tables produced by the parser."""
    pipeline, load_info = _run_uniprot_duckdb_pipeline(duckdb_uniprot_settings, tmp_path)

    assert not load_info.has_failed_jobs

    # tables that the parser always emits at least one row for, per entry
    for table_name in ("entity", "identifier", "name", "protein", "entity_x_source_file"):
        with (
            pipeline.sql_client() as client,
            client.execute_query(f"SELECT COUNT(*) FROM {table_name}") as cur,  # noqa: S608
        ):
            (count,) = cur.fetchone()
        assert count > 0, f"expected rows in table '{table_name}'"

    # entity_x_source_file must reference the fixture files we loaded from
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT DISTINCT source_file FROM entity_x_source_file") as cur,
    ):
        source_files = {row[0] for row in cur.fetchall()}

    assert source_files
    assert all(sf.endswith(".xml") for sf in source_files)


def test_integration_cli_uniprot_pipeline_output_validated(
    duckdb_uniprot_settings: UniProtSettings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the real ``cli()`` wiring end-to-end against a DuckDB destination."""
    monkeypatch.setattr(uniprot_module, "UniProtSettings", MagicMock(return_value=duckdb_uniprot_settings))

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        *,
        settings: UniProtSettings,
        resource: Any,
        pipeline_kwargs: dict[str, Any],
    ) -> None:
        """Replacement for core.run_pipeline that runs the resource through DuckDB."""
        assert settings is duckdb_uniprot_settings
        assert pipeline_kwargs == {"pipeline_name": "uniprot_kb", "dataset_name": "uniprot_kb"}
        pipeline = dlt.pipeline(
            pipeline_name="test_uniprot_cli_pipeline",
            destination="duckdb",
            dataset_name="test_uniprot_cli",
            pipelines_dir=str(tmp_path / "pipelines"),
        )
        captured["load_info"] = pipeline.run(resource)
        captured["pipeline"] = pipeline

    monkeypatch.setattr(uniprot_module, "run_pipeline", fake_run_pipeline)

    cli()

    uniprot_module.UniProtSettings.assert_called_once_with()

    load_info = captured["load_info"]
    pipeline = captured["pipeline"]
    assert not load_info.has_failed_jobs

    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT COUNT(*) FROM entity") as cur,
    ):
        (entity_count,) = cur.fetchone()

    assert entity_count > 0
