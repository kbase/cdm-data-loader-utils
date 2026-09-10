"""Tests for the UniRef DLT pipeline."""

from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import dlt
import pytest
from frozendict import frozendict

from cdm_data_loaders.parsers.uniprot.uniref import UNIREF_VARIANTS
from cdm_data_loaders.pipelines import uniref as uniref_module
from cdm_data_loaders.pipelines.uniref import (
    VARIANT,
    UnirefSettings,
    cli,
    parse_uniref,
)
from tests.cdm_data_loaders.core.conftest import (
    TEST_BATCH_FILE_SETTINGS,
    TEST_BATCH_FILE_SETTINGS_RECONCILED,
    make_settings_autofill_config,
)

START_AT_VALUE = 25
START_AT_STRING = "25"

TEST_DEFAULT_UNIREF_VARIANT = "50"

UNIREF_FIXTURE_DIR = Path("tests") / "data" / "uniprot" / "uniref"


TEST_SETTINGS = frozendict(
    {**TEST_BATCH_FILE_SETTINGS, VARIANT: TEST_DEFAULT_UNIREF_VARIANT},
)

TEST_SETTINGS_RECONCILED = frozendict(
    {**TEST_BATCH_FILE_SETTINGS_RECONCILED, VARIANT: TEST_DEFAULT_UNIREF_VARIANT},
)

UNIREF_VARIANT_ALIASES = ["v", "variant"]


@pytest.fixture(params=UNIREF_VARIANTS)
def uniref_variant_value(request: pytest.FixtureRequest) -> str:
    """Parametrized fixture over all valid uniref variants."""
    return request.param


# Integration test for the UniRef pipeline
#
# Data output goes to a local DuckDB destination.


@pytest.fixture
def duckdb_uniref_settings(tmp_path: Path) -> UnirefSettings:
    """Provide UnirefSettings pointing at the real UniRef XML fixtures."""
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    return make_settings_autofill_config(  # type: ignore[reportReturnType]
        UnirefSettings,
        {
            VARIANT: TEST_DEFAULT_UNIREF_VARIANT,
            "input_dir": str(UNIREF_FIXTURE_DIR),
            "output_dir": str(output_dir),
            "use_destination": "local_fs",
            "use_output_dir_for_pipeline_metadata": False,
        },
    )


def _run_uniref_duckdb_pipeline(settings: UnirefSettings, tmp_path: Path) -> tuple[Any, Any]:
    """Run ``parse_uniref`` through a real DuckDB pipeline and return (pipeline, load_info)."""
    pipeline = dlt.pipeline(
        pipeline_name="test_uniref_pipeline",
        destination="duckdb",
        dataset_name="test_uniref",
        pipelines_dir=str(tmp_path / "pipelines"),
    )
    load_info = pipeline.run(parse_uniref(settings))
    return pipeline, load_info


def test_integration_uniref_pipeline_loads_into_duckdb(
    duckdb_uniref_settings: UnirefSettings,
    tmp_path: Path,
) -> None:
    """Full integration test: parse the real fixture files into a DuckDB destination.

    Confirms that the pipeline runs end-to-end without failed jobs and that the
    expected CDM tables are populated with data queried back out of DuckDB.
    """
    pipeline, load_info = _run_uniref_duckdb_pipeline(duckdb_uniref_settings, tmp_path)

    assert not load_info.has_failed_jobs

    # the parser emits an "entity" row per UniRef cluster; the fixtures contain
    # several clusters across the chunk files, so this table must not be empty.
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT COUNT(*) FROM entity") as cur,
    ):
        (entity_count,) = cur.fetchone()

    assert entity_count > 0

    # every entity must have a "uniref:"-prefixed entity_id and be a Cluster.
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT entity_id, entity_type FROM entity") as cur,
    ):
        rows = cur.fetchall()

    assert all(entity_id.startswith("uniref:") for entity_id, _ in rows)
    assert all(entity_type == "Cluster" for _, entity_type in rows)


def test_integration_uniref_pipeline_populates_related_tables(
    duckdb_uniref_settings: UnirefSettings,
    tmp_path: Path,
) -> None:
    """The pipeline should populate the related CDM tables produced by the parser."""
    pipeline, load_info = _run_uniref_duckdb_pipeline(duckdb_uniref_settings, tmp_path)

    assert not load_info.has_failed_jobs

    # tables that the parser always emits at least one row for, per entry
    for table_name in ("entity", "cluster", "clustermember", "entity_x_source_file"):
        with (
            pipeline.sql_client() as client,
            client.execute_query(f"SELECT COUNT(*) FROM {table_name}") as cur,  # noqa: S608
        ):
            (count,) = cur.fetchone()
        assert count > 0, f"expected rows in table '{table_name}'"

    # every cluster must carry the injected "UniRef <variant>" protocol label
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT DISTINCT protocol FROM cluster") as cur,
    ):
        protocols = {row[0] for row in cur.fetchall()}

    assert protocols == {f"UniRef {duckdb_uniref_settings.variant}"}

    # entity_x_source_file must reference the fixture files we loaded from
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT DISTINCT source_file FROM entity_x_source_file") as cur,
    ):
        source_files = {row[0] for row in cur.fetchall()}

    assert source_files
    assert all(sf.endswith(".xml") for sf in source_files)


def test_integration_cli_uniref_pipeline_output_validated(
    duckdb_uniref_settings: UnirefSettings,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Exercise the real ``cli()`` wiring end-to-end against a DuckDB destination."""
    monkeypatch.setattr(uniref_module, "UnirefSettings", MagicMock(return_value=duckdb_uniref_settings))

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        *,
        settings: UnirefSettings,
        resource: Any,
        pipeline_kwargs: dict[str, Any],
    ) -> None:
        """Replacement for core.run_pipeline that runs the resource through DuckDB."""
        assert settings is duckdb_uniref_settings
        assert pipeline_kwargs == {
            "pipeline_name": f"uniref_{duckdb_uniref_settings.variant}",
            "dataset_name": "uniprot_kb",
        }
        pipeline = dlt.pipeline(
            pipeline_name="test_uniref_cli_pipeline",
            destination="duckdb",
            dataset_name="test_uniref_cli",
            pipelines_dir=str(tmp_path / "pipelines"),
        )
        captured["load_info"] = pipeline.run(resource)
        captured["pipeline"] = pipeline

    monkeypatch.setattr(uniref_module, "run_pipeline", fake_run_pipeline)

    cli()

    # UnirefSettings was constructed with the dlt config coming from core
    uniref_module.UnirefSettings.assert_called_once_with()

    load_info = captured["load_info"]
    pipeline = captured["pipeline"]
    assert not load_info.has_failed_jobs

    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT COUNT(*) FROM entity") as cur,
    ):
        (entity_count,) = cur.fetchone()

    assert entity_count > 0
