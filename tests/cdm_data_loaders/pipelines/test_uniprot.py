"""Tests for the UniProt DLT pipeline."""

from collections.abc import Callable
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import dlt
import pytest
from frozendict import frozendict

from cdm_data_loaders.parsers.uniprot.uniprot_kb import ENTRY_XML_TAG
from cdm_data_loaders.pipelines import uniprot_kb as uniprot_module
from cdm_data_loaders.pipelines.uniprot_kb import (
    UNIPROT_LOG_INTERVAL,
    UniProtSettings,
    cli,
    parse_uniprot,
    run_uniprot_pipeline,
)
from tests.cdm_data_loaders.core.conftest import (
    TEST_BATCH_FILE_SETTINGS,
    TEST_BATCH_FILE_SETTINGS_RECONCILED,
    check_settings,
    make_settings_autofill_config,
)
from tests.conftest import TEST_DATA_DIR
from tests.helpers import assert_cli_field_roundtrips, assert_no_cli_clashes


@pytest.fixture
def test_settings() -> UniProtSettings:
    """Provide a minimal valid UniProtSettings object."""
    return make_settings_autofill_config(UniProtSettings)  # type: ignore[reportReturnType]


TEST_SETTINGS = frozendict(
    {**TEST_BATCH_FILE_SETTINGS, "log_interval": UNIPROT_LOG_INTERVAL},
)

TEST_SETTINGS_RECONCILED = frozendict({**TEST_BATCH_FILE_SETTINGS_RECONCILED, "log_interval": UNIPROT_LOG_INTERVAL})


def test_uniprot_settings_all_params_set() -> None:
    """Ensure that settings are set correctly when all args are specified.

    Note that TEST_SETTINGS includes a value for pipeline_dir.
    """
    s = make_settings_autofill_config(UniProtSettings, TEST_SETTINGS)  # type: ignore[reportReturnType]
    check_settings(s, TEST_SETTINGS_RECONCILED)


# model_fields
def test_no_cli_param_clashes() -> None:
    """Ensure that the settings object does not have any internal alias collisions."""
    assert_no_cli_clashes(UniProtSettings)


@pytest.mark.parametrize("field_name", list(TEST_SETTINGS))
def test_cli_fields_parse_correctly(field_name: str) -> None:
    """Ensure that all fields and variants are parsed correctly."""
    assert_cli_field_roundtrips(
        UniProtSettings, field_name, TEST_SETTINGS[field_name], TEST_SETTINGS_RECONCILED[field_name]
    )


def test_cli_passes_settings_class_to_run_cli() -> None:
    """Ensure that cli() calls run_cli with UniProtSettings and the UniProt log_interval default."""
    with patch.object(uniprot_module, "run_cli") as mock_run_cli:
        cli()

    mock_run_cli.assert_called_once()
    assert mock_run_cli.call_args[0] == (UniProtSettings, run_uniprot_pipeline)
    assert mock_run_cli.call_args.kwargs == {}


def test_cli_calls_run_uniprot_pipeline(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that cli() calls run_uniprot_pipeline with the test_settings."""
    mock_settings_instance = MagicMock()
    mock_settings_cls = MagicMock(return_value=mock_settings_instance)
    mock_run_uniprot_pipeline = MagicMock()

    monkeypatch.setattr(uniprot_module, "UniProtSettings", mock_settings_cls)
    monkeypatch.setattr(uniprot_module, "run_uniprot_pipeline", mock_run_uniprot_pipeline)

    cli()

    mock_settings_cls.assert_called_once_with()
    mock_run_uniprot_pipeline.assert_called_once_with(mock_settings_instance)


# Tests for running the pipeline itself
def test_run_uniprot_pipeline_args_set_correctly(test_settings: UniProtSettings) -> None:
    """Ensure that the pipeline arguments are set correctly, and each pipeline has a different name."""
    with patch.object(uniprot_module, "run_pipeline") as mock_run_pipeline:
        run_uniprot_pipeline(test_settings)

    assert mock_run_pipeline.call_count == 1
    _, kwargs = mock_run_pipeline.call_args
    assert kwargs.keys() == {"settings", "resource", "pipeline_kwargs"}
    assert kwargs["pipeline_kwargs"] == {
        "pipeline_name": "uniprot_kb",
        "dataset_name": "uniprot_kb",
    }
    assert kwargs["settings"] == test_settings
    assert isinstance(kwargs["resource"], Callable)


def test_run_uniprot_pipeline_sets_core_run_pipeline_args_correctly(
    test_settings: UniProtSettings, mock_dlt: MagicMock, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure that run_uniprot_pipeline calls core.run_pipeline with the correct args."""
    mock_parse_uniprot = MagicMock()
    monkeypatch.setattr(uniprot_module, "parse_uniprot", mock_parse_uniprot)

    run_uniprot_pipeline(test_settings)

    # parse_uniprot was called once with the test_settings to produce the resource
    mock_parse_uniprot.assert_called_once_with(test_settings)

    # the return value of parse_uniprot(test_settings) is what gets passed to pipeline.run
    expected_resource = mock_parse_uniprot.return_value

    mock_dlt.destination.assert_called_once_with(test_settings.use_destination)
    mock_dlt.pipeline.assert_called_once_with(
        destination=mock_dlt.destination.return_value,
        pipeline_name="uniprot_kb",
        dataset_name="uniprot_kb",
    )
    mock_dlt.pipeline.return_value.run.assert_called_once_with(expected_resource)


def test_parse_uniprot_resource(test_settings: UniProtSettings) -> None:
    """Ensure that parse_uniprot calls process_xml_file_batches with the namespaced UniProt XML tag."""
    with patch.object(uniprot_module, "process_xml_file_batches") as mock_stream:
        mock_stream.return_value = iter([])
        list(parse_uniprot(test_settings))

    assert mock_stream.call_count == 1
    kwargs = mock_stream.call_args.kwargs
    assert kwargs.keys() == {"settings", "xml_tag", "parse_fn"}
    assert kwargs["xml_tag"] == ENTRY_XML_TAG
    assert kwargs["settings"] == test_settings
    assert isinstance(kwargs["parse_fn"], Callable)


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
    """Provide UniProtSettings pointing at the real UniProt XML fixtures.

    ``input_dir`` points at the fixture directory containing ``chunk_0000N.xml``
    files, and ``output_dir`` is a local directory inside ``tmp_path`` so the run is
    fully isolated.
    """
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
    """Exercise the real ``cli()`` wiring end-to-end against a DuckDB destination.

    ``cli()`` builds a ``UniProtSettings`` from the dlt config / CLI args and then
    runs the pipeline via ``run_uniprot_pipeline`` -> ``core.run_pipeline``. We
    stub out ``UniProtSettings`` construction to return our fixture-backed
    settings, and redirect ``core.run_pipeline`` to a real DuckDB pipeline so the
    full flow (settings -> resource -> pipeline.run -> loaded data) is validated.
    """
    monkeypatch.setattr(uniprot_module, "UniProtSettings", MagicMock(return_value=duckdb_uniprot_settings))

    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        *,
        settings: UniProtSettings,
        resource: Any,  # noqa: ANN401
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
