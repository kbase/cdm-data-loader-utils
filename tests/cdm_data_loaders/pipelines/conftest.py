"""Shared fixtures for pipelines tests."""

import copy
import logging
import os
from collections.abc import Callable, Iterator
from itertools import count
from pathlib import Path
from typing import Any, Final
from unittest.mock import MagicMock, patch

import dlt
import pytest
from dlt.common.pipeline import LoadInfo

from cdm_data_loaders.core.fields import LOCAL_FS, S3
from cdm_data_loaders.pipelines import core
from cdm_data_loaders.pipelines.core import (
    DISABLE_COMPRESSION_ENV_VAR,
    LOAD_INFO_TABLE_NAME,
)
from cdm_data_loaders.utils import cdm_logger

DEFAULT_DLT_TABLES = {"_dlt_version", "_dlt_loads", "_dlt_pipeline_state", LOAD_INFO_TABLE_NAME}

HOOK_ENV_VAR: Final[str] = "RUNTIME__SLACK_INCOMING_HOOK"
HOOK_PART_ENV_VARS: Final[tuple[str, ...]] = ("VARIABLE_B", "VARIABLE_T", "CHAR_STR")
VALID_DESTINATIONS = [LOCAL_FS, S3]

TEST_LOG_CONFIG_FILE: Final[Path] = Path("tests") / "data" / "pipelines" / "log_config.json"


DLT_CONFIG: Final[dict[str, Any]] = {
    "destination": {
        LOCAL_FS: {"bucket_url": "/configured/out/"},
        "s3": {"bucket_url": "s3://bucket/prefix"},
    }
}


def duckdb_pipeline(tmp_path: Path, pipeline_name: str, dataset_name: str, db_name: str | None = None) -> dlt.Pipeline:
    """Return a dlt pipeline with a duckdb destination in tmp_path.

    The database file defaults to a name derived from the pipeline name so that the DuckDB
    catalog differs from the dataset schema; equal names trigger a DuckDB binder error.
    """
    db_file = f"{db_name or pipeline_name}_db.duckdb"
    return dlt.pipeline(
        pipeline_name=pipeline_name,
        dataset_name=dataset_name,
        destination=dlt.destinations.duckdb(str(tmp_path / db_file)),
        pipelines_dir=str(tmp_path / "pipelines"),
    )


@pytest.fixture(autouse=True)
def clean_env() -> Iterator[None]:
    """Remove env vars that core reads or writes, and restore the whole environment afterwards."""
    with patch.dict(os.environ):
        for name in (*HOOK_PART_ENV_VARS, HOOK_ENV_VAR, DISABLE_COMPRESSION_ENV_VAR):
            os.environ.pop(name, None)
        for name in [key for key in os.environ if key.startswith("CDL_")]:
            del os.environ[name]
        os.environ["RUNTIME__DLTHUB_TELEMETRY"] = "false"
        yield


@pytest.fixture
def slack_mock(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace the dlt Slack sender used by send_slack_message_carefully."""
    mock = MagicMock()
    monkeypatch.setattr(core, "send_slack_message", mock)
    return mock


@pytest.fixture
def careful_slack_mock(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace send_slack_message_carefully as called by run_pipeline."""
    mock = MagicMock()
    monkeypatch.setattr(core, "send_slack_message_carefully", mock)
    return mock


@pytest.fixture
def mock_init_logger(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Stop run_cli reconfiguring logging and point dlt.config at a plain dict."""
    mock = MagicMock()
    monkeypatch.setattr(core, "init_logger", mock)
    monkeypatch.setattr(dlt, "config", copy.deepcopy(DLT_CONFIG))
    return mock


@pytest.fixture
def mock_pipeline(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace dlt.pipeline with a factory whose pipeline has no Slack hook and returns no load info."""
    pipeline = MagicMock()
    pipeline.runtime_config.slack_incoming_hook = None
    pipeline.run.return_value = None
    factory = MagicMock(return_value=pipeline)
    monkeypatch.setattr(dlt, "pipeline", factory)
    return factory


@pytest.fixture
def mock_dlt(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace dlt.destination and dlt.pipeline with a single mock, as run_ncbi_pipeline uses both."""
    mock = MagicMock()
    mock.pipeline.return_value.runtime_config.slack_incoming_hook = None
    mock.pipeline.return_value.run.return_value = MagicMock()
    monkeypatch.setattr(dlt, "destination", mock.destination)
    monkeypatch.setattr(dlt, "pipeline", mock.pipeline)
    return mock


@pytest.fixture
def mock_build_destination(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Replace build_destination as called by run_pipeline."""
    mock = MagicMock()
    monkeypatch.setattr(core, "build_destination", mock)
    return mock


def isolated_run(
    tmp_path: Path,
    dataset_prefix: str,
    run_fn: Callable[..., LoadInfo | None],
    base_kwargs: dict[str, Any],
) -> Callable[..., tuple[LoadInfo | None, Path]]:
    """Return a factory that runs a pipeline with per-run isolated output state.

    Each call gets a unique output_dir, dataset_name, and log config file under
    tmp_path, so repeated runs never collide or accumulate rows.

    :param tmp_path: pytest tmp_path for this test
    :type tmp_path: Path
    :param dataset_prefix: prefix for each run's unique dataset name
    :type dataset_prefix: str
    :param run_fn: the pipeline run function to invoke with built settings
    :type run_fn: Callable[..., LoadInfo | None]
    :param base_kwargs: base settings kwargs; keys open to override per call
    :type base_kwargs: dict[str, Any]
    :return: factory(chunk_dir=None, **overrides) -> (load_info, output_dir)
    :rtype: Callable[..., tuple[LoadInfo | None, Path]]
    """
    run_counter = count()

    def _factory(**overrides: Any) -> tuple[LoadInfo | None, Path]:
        run_index = next(run_counter)
        output_dir = tmp_path / f"output_{run_index}"
        output_dir.mkdir()
        log_config_file = tmp_path / f"logging_{run_index}.conf"
        log_config_file.touch()
        kwargs: dict[str, Any] = {
            "dataset_name": f"{dataset_prefix}_{run_index}",
            "log_config_file": str(log_config_file),
            "output_dir": str(output_dir),
            **base_kwargs,
        }
        kwargs.update(overrides)
        settings_cls = kwargs.pop("settings_cls")
        load_info = run_fn(settings_cls(**kwargs))
        return load_info, output_dir

    return _factory


@pytest.fixture
def isolated_run_factory(tmp_path: Path) -> Callable[..., Callable[..., tuple[LoadInfo | None, Path]]]:
    """Return a builder for isolated pipeline run factories bound to this test's tmp_path."""

    def _build(
        dataset_prefix: str, run_fn: Callable[..., LoadInfo | None], base_kwargs: dict[str, Any]
    ) -> Callable[..., tuple[LoadInfo | None, Path]]:
        return isolated_run(tmp_path, dataset_prefix, run_fn, base_kwargs)

    return _build


@pytest.fixture
def restore_logging() -> Iterator[None]:
    """Undo the global logging changes made by init_logger, so caplog and later tests are unaffected."""
    root = logging.getLogger()
    handlers, level = root.handlers[:], root.level
    applied_key = cdm_logger.applied_key
    yield
    for handler in root.handlers[:]:
        root.removeHandler(handler)
    for handler in handlers:
        root.addHandler(handler)
    root.setLevel(level)
    cdm_logger.applied_key = applied_key
