"""Shared fixtures for pipelines tests."""

from itertools import batched
from pathlib import Path
from typing import Any, Final
from unittest.mock import MagicMock

import pytest

from cdm_data_loaders.pipelines import core

START_AT_VALUE: Final[int] = 50
START_AT_STRING: Final[str] = "50"
TEST_LOG_CONFIG_FILE: Final[str] = "log_conf.json"


def make_batcher(files: list[Path], batch_size: int = 5) -> MagicMock:
    """Return a mock NumericFileSequenceBatcher that yields ``files`` in batches then an empty list."""
    batches = [list(b) for b in batched(files, batch_size, strict=False)]
    mock_batcher = MagicMock()
    mock_batcher.get_batch.side_effect = [*batches, []]
    return mock_batcher


@pytest.fixture
def fake_files() -> list[Path]:
    """List of five files, used for testing."""
    return [Path(f"/fake/input/part_{n}.xml") for n in [1, 2, 3, 4, 5]]


@pytest.fixture(autouse=True)
def mock_init_logger(monkeypatch: pytest.MonkeyPatch) -> None:
    """Mock the init_logger call in core to prevent the logger from trying to initialise itself every time."""
    monkeypatch.setattr(core, "init_logger", MagicMock())


@pytest.fixture(autouse=True)
def mock_send_slack_message(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Patch send_slack_message in core to prevent undue slack notifications."""
    slack_mock = MagicMock()
    monkeypatch.setattr(core, "send_slack_message", slack_mock)
    return slack_mock


@pytest.fixture(autouse=True)
def patch_dlt_config(dlt_config: dict[str, Any], monkeypatch: pytest.MonkeyPatch) -> None:
    """Monkeypatch the dlt config in all tests."""
    monkeypatch.setattr(core.dlt, "config", dlt_config)


@pytest.fixture
def mock_dlt(monkeypatch: pytest.MonkeyPatch, dlt_config: dict[str, Any]) -> MagicMock:
    """Patch dlt in core, wiring pipeline.return_value to a fresh MagicMock."""
    dlt_mock = MagicMock()
    # patch the slack_incoming_hook config value so that tests do not send slack notifications
    dlt_mock.pipeline.return_value.runtime_config.slack_incoming_hook = None
    # patch the config in case dlt is used when initialising a settings object
    dlt_mock.config = dlt_config
    monkeypatch.setattr(core, "dlt", dlt_mock)
    return dlt_mock


@pytest.fixture
def patched_io(monkeypatch: pytest.MonkeyPatch) -> tuple[MagicMock, MagicMock]:
    """Patch NumericFileSequenceBatcher and stream_xml_file inside core; return (mock_batcher_cls, mock_stream)."""
    mock_batcher_cls = MagicMock()
    mock_stream = MagicMock(return_value=[])
    monkeypatch.setattr(core, "NumericFileSequenceBatcher", mock_batcher_cls)
    monkeypatch.setattr(core, "stream_xml_file", mock_stream)
    return mock_batcher_cls, mock_stream


@pytest.fixture
def patched_io_empty_batcher(monkeypatch: pytest.MonkeyPatch) -> tuple[MagicMock, MagicMock]:
    """Like patched_io but NumericFileSequenceBatcher immediately returns an empty batch."""
    mock_batcher_cls = MagicMock()
    mock_stream = MagicMock(return_value=[])
    mock_batcher_cls.return_value = make_batcher([])
    monkeypatch.setattr(core, "NumericFileSequenceBatcher", mock_batcher_cls)
    monkeypatch.setattr(core, "stream_xml_file", mock_stream)
    return mock_batcher_cls, mock_stream
