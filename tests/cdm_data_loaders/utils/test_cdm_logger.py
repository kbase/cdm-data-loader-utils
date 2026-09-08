"""Tests for the CDM logger."""

import json
import logging
from collections.abc import Generator
from contextlib import nullcontext
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
import yaml
from frozendict import frozendict

import cdm_data_loaders.utils.cdm_logger as cdm_logger_module
from cdm_data_loaders.utils.cdm_logger import (
    LoggerSettings,
    _load_config_from_path,
    configure_root_logger_from_dlt,
    get_cdm_logger,
)


def get_logger_dict() -> dict[str, Any]:
    """Get the dictionary of existing loggers."""
    return logging.Logger.manager.loggerDict


@pytest.fixture(autouse=True)
def reset_module_state(monkeypatch: pytest.MonkeyPatch) -> Generator[None, Any]:
    """Reset global module state and root logging config between every test."""
    # Reset the guard flag so each test starts from a clean slate
    monkeypatch.setattr(cdm_logger_module, "ROOT_LOGGER_CONFIGURED", False)

    # Clear root logger handlers so dictConfig calls don't accumulate
    root_logger: logging.Logger = logging.getLogger("root")
    root_logger.handlers.clear()

    yield

    # Teardown: reset again after test
    root_logger.handlers.clear()


@pytest.fixture
def mock_config_root_logger() -> Generator[MagicMock, Any]:
    """Patch dlt's config_root_logger so tests are hermetic."""
    with patch("cdm_data_loaders.utils.cdm_logger.config_root_logger") as mock:
        yield mock


@pytest.fixture
def mock_dlt_logger(monkeypatch: pytest.MonkeyPatch) -> MagicMock:
    """Create a mock DLT logger with everyone's fave mocker, MagicMock."""
    dlt_mock = MagicMock()
    dlt_mock.propagate = False
    dlt_mock.handlers = [1, 2, 3, 4]

    monkeypatch.setitem(logging.Logger.manager.loggerDict, "dlt", dlt_mock)
    return dlt_mock


@pytest.fixture
def clean_env(monkeypatch: pytest.MonkeyPatch) -> None:
    """Remove logging-related env vars during tests."""
    monkeypatch.delenv("LOG_CONFIG_FILE", raising=False)
    monkeypatch.delenv("LOG-CONFIG-FILE", raising=False)


@pytest.fixture
def empty_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Remove LOG_CONFIG_FILE env vars and chdir to an empty temporary directory with no config file."""
    monkeypatch.delenv("LOG_CONFIG_FILE", raising=False)
    monkeypatch.delenv("LOG-CONFIG-FILE", raising=False)
    monkeypatch.chdir(tmp_path)


YAML_CONFIG = frozendict(
    {
        "version": 1,
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "stream": "ext://sys.stdout",
            }
        },
        "loggers": {"root": {"level": "DEBUG", "handlers": ["console"]}},
    }
)


JSON_CONFIG = frozendict(
    {
        "version": 1,
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "stream": "ext://sys.stdout",
            }
        },
        "loggers": {"root": {"level": "WARNING", "handlers": ["console"]}},
    }
)


@pytest.fixture
def yaml_config_file(tmp_path: Path) -> Path:
    """Write a minimal valid YAML logging config and return its path."""
    p = tmp_path / "logging.yaml"
    p.write_text(yaml.dump(dict(YAML_CONFIG)))
    return p


@pytest.fixture
def json_config_file(tmp_path: Path) -> Path:
    """Write a minimal valid JSON logging config and return its path."""
    p = tmp_path / "logging.json"
    p.write_text(json.dumps(dict(JSON_CONFIG)))
    return p


@pytest.fixture
def settings_with_yaml(yaml_config_file: Path) -> LoggerSettings:
    """Logger settings object for the YAML config file."""
    return LoggerSettings(log_config_file=str(yaml_config_file))


@pytest.fixture
def settings_with_json(json_config_file: Path) -> LoggerSettings:
    """Logger settings object for the json config file."""
    return LoggerSettings(log_config_file=str(json_config_file))


@pytest.fixture
def settings_without_config() -> LoggerSettings:
    """Logger settings object with no log config file."""
    return LoggerSettings(log_config_file=None)


FILE_CONTENT = frozendict(
    {
        "valid.json": json.dumps(dict(JSON_CONFIG)),
        "valid.yml": yaml.dump(dict(YAML_CONFIG)),
        "valid.yaml": yaml.dump({"version": 1, "loggers": {}}),
        "invalid.json": "{not valid json",
        "invalid.yaml": ": : invalid: yaml: : :",
    }
)

# LoggerSettings


def test_logger_settings_with_log_config_file(yaml_config_file: Path) -> None:
    """LoggerSettings stores the log_config_file path passed directly to the constructor."""
    settings = LoggerSettings(log_config_file=str(yaml_config_file))
    assert settings.log_config_file == str(yaml_config_file)


@pytest.mark.parametrize("log_config_file_is_none", [True, False])
def test_logger_settings_with_without_params(log_config_file_is_none: bool) -> None:
    """LoggerSettings.log_config_file is None when no value is provided."""
    settings = LoggerSettings(log_config_file=None) if log_config_file_is_none else LoggerSettings()  # pyright: ignore[reportCallIssue]
    assert settings.log_config_file is None


# _load_config_from_path


@pytest.mark.parametrize(
    ("file_name", "expected"),
    [
        ("logger_config", pytest.raises(ValueError, match=r"Unsupported config file format: logger_config")),
        ("logger.toml", pytest.raises(ValueError, match=r"Unsupported config file format: logger.toml")),
        ("path/to/logger.cfg", pytest.raises(ValueError, match=r"Unsupported config file format: logger.cfg")),
        ("path/to/file.json", pytest.raises(FileNotFoundError, match="No such file or directory")),
        ("valid.json", nullcontext(dict(JSON_CONFIG))),
        ("valid.yaml", nullcontext({"version": 1, "loggers": {}})),
        ("valid.yml", nullcontext(dict(YAML_CONFIG))),
        (
            "invalid.json",
            pytest.raises(json.decoder.JSONDecodeError, match="Expecting property name enclosed in double quotes"),
        ),
        ("invalid.yaml", pytest.raises(yaml.YAMLError, match=r"expected <block end>, but found \':\'\n")),
    ],
)
def test_load_config_from_path_pass_fail(
    tmp_path: Path, file_name: str, expected: nullcontext | pytest.RaisesExc
) -> None:
    """Test the loading of a config file, both successfully and with failure states."""
    p = tmp_path / file_name
    if file_name in FILE_CONTENT:
        p.write_text(FILE_CONTENT[file_name])

    with expected as e:
        assert _load_config_from_path(p) == e


# configure_root_logger_from_dlt
def test_configure_root_logger_from_dlt_configures_dlt_on_first_call(
    mock_config_root_logger: MagicMock, mock_dlt_logger: MagicMock
) -> None:
    """The first call sets dlt logger propagation, invokes config_root_logger, and sets the guard."""
    # check what is initialised in the way of loggers
    assert "dlt" in logging.Logger.manager.loggerDict

    configure_root_logger_from_dlt()

    mock_config_root_logger.assert_called_once()
    dlt_logger = logging.getLogger("dlt")
    assert dlt_logger == mock_dlt_logger
    assert dlt_logger.propagate is True
    assert dlt_logger.handlers == []
    assert cdm_logger_module.ROOT_LOGGER_CONFIGURED is True


def test_configure_root_logger_from_dlt_does_not_reconfigure_on_subsequent_calls(
    mock_config_root_logger: MagicMock,
) -> None:
    """Calling configure_root_logger_from_dlt a second time is a no-op; config_root_logger is only called once."""
    configure_root_logger_from_dlt()
    configure_root_logger_from_dlt()

    mock_config_root_logger.assert_called_once()


def test_configure_root_logger_skips_when_guard_already_set(
    mock_config_root_logger: MagicMock, monkeypatch: pytest.MonkeyPatch
) -> None:
    """When ROOT_LOGGER_CONFIGURED is already True, config_root_logger is never called."""
    monkeypatch.setattr(cdm_logger_module, "ROOT_LOGGER_CONFIGURED", True)

    configure_root_logger_from_dlt()

    mock_config_root_logger.assert_not_called()


# get_cdm_logger
def test_get_cdm_logger_calls_init_logger_with_args() -> None:
    """get_cdm_logger should call init_logger to ensure logging is initialised."""
    settings = LoggerSettings()  # pyright: ignore[reportCallIssue]

    with patch("cdm_data_loaders.utils.cdm_logger.init_logger") as mock_init:
        logger = get_cdm_logger("test.logger", settings=settings)

    mock_init.assert_called_once_with(settings)
    assert isinstance(logger, logging.Logger)
    assert logger.name == "test.logger"


def test_get_cdm_logger_without_settings_passes_none_to_init_logger() -> None:
    """get_cdm_logger called with no settings argument should pass None to init_logger."""
    with patch("cdm_data_loaders.utils.cdm_logger.init_logger") as mock_init:
        logger = get_cdm_logger("another.logger")

    mock_init.assert_called_once_with(None)
    assert isinstance(logger, logging.Logger)
    assert logger.name == "another.logger"


@pytest.mark.usefixtures("mock_config_root_logger")
def test_get_cdm_logger_returns_logger_with_correct_name(settings_with_yaml: LoggerSettings) -> None:
    """The returned object is a Logger whose name matches the argument passed in."""
    result = get_cdm_logger("my.service", settings=settings_with_yaml)

    assert isinstance(result, logging.Logger)
    assert result.name == "my.service"


@pytest.mark.usefixtures("mock_config_root_logger")
def test_get_cdm_logger_uses_yaml_config_when_provided(
    settings_with_yaml: LoggerSettings,
) -> None:
    """When settings point to a YAML file, dictConfig receives the parsed YAML content."""
    with patch("cdm_data_loaders.utils.cdm_logger.logging.config.dictConfig") as mock_dict_config:
        get_cdm_logger("svc", settings=settings_with_yaml)

    applied_config = mock_dict_config.call_args[0][0]
    assert applied_config["disable_existing_loggers"] is False
    assert applied_config["loggers"]["root"]["level"] == "DEBUG"


@pytest.mark.usefixtures("mock_config_root_logger")
def test_get_cdm_logger_uses_json_config_when_provided(
    settings_with_json: LoggerSettings,
) -> None:
    """When settings point to a JSON file, dictConfig receives the parsed JSON content."""
    with patch("cdm_data_loaders.utils.cdm_logger.logging.config.dictConfig") as mock_dict_config:
        get_cdm_logger("svc", settings=settings_with_json)

    applied_config = mock_dict_config.call_args[0][0]
    assert applied_config["loggers"]["root"]["level"] == "WARNING"


@pytest.mark.usefixtures("mock_config_root_logger")
def test_get_cdm_logger_falls_back_to_builtin_config_when_no_settings(caplog: pytest.LogCaptureFixture) -> None:
    """Passing settings=None applies the built-in INFO config and emits a warning."""
    with (
        patch("cdm_data_loaders.utils.cdm_logger.logging.config.dictConfig") as mock_dict_config,
        caplog.at_level(logging.WARNING, logger="cdm_data_loaders.utils.cdm_logger"),
    ):
        get_cdm_logger("svc", settings=None)

    applied_config = mock_dict_config.call_args[0][0]
    assert applied_config["loggers"]["root"]["level"] == "INFO"
    assert any("Falling back to built-in config" in m for m in caplog.messages)


@pytest.mark.usefixtures("mock_config_root_logger")
def test_get_cdm_logger_falls_back_to_builtin_config_when_log_config_file_is_none(
    settings_without_config: LoggerSettings, caplog: pytest.LogCaptureFixture
) -> None:
    """Settings with log_config_file=None applies the built-in INFO config and emits a warning."""
    with (
        patch("cdm_data_loaders.utils.cdm_logger.logging.config.dictConfig") as mock_dict_config,
        caplog.at_level(logging.WARNING, logger="cdm_data_loaders.utils.cdm_logger"),
    ):
        get_cdm_logger("svc", settings=settings_without_config)

    applied_config = mock_dict_config.call_args[0][0]
    assert applied_config["loggers"]["root"]["level"] == "INFO"
    assert any("Falling back to built-in config" in m for m in caplog.messages)


@pytest.mark.usefixtures("mock_config_root_logger")
def test_get_cdm_logger_disable_existing_loggers_always_set_false(
    settings_with_yaml: LoggerSettings,
) -> None:
    """disable_existing_loggers is always injected as False to prevent silently killing pre-existing loggers."""
    with patch("cdm_data_loaders.utils.cdm_logger.logging.config.dictConfig") as mock_dict_config:
        get_cdm_logger("svc", settings=settings_with_yaml)

    applied_config = mock_dict_config.call_args[0][0]
    assert applied_config["disable_existing_loggers"] is False


@pytest.mark.usefixtures("mock_config_root_logger")
def test_get_cdm_logger_calls_configure_root_logger_from_dlt(
    settings_with_yaml: LoggerSettings,
) -> None:
    """get_cdm_logger always delegates to configure_root_logger_from_dlt after applying config."""
    with patch("cdm_data_loaders.utils.cdm_logger.configure_root_logger_from_dlt") as mock_configure:
        get_cdm_logger("svc", settings=settings_with_yaml)

    mock_configure.assert_called_once()


def test_get_cdm_logger_dlt_already_initialised_does_not_reconfigure(
    settings_with_yaml: LoggerSettings, mock_config_root_logger: MagicMock
) -> None:
    """When the dlt logger is already initialised, repeated calls to get_cdm_logger do not re-invoke config_root_logger."""
    cdm_logger_module.ROOT_LOGGER_CONFIGURED = True

    get_cdm_logger("svc", settings=settings_with_yaml)
    get_cdm_logger("svc", settings=settings_with_yaml)

    mock_config_root_logger.assert_not_called()


def test_get_cdm_logger_dlt_not_yet_initialised_configures_on_first_call(
    settings_with_yaml: LoggerSettings, mock_config_root_logger: MagicMock
) -> None:
    """When the dlt logger has not yet been initialised, the first call to get_cdm_logger configures it and sets the guard."""
    assert cdm_logger_module.ROOT_LOGGER_CONFIGURED is False

    get_cdm_logger("svc", settings=settings_with_yaml)

    mock_config_root_logger.assert_called_once()
    assert cdm_logger_module.ROOT_LOGGER_CONFIGURED is True
