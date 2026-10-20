"""Tests for cdm_data_loaders.utils.cdm_logger."""

import copy
import json
import logging
import os
import re
import sys
from pathlib import Path
from types import TracebackType
from typing import Any

import pytest
import yaml
from dlt.common.runtime.json_logging import JSONLogFormatter
from frozendict import frozendict

from cdm_data_loaders.core.settings import LoggerSettings
from cdm_data_loaders.utils.cdm_logger import (
    JSON_FORMATTER_NAME,
    _default_config,
    _enforce_json_formatter,
    _load_config_from_path,
    adopt_dlt_logger,
    init_logger,
    log_exception,
)
from tests.helpers import jsonlines_text_to_objects, write_json_or_yaml

pytestmark = pytest.mark.usefixtures("isolated_logging")

NAMED_LOGGER = "cdm_test.named"
TEN_SECONDS_NS = 10_000_000_000

FALLBACK_RECORD = frozendict(
    {
        "type": "log",
        "logger": "cdm_data_loaders.utils.cdm_logger",
        "level": "WARNING",
        "msg": "No logging config file found. Falling back to built-in config.",
    }
)

ExcInfo = tuple[type[BaseException], BaseException, TracebackType | None]


@pytest.fixture
def default_excepthook_calls(monkeypatch: pytest.MonkeyPatch) -> list[ExcInfo]:
    """Replace sys.__excepthook__ with a recorder and return the calls it receives."""
    calls: list[ExcInfo] = []
    monkeypatch.setattr(sys, "__excepthook__", lambda *args: calls.append(args))
    return calls


def make_config(level: str = "DEBUG") -> dict[str, Any]:
    """Build a logging config with a plain-text root handler and a named logger, both writing to stdout."""
    return {
        "version": 1,
        "formatters": {"plain": {"format": "%(levelname)s %(message)s"}},
        "handlers": {
            "out": {"class": "logging.StreamHandler", "stream": "ext://sys.stdout", "formatter": "plain"},
            "named_out": {"class": "logging.StreamHandler", "stream": "ext://sys.stdout"},
        },
        "loggers": {NAMED_LOGGER: {"level": "WARNING", "handlers": ["named_out"], "propagate": False}},
        "root": {"level": level, "handlers": ["out"]},
    }


def settings_for(path: Path | str | None) -> LoggerSettings:
    """Build LoggerSettings pointing at ``path``, or with no config file if ``path`` is None."""
    return LoggerSettings(log_config_file=None if path is None else str(path))


def set_mtime(path: Path, mtime_ns: int) -> None:
    """Set the access and modification times of ``path``."""
    os.utime(path, ns=(mtime_ns, mtime_ns))


def summarise(record: dict[str, Any]) -> dict[str, Any]:
    """Keep only the fields of a JSON log record that are the same on every run."""
    return {key: record[key] for key in ("type", "logger", "level", "msg")}


def raised(exc: BaseException) -> ExcInfo:
    """Raise and catch ``exc`` so that it has a real traceback."""
    try:
        raise exc  # noqa: TRY301
    except BaseException as caught:  # noqa: BLE001
        return type(caught), caught, caught.__traceback__


def test_default_config_pass_returns_fresh_copy() -> None:
    """The built-in config has the expected content, and changing one copy does not affect the next."""
    expected = {
        "version": 1,
        "handlers": {"console": {"class": "logging.StreamHandler", "stream": "ext://sys.stdout"}},
        "root": {"level": "WARNING", "handlers": ["console"]},
        "disable_existing_loggers": False,
    }
    config = _default_config()
    assert config == expected

    config["handlers"]["console"]["stream"] = "ext://sys.stderr"
    config["root"]["handlers"].append("other")
    assert _default_config() == expected


@pytest.mark.parametrize("as_str", [pytest.param(False, id="path_object"), pytest.param(True, id="str_path")])
@pytest.mark.parametrize(
    "file_name",
    [
        pytest.param("logging.json", id="json"),
        pytest.param("logging.yaml", id="yaml"),
        pytest.param("logging.yml", id="yml"),
        pytest.param("LOGGING.YAML", id="uppercase_suffix"),
    ],
)
def test_load_config_from_path_pass_valid_file(tmp_path: Path, file_name: str, as_str: bool) -> None:
    """Valid JSON and YAML files are parsed into the config they contain."""
    path = write_json_or_yaml(tmp_path / file_name, make_config())
    assert _load_config_from_path(str(path) if as_str else path) == make_config()


@pytest.mark.parametrize(
    "file_name",
    [
        pytest.param("logger_config", id="no_suffix"),
        pytest.param("logger.toml", id="toml"),
        pytest.param("nested/logger.cfg", id="cfg_in_subdir"),
        pytest.param("logging.yaml.bak", id="backup_suffix"),
        pytest.param(".yaml", id="dotfile_named_yaml"),
    ],
)
def test_load_config_from_path_fail_unsupported_format(tmp_path: Path, file_name: str) -> None:
    """Unsupported suffixes are rejected before the file is opened, so the file doesn't need to exist."""
    expected = f"Unsupported config file format: {Path(file_name).name}"
    with pytest.raises(ValueError, match=re.escape(expected)):
        _load_config_from_path(tmp_path / file_name)


@pytest.mark.parametrize(
    "file_name", [pytest.param("missing.json", id="json"), pytest.param("missing.yaml", id="yaml")]
)
def test_load_config_from_path_fail_missing_file(tmp_path: Path, file_name: str) -> None:
    """A supported suffix with no file behind it raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        _load_config_from_path(tmp_path / file_name)


@pytest.mark.parametrize(
    ("file_name", "content", "expected"),
    [
        pytest.param("invalid.json", "{not valid json", json.JSONDecodeError, id="json"),
        pytest.param("invalid.yaml", ": : invalid: yaml: : :", yaml.YAMLError, id="yaml"),
    ],
)
def test_load_config_from_path_fail_invalid_syntax(
    tmp_path: Path, file_name: str, content: str, expected: type[Exception]
) -> None:
    """Files that can't be parsed raise the parser's error."""
    path = tmp_path / file_name
    path.write_text(content, encoding="utf-8")
    with pytest.raises(expected):
        _load_config_from_path(path)


@pytest.mark.parametrize(
    ("file_name", "content", "type_name"),
    [
        pytest.param("empty.yaml", "", "NoneType", id="empty_yaml"),
        pytest.param("list.json", "[]", "list", id="json_list"),
        pytest.param("scalar.yml", "just a string", "str", id="yaml_scalar"),
    ],
)
def test_load_config_from_path_fail_not_a_mapping(tmp_path: Path, file_name: str, content: str, type_name: str) -> None:
    """Files that parse to something other than a mapping raise a ValueError naming the file and type."""
    path = tmp_path / file_name
    path.write_text(content, encoding="utf-8")
    expected = f"Logging config in {file_name} must be a mapping, not {type_name}"
    with pytest.raises(TypeError, match=re.escape(expected)):
        _load_config_from_path(path)


@pytest.mark.parametrize(
    ("config", "expected"),
    [
        pytest.param(
            {"version": 1},
            {
                "version": 1,
                "formatters": {JSON_FORMATTER_NAME: {"()": JSONLogFormatter}},
                "disable_existing_loggers": False,
            },
            id="minimal_config",
        ),
        pytest.param(
            {"version": 1, "handlers": {}, "disable_existing_loggers": False},
            {
                "version": 1,
                "handlers": {},
                "formatters": {JSON_FORMATTER_NAME: {"()": JSONLogFormatter}},
                "disable_existing_loggers": False,
            },
            id="no_handlers",
        ),
        pytest.param(
            {
                "version": 1,
                "formatters": {"plain": {"format": "%(message)s"}},
                "handlers": {
                    "out": {"class": "logging.StreamHandler", "formatter": "plain"},
                    "null": {"class": "logging.NullHandler"},
                },
                "root": {"level": "INFO", "handlers": ["out", "null"]},
                "disable_existing_loggers": True,
            },
            {
                "version": 1,
                "formatters": {"plain": {"format": "%(message)s"}, JSON_FORMATTER_NAME: {"()": JSONLogFormatter}},
                "handlers": {
                    "out": {"class": "logging.StreamHandler", "formatter": JSON_FORMATTER_NAME},
                    "null": {"class": "logging.NullHandler", "formatter": JSON_FORMATTER_NAME},
                },
                "root": {"level": "INFO", "handlers": ["out", "null"]},
                "disable_existing_loggers": False,
            },
            id="handlers_with_and_without_formatter",
        ),
        pytest.param(
            {"version": 1, "formatters": {JSON_FORMATTER_NAME: {"format": "%(message)s"}}},
            {
                "version": 1,
                "formatters": {JSON_FORMATTER_NAME: {"()": JSONLogFormatter}},
                "disable_existing_loggers": False,
            },
            id="replaces_formatter_with_same_name",
        ),
    ],
)
def test_enforce_json_formatter_pass(config: dict[str, Any], expected: dict[str, Any]) -> None:
    """Every handler gets the JSON formatter, disable_existing_loggers is forced off, and the input is unchanged."""
    original = copy.deepcopy(config)
    assert _enforce_json_formatter(config) == expected
    assert config == original


@pytest.mark.parametrize(
    ("n_handlers", "propagate"),
    [
        pytest.param(2, False, id="own_handlers_no_propagation"),
        pytest.param(0, True, id="already_adopted"),
    ],
)
def test_adopt_dlt_logger_pass(n_handlers: int, propagate: bool) -> None:
    """The dlt logger loses its handlers and propagates to root, and calling it again changes nothing."""
    dlt_logger = logging.getLogger("dlt")
    dlt_logger.handlers.clear()
    for _ in range(n_handlers):
        dlt_logger.addHandler(logging.NullHandler())
    dlt_logger.propagate = propagate

    for _ in range(2):
        adopt_dlt_logger()
        assert dlt_logger.handlers == []
        assert dlt_logger.propagate is True


@pytest.mark.parametrize(
    "pass_settings", [pytest.param(False, id="settings_none"), pytest.param(True, id="log_config_file_none")]
)
def test_init_logger_pass_fallback_config(capsys: pytest.CaptureFixture[str], pass_settings: bool) -> None:
    """With no config file, root gets one JSON-formatted stdout handler at INFO and a JSON warning is logged."""
    init_logger(settings_for(None) if pass_settings else None)

    root = logging.getLogger()
    assert root.level == logging.WARNING
    assert len(root.handlers) == 1
    handler = root.handlers[0]
    assert isinstance(handler, logging.StreamHandler)
    assert isinstance(handler.formatter, JSONLogFormatter)
    assert handler.stream is sys.stdout
    assert [summarise(r) for r in jsonlines_text_to_objects(capsys.readouterr().out)] == [FALLBACK_RECORD]


@pytest.mark.parametrize(
    "file_name",
    [
        pytest.param("logging.yaml", id="yaml"),
        pytest.param("logging.yml", id="yml"),
        pytest.param("logging.json", id="json"),
    ],
)
def test_init_logger_pass_config_file(tmp_path: Path, capsys: pytest.CaptureFixture[str], file_name: str) -> None:
    """A config file is applied, its handlers are switched to JSON, and no fallback warning is logged."""
    init_logger(settings_for(write_json_or_yaml(tmp_path / file_name, make_config("DEBUG"))))

    root = logging.getLogger()
    named = logging.getLogger(NAMED_LOGGER)
    assert root.level == logging.DEBUG
    assert (named.level, named.propagate) == (logging.WARNING, False)
    assert [type(h.formatter) for h in [*root.handlers, *named.handlers]] == [JSONLogFormatter, JSONLogFormatter]

    logging.getLogger("cdm_test.child").debug("debug message")
    assert [summarise(r) for r in jsonlines_text_to_objects(capsys.readouterr().out)] == [
        {"type": "log", "logger": "cdm_test.child", "level": "DEBUG", "msg": "debug message"}
    ]


def test_init_logger_pass_reads_config_file_from_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """With no settings passed, the config file named in CDL_LOG_CONFIG_FILE is applied."""
    path = write_json_or_yaml(tmp_path / "logging.yaml", make_config("DEBUG"))
    monkeypatch.setenv("CDL_LOG_CONFIG_FILE", str(path))

    init_logger()

    assert logging.getLogger().level == logging.DEBUG


@pytest.mark.parametrize("use_file", [pytest.param(False, id="fallback"), pytest.param(True, id="config_file")])
def test_init_logger_pass_skips_repeat_call(tmp_path: Path, capsys: pytest.CaptureFixture[str], use_file: bool) -> None:
    """A repeat call with the same config keeps the existing handlers and logs nothing."""
    path = write_json_or_yaml(tmp_path / "logging.yaml", make_config()) if use_file else None
    init_logger(settings_for(path))
    root = logging.getLogger()
    handlers = list(root.handlers)
    capsys.readouterr()

    init_logger(settings_for(path))

    assert root.handlers == handlers
    assert capsys.readouterr().out == ""


def test_init_logger_pass_same_file_different_path_spelling(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Relative and absolute paths to the same file count as the same config."""
    monkeypatch.chdir(tmp_path)
    path = write_json_or_yaml(tmp_path / "logging.yaml", make_config())
    init_logger(settings_for("logging.yaml"))
    root = logging.getLogger()
    handlers = list(root.handlers)

    init_logger(settings_for(path))

    assert root.handlers == handlers


def test_init_logger_pass_reapplies_when_mtime_changes(tmp_path: Path) -> None:
    """Editing the config file (which changes its mtime) causes the next call to reapply it."""
    path = write_json_or_yaml(tmp_path / "logging.yaml", make_config("INFO"))
    original_mtime = path.stat().st_mtime_ns
    init_logger(settings_for(path))
    root = logging.getLogger()
    handlers = list(root.handlers)

    write_json_or_yaml(path, make_config("DEBUG"))
    set_mtime(path, original_mtime + TEN_SECONDS_NS)
    init_logger(settings_for(path))

    assert root.level == logging.DEBUG
    assert len(root.handlers) == 1
    assert root.handlers[0] not in handlers


def test_init_logger_pass_same_mtime_skipped_unless_forced(tmp_path: Path) -> None:
    """A content change that keeps the mtime is skipped, but force=True applies it."""
    path = write_json_or_yaml(tmp_path / "logging.yaml", make_config("INFO"))
    original_mtime = path.stat().st_mtime_ns
    init_logger(settings_for(path))
    root = logging.getLogger()
    handlers = list(root.handlers)

    write_json_or_yaml(path, make_config("DEBUG"))
    set_mtime(path, original_mtime)
    init_logger(settings_for(path))
    assert root.level == logging.INFO
    assert root.handlers == handlers

    init_logger(settings_for(path), force=True)
    assert root.level == logging.DEBUG
    assert len(root.handlers) == 1
    assert root.handlers[0] not in handlers


def test_init_logger_pass_force_relogs_fallback_warning(capsys: pytest.CaptureFixture[str]) -> None:
    """Forcing the fallback config rebuilds the handler and logs the warning again."""
    init_logger(settings_for(None))
    root = logging.getLogger()
    handlers = list(root.handlers)

    init_logger(settings_for(None), force=True)

    assert len(root.handlers) == 1
    assert root.handlers[0] not in handlers
    assert [summarise(r) for r in jsonlines_text_to_objects(capsys.readouterr().out)] == [
        FALLBACK_RECORD,
        FALLBACK_RECORD,
    ]


def test_init_logger_pass_reapplies_after_switching_config(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    """Switching fallback -> file -> fallback reapplies the config each time, with only the latest remembered."""
    path = write_json_or_yaml(tmp_path / "logging.yaml", make_config("DEBUG"))
    root = logging.getLogger()
    levels: list[int] = []
    seen: list[logging.Handler] = []

    for settings in (settings_for(None), settings_for(path), settings_for(None)):
        init_logger(settings)
        assert all(h not in seen for h in root.handlers)
        seen.extend(root.handlers)
        levels.append(root.level)

    assert levels == [logging.WARNING, logging.DEBUG, logging.WARNING]
    assert [summarise(r) for r in jsonlines_text_to_objects(capsys.readouterr().out)] == [
        FALLBACK_RECORD,
        FALLBACK_RECORD,
    ]


def test_init_logger_pass_adopts_dlt_logger_on_skipped_call() -> None:
    """A skipped call still routes the dlt logger through root if dlt has reset it."""
    init_logger(settings_for(None))
    root = logging.getLogger()
    handlers = list(root.handlers)
    dlt_logger = logging.getLogger("dlt")
    dlt_logger.addHandler(logging.NullHandler())
    dlt_logger.propagate = False

    init_logger(settings_for(None))

    assert root.handlers == handlers
    assert dlt_logger.handlers == []
    assert dlt_logger.propagate is True


def test_init_logger_pass_keeps_caplog_handler_on_skipped_call(caplog: pytest.LogCaptureFixture) -> None:
    """A skipped call does not remove handlers added after the config was applied."""
    init_logger(settings_for(None))
    logging.getLogger().addHandler(caplog.handler)

    init_logger(settings_for(None))
    logging.getLogger("cdm_test.caplog").warning("still captured")

    assert caplog.messages == ["still captured"]


@pytest.mark.parametrize(
    ("file_name", "content", "expected", "match"),
    [
        pytest.param("missing.yaml", None, FileNotFoundError, None, id="missing_file"),
        pytest.param("logging.toml", "x = 1", ValueError, "Unsupported config file format: logging.toml", id="toml"),
        pytest.param("empty.yaml", "", TypeError, "must be a mapping", id="empty_yaml"),
        pytest.param("invalid.yaml", ": : invalid: yaml: : :", yaml.YAMLError, None, id="invalid_yaml"),
    ],
)
def test_init_logger_fail_unloadable_config(
    tmp_path: Path, file_name: str, content: str | None, expected: type[Exception], match: str | None
) -> None:
    """A config that can't be loaded raises, leaves root untouched, and is retried on the next call."""
    path = tmp_path / file_name
    if content is not None:
        path.write_text(content, encoding="utf-8")
    root = logging.getLogger()
    handlers, level = list(root.handlers), root.level

    for _ in range(2):
        with pytest.raises(expected, match=match):
            init_logger(settings_for(path))

    assert root.handlers == handlers
    assert root.level == level


def test_init_logger_fail_config_rejected_by_dictconfig(tmp_path: Path) -> None:
    """A config that dictConfig rejects raises every time, rather than being recorded as applied."""
    config = make_config()
    config["handlers"]["out"]["class"] = "logging.NoSuchHandler"
    path = write_json_or_yaml(tmp_path / "logging.yaml", config)

    for _ in range(2):
        with pytest.raises(ValueError, match=re.escape("Unable to configure handler 'out'")):
            init_logger(settings_for(path))


def test_init_logger_pass_formats_all_output_as_json(capsys: pytest.CaptureFixture[str]) -> None:
    """After init, app, dlt and exception messages are written to stdout as JSON, and DEBUG is filtered out."""
    dlt_logger = logging.getLogger("dlt")
    dlt_logger.disabled = False
    dlt_logger.setLevel(logging.NOTSET)
    init_logger(settings_for(None))
    capsys.readouterr()

    log = logging.getLogger("cdm_test.output")
    log.info("hello %s", "world")
    log.debug("not shown")
    dlt_logger.warning("from dlt")
    try:
        err_msg = "Oh crap!"
        raise ValueError(err_msg)  # noqa: TRY301
    except ValueError:
        log.exception("failed")

    records = jsonlines_text_to_objects(capsys.readouterr().out)
    assert [summarise(r) for r in records] == [
        {"type": "log", "logger": "dlt", "level": "WARNING", "msg": "from dlt"},
        {"type": "log", "logger": "cdm_test.output", "level": "ERROR", "msg": "failed"},
    ]
    assert "ValueError: Oh crap!" in records[-1]["exc_info"]


@pytest.mark.parametrize(
    "exc",
    [
        pytest.param(ValueError("bad value"), id="value_error"),
        pytest.param(SystemExit(1), id="system_exit"),
    ],
)
def test_log_exception_pass_logs_exception(
    caplog: pytest.LogCaptureFixture, default_excepthook_calls: list[ExcInfo], exc: BaseException
) -> None:
    """Exceptions other than KeyboardInterrupt are logged on the root logger with their exc_info."""
    exc_info = raised(exc)

    log_exception(*exc_info)

    assert [(r.name, r.levelno, r.getMessage(), r.exc_info) for r in caplog.records] == [
        ("root", logging.ERROR, "Uncaught exception", exc_info)
    ]
    assert default_excepthook_calls == []


class CustomInterrupt(KeyboardInterrupt):
    """KeyboardInterrupt subclass for testing."""


@pytest.mark.parametrize(
    "exc",
    [pytest.param(KeyboardInterrupt(), id="keyboard_interrupt"), pytest.param(CustomInterrupt(), id="subclass")],
)
def test_log_exception_pass_delegates_keyboard_interrupt(
    caplog: pytest.LogCaptureFixture, default_excepthook_calls: list[ExcInfo], exc: BaseException
) -> None:
    """KeyboardInterrupt and its subclasses go to the default excepthook and are not logged."""
    exc_info = raised(exc)

    log_exception(*exc_info)

    assert default_excepthook_calls == [exc_info]
    assert caplog.records == []


def test_log_exception_pass_installed_as_excepthook() -> None:
    """Importing the module installs log_exception as sys.excepthook."""
    assert sys.excepthook is log_exception
