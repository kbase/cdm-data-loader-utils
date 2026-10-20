# cdm_data_loaders/utils/cdm_logger.py
"""Provides structured logging with contextual metadata for CDM data import pipelines."""

import copy
import json
import logging
import logging.config
import sys
from pathlib import Path
from types import TracebackType
from typing import Any

import yaml
from dlt.common.runtime.json_logging import JSONLogFormatter

from cdm_data_loaders.core.settings import LoggerSettings

DLT_LOGGER_NAME = "dlt"
JSON_FORMATTER_NAME = "dlt_json"
FALLBACK_KEY = ("<built-in>", 0)
applied_key: tuple[str, int] | None = None


# There is no default logger config file (see cdm_data_loaders.core.fields.DEFAULTS[LOG_CONFIG_FILE])
# so the logger falls back to this config if an external config file is not found.
def _default_config() -> dict[str, Any]:
    """Default config for the logger: logs WARNING and above to STDOUT."""
    return {
        "version": 1,
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "stream": "ext://sys.stdout",
            },
        },
        "root": {"level": "WARNING", "handlers": ["console"]},
        "disable_existing_loggers": False,
    }


def _load_config_from_path(path: str | Path) -> dict[str, Any]:
    """Attempt to load and parse a logging config from the given path.

    :param path: path to a JSON or YAML logging config file
    :type  path: str | Path
    :return: parsed config dict
    :rtype : dict[str, Any]
    """
    path = Path(path)

    suffix = path.suffix.lower()

    if suffix in (".yml", ".yaml"):
        with path.open() as f:
            config = yaml.safe_load(f)
    elif suffix == ".json":
        with path.open() as f:
            config = json.load(f)
    else:
        err_msg = f"Unsupported config file format: {path.name}"
        raise ValueError(err_msg)

    if not isinstance(config, dict):
        err_msg = f"Logging config in {path.name} must be a mapping, not {type(config).__name__}"
        raise TypeError(err_msg)

    return config


def _enforce_json_formatter(config: dict[str, Any]) -> dict[str, Any]:
    """Return a copy of ``config`` in which every handler uses dlt's JSONLogFormatter.

    :param config: logger config
    :type  config: dict[str, Any]
    :return: config copy with formatters set to dlt's JSONLogFormatter
    :rtype : dict[str, Any]
    """
    config = copy.deepcopy(config)
    # dictConfig accepts a callable for the "()" factory key
    config.setdefault("formatters", {})[JSON_FORMATTER_NAME] = {"()": JSONLogFormatter}
    for handler_cfg in config.get("handlers", {}).values():
        handler_cfg["formatter"] = JSON_FORMATTER_NAME
    config["disable_existing_loggers"] = False
    return config


def adopt_dlt_logger() -> None:
    """Route dlt's logger through the root logger's (JSON-formatted) handlers.

    Safe to call repeatedly.
    """
    dlt_logger = logging.getLogger(DLT_LOGGER_NAME)
    # pass messages from the dlt logger up to the root
    dlt_logger.propagate = True
    # clear the dlt logger handlers and let everything be handled by the root logger instead.
    dlt_logger.handlers.clear()


def init_logger(settings: LoggerSettings | None = None, *, force: bool = False) -> None:
    """Initialise JSON logging; repeat calls are skipped unless the config file changes.

    :param settings: settings for the logger, defaults to None
    :type  settings: LoggerSettings | None, optional
    :param force: force the logger to be reconfigured, defaults to False
    :type  force: bool, optional
    """
    global applied_key  # noqa: PLW0603

    if settings is None:
        settings = LoggerSettings()  # pyright: ignore[reportCallIssue]

    path = Path(settings.log_config_file).resolve() if settings.log_config_file else None
    key = (str(path), path.stat().st_mtime_ns) if path else FALLBACK_KEY

    if force or key != applied_key:
        config = _load_config_from_path(path) if path else _default_config()
        logging.config.dictConfig(_enforce_json_formatter(config))
        applied_key = key
        if path is None:
            logging.getLogger(__name__).warning("No logging config file found. Falling back to built-in config.")

    adopt_dlt_logger()


def log_exception(exc_type: type[BaseException], exc_value: BaseException, exc_traceback: TracebackType | None) -> None:
    """Log uncaught exceptions using the exception handler."""
    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logging.getLogger().error("Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback))


sys.excepthook = log_exception
