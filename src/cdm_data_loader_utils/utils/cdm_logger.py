"""
Provides structured logging with contextual metadata for CDM data import pipelines.
"""

import logging
import logging.handlers
import os
import sys
from pathlib import Path

DEFAULT_LOGGER_NAME = "cdm_data_loader"

GENERIC_ERROR_MESSAGE = "An error of unknown origin occurred."

LOG_FILENAME = "cdm_data_loader.log"
MAX_LOG_FILE_SIZE = 2**30  # 1 GiB
MAX_LOG_BACKUPS = 5

__LOGGER = None

# TODO: adopt logging config, set just once
LOGGING_CONFIG = {
    "root": {"name": "cdm_data_loader", "level": "INFO", "handlers": ["console", "file"]},
    "version": 1,
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "json",
            "level": "INFO",
            "stream": "ext://sys.stdout",
        },
        "file": {
            "class": "logging.handlers.RotatingFileHandler",
            "formatter": "json",
            "filename": LOG_FILENAME,
            "maxBytes": MAX_LOG_FILE_SIZE,
            "backupCount": MAX_LOG_BACKUPS,
        },
    },
    "formatters": {
        "json": {
            "format": '{"time": "%(asctime)s", "level": "%(levelname)s", "module": "%(module)s", "msg": "%(message)s"}'
        }
    },
}


def get_cdm_logger(
    logger_name: str | None = None, log_level: str | None = None, log_dir: str | None = None
) -> logging.Logger:
    """Retrieve the logger, initialising it if necessary.

    If the logger name is not set, the default name "cdm_data_loader" will be used.

    :param logger_name: name for the logger, defaults to None
    :type logger_name: str | None, optional
    :param log_level: logger level, defaults to None
    :type log_level: str | None, optional
    :param log_dir: directory to save log files to, optional. If no directory is specified, logs will just be emitted to the console.
    :type log_dir: str | None
    :return: initialised logger
    :rtype: logging.Logger
    """
    global __LOGGER
    if not __LOGGER:
        __LOGGER = init_logger(logger_name, log_level, log_dir)
    return __LOGGER


def init_logger(
    logger_name: str | None = None, log_level: str | None = None, log_dir: str | None = None
) -> logging.Logger:
    """Initialise the logger for the module.

    If the logger name is not set, the default name "cdm_data_loader" will be used.

    :param logger_name: name for the logger, defaults to None
    :type logger_name: str | None, optional
    :param log_level: logger level, defaults to None
    :type log_level: str | None, optional
    :param log_dir: directory to save log files to, optional. If no directory is specified, logs will just be emitted to the console.
    :type log_dir: str | None
    :return: initialised logger
    :rtype: logging.Logger
    """
    if not logger_name:
        logger_name = DEFAULT_LOGGER_NAME

    # Always get the same logger by name
    logger = logging.getLogger(logger_name)

    # Determine log level (argument > env var > default)
    effective_log_level = (log_level or os.getenv("LOG_LEVEL", "INFO")).upper()
    logger.setLevel(getattr(logging, effective_log_level, logging.DEBUG))

    # JSON-style structured formatter
    formatter = logging.Formatter(
        '{"time": "%(asctime)s", "level": "%(levelname)s", "module": "%(module)s", "msg": "%(message)s"}'
    )

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    if log_dir:
        log_dir_path = Path(log_dir)
        if not log_dir_path.exists() and log_dir_path.is_dir():
            msg = f"{log_dir} does not exist or is not a directory."
            raise FileNotFoundError(msg)
        # Add the log message handler to the logger
        file_handler = logging.handlers.RotatingFileHandler(
            LOG_FILENAME, maxBytes=MAX_LOG_FILE_SIZE, backupCount=MAX_LOG_BACKUPS
        )
        logger.addHandler(file_handler)
    return logger


def log_and_die(error_msg: str, error_class: type[Exception], logger_name: str | None = None) -> None:
    """Log an error message and then raise the error.

    :param error_msg: error message string
    :type error_msg: str
    :param error_class: class of error to throw
    :type error_class: type[Exception]
    :param logger_name: name of the logger to use, defaults to None
    :type logger_name: str | None, optional
    """
    logger = get_cdm_logger(logger_name)

    if not error_msg:
        logger.warning("No error supplied to log_and_die. Using generic error message.")
        error_msg = GENERIC_ERROR_MESSAGE

    if not isinstance(error_class, type) or not issubclass(error_class, BaseException):
        error_class = RuntimeError

    logger.error(error_msg)
    raise error_class(error_msg)
