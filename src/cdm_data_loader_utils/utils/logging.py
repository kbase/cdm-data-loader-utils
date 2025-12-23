"""
Provides structured logging with contextual metadata for CDM data import pipelines.
"""

import logging
import os


def get_logger(logger_name: str | None = None, log_level: str | None = None) -> logging.Logger:
    """Initialise the logger for the module.

    If the logger name is not set, the default name "cdm_data_loader" will be used.

    :param logger_name: name for the logger, defaults to None
    :type logger_name: str | None, optional
    :param log_level: logger level, defaults to None
    :type log_level: str | None, optional
    :return: initialised logger
    :rtype: logging.Logger
    """
    if not logger_name:
        logger_name = "cdm_data_loader"
    # Always get the same logger by name
    logger = logging.getLogger(logger_name)

    # Determine log level (argument > env var > default)
    effective_log_level = (log_level or os.getenv("LOG_LEVEL", "INFO")).upper()
    logger.setLevel(getattr(logging, effective_log_level, logging.DEBUG))

    return logger
