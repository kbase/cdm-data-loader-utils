"""Shared fixtures for pipelines tests."""

from typing import Any, Final

import dlt
import pytest
from frozendict import frozendict

from cdm_data_loaders.core.fields import (
    BUFFER_SIZE,
    DEV_MODE,
    INPUT_DIR,
    LOG_CONFIG_FILE,
    LOG_INTERVAL,
    OUTPUT_DIR,
    START_AT,
    USE_DESTINATION,
    USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
)
from cdm_data_loaders.core.settings import (
    DEFAULT_BATCH_FILE_SETTINGS,
    DEFAULT_CTS_SETTINGS,
    CtsSettings,
)
from tests.conftest import TEST_DLT_CONFIG

CASSETTES_DIR = "tests/cassettes"

START_AT_VALUE: Final[int] = 50
START_AT_STRING: Final[str] = "50"
TEST_LOG_CONFIG_FILE: Final[str] = "log_conf.json"

DESTINATION_TO_OUTPUT = frozendict(
    {
        "local_fs": TEST_DLT_CONFIG["destination.local_fs.bucket_url"],
        "s3": TEST_DLT_CONFIG["destination.s3.bucket_url"],
    }
)

DESTINATION_OUTPUT: Final[str] = DESTINATION_TO_OUTPUT[DEFAULT_CTS_SETTINGS["use_destination"]]

DEFAULT_CTS_SETTINGS_RECONCILED = frozendict(
    {
        **DEFAULT_CTS_SETTINGS,
        OUTPUT_DIR: DESTINATION_OUTPUT,
        "raw_data_dir": f"{DESTINATION_OUTPUT}/raw_data",
        "pipeline_dir": None,
    }
)

DEFAULT_BATCH_FILE_SETTINGS_RECONCILED = frozendict({**DEFAULT_BATCH_FILE_SETTINGS, **DEFAULT_CTS_SETTINGS_RECONCILED})

TEST_CTS_SETTINGS = frozendict(
    {
        DEV_MODE: "false",
        INPUT_DIR: "/dir/path",
        LOG_CONFIG_FILE: "some/path",
        OUTPUT_DIR: "/some/dir",
        USE_DESTINATION: "local_fs",
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: "true",
    }
)

TEST_CTS_SETTINGS_RECONCILED = frozendict(
    {
        **TEST_CTS_SETTINGS,
        DEV_MODE: False,
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: True,
        "pipeline_dir": "/some/dir/.dlt_conf",
        "raw_data_dir": "/some/dir/raw_data",
    }
)

TEST_BATCH_FILE_SETTINGS = frozendict(
    **TEST_CTS_SETTINGS,
    start_at=START_AT_STRING,
    log_interval=5000,
    buffer_size=25,
)

TEST_BATCH_FILE_SETTINGS_RECONCILED = frozendict(
    {
        **TEST_CTS_SETTINGS_RECONCILED,
        BUFFER_SIZE: 25,
        LOG_INTERVAL: 5000,
        START_AT: START_AT_VALUE,
        "pipeline_dir": "/some/dir/.dlt_conf",
        "raw_data_dir": "/some/dir/raw_data",
    }
)


def make_settings(
    settings_cls: type[CtsSettings],
    dlt_config: dict[str, Any] | None = None,
    kwargs: dict[str, Any] | frozendict[str, Any] | None = None,
) -> CtsSettings:  # CtsSettings | BatchedFileInputSettings | NcbiRestApiSettings | AtbSettings:
    """Generate a validated Settings object."""
    with pytest.MonkeyPatch.context() as m:
        m.setattr(dlt, "config", dlt_config)
        return settings_cls(**(kwargs or {}))  # pyright: ignore[reportArgumentType]


def make_settings_autofill_config(
    settings_cls: type[CtsSettings],
    kwargs: dict[str, Any] | frozendict[str, Any] | None = None,
) -> (
    CtsSettings
):  # CtsSettings | BatchedFileInputSettings | NcbiRestApiSettings | AtbSettings | UniProtSettings | UnirefSettings:
    """Generate a validated Settings object, supplying the dlt_config if necessary."""
    return settings_cls(**(kwargs or {}))  # pyright: ignore[reportArgumentType]


def check_settings(
    settings_object: CtsSettings,
    expected: dict[str, Any] | frozendict[str, Any],
) -> None:
    """Check that the settings object has the expected values."""
    assert settings_object._dlt_config is not None  # noqa: SLF001
    assert settings_object.model_dump() == expected

    # make sure we have both raw_data_dir and pipeline_dir
    assert "raw_data_dir" in expected
    assert "pipeline_dir" in expected
    for attr, value in expected.items():
        assert getattr(settings_object, attr) == value
