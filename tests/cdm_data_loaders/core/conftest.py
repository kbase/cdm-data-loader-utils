"""Shared fixtures and constants for core and pipelines tests."""

from collections.abc import Callable, Mapping
from typing import Any, Final

import pytest
from frozendict import frozendict
from pydantic_settings import BaseSettings, CliApp

from cdm_data_loaders.core.fields import (
    BUFFER_SIZE,
    DLT_DEV_MODE,
    INPUT_DIR,
    LOCAL_FS,
    LOG_CONFIG_FILE,
    LOG_INTERVAL,
    OUTPUT_DIR,
    S3,
    USE_DESTINATION,
    USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
)
from cdm_data_loaders.core.settings import (
    DEFAULT_CTS_SETTINGS,
    CtsSettings,
)
from tests.conftest import TEST_DLT_CONFIG

TEST_LOG_CONFIG_FILE: Final[str] = "log_conf.json"

DESTINATION_TO_OUTPUT = frozendict(
    {
        LOCAL_FS: TEST_DLT_CONFIG["destination.local_fs.bucket_url"],
        S3: TEST_DLT_CONFIG["destination.s3.bucket_url"],
    }
)

DESTINATION_OUTPUT: Final[str] = DESTINATION_TO_OUTPUT[DEFAULT_CTS_SETTINGS["use_destination"]]

DEFAULT_CTS_SETTINGS_RECONCILED = frozendict(
    {
        **DEFAULT_CTS_SETTINGS,
        OUTPUT_DIR: DESTINATION_OUTPUT,
        "output_is_local": True,
        "raw_data_dir": f"{DESTINATION_OUTPUT}/raw_data",
        "pipeline_dir": None,
    }
)

TEST_CTS_SETTINGS = frozendict(
    {
        DLT_DEV_MODE: "false",
        INPUT_DIR: "/dir/path",
        LOG_CONFIG_FILE: "some/path",
        OUTPUT_DIR: "/some/dir",
        USE_DESTINATION: LOCAL_FS,
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: "true",
        LOG_INTERVAL: "5000",
        BUFFER_SIZE: 25,
    }
)

TEST_CTS_SETTINGS_RECONCILED = frozendict(
    {
        **TEST_CTS_SETTINGS,
        DLT_DEV_MODE: False,
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: True,
        LOG_INTERVAL: 5000,
        "output_is_local": True,
        "pipeline_dir": "/some/dir/.dlt_conf",
        "raw_data_dir": "/some/dir/raw_data",
    }
)


def check_settings(
    settings_object: CtsSettings,
    expected: dict[str, Any] | frozendict[str, Any],
) -> None:
    """Check that the settings object has the expected values."""
    assert settings_object.model_dump() == expected

    # make sure we have both raw_data_dir and pipeline_dir
    assert "raw_data_dir" in expected
    assert "pipeline_dir" in expected
    for attr, value in expected.items():
        assert getattr(settings_object, attr) == value


# ways of supplying values to a settings class
INIT: Final[str] = "init"
ENV: Final[str] = "env"
CLI: Final[str] = "cli"
SETTINGS_SOURCES: Final[tuple[str, ...]] = (INIT, ENV, CLI)

type SettingsFactory = Callable[[type[BaseSettings], str, Mapping[str, str]], BaseSettings]


@pytest.fixture
def make_settings(monkeypatch: pytest.MonkeyPatch) -> SettingsFactory:
    """Build a settings object from string values supplied as init kwargs, env vars, or CLI args."""

    def _make(settings_cls: type[BaseSettings], source: str, values: Mapping[str, str]) -> BaseSettings:
        if source == INIT:
            return settings_cls(**values)  # pyright: ignore[reportArgumentType]
        if source == ENV:
            prefix = settings_cls.model_config.get("env_prefix", "")
            for name, value in values.items():
                monkeypatch.setenv(f"{prefix}{name}".upper(), value)
            return settings_cls()
        if source == CLI:
            cli_args = [arg for name, value in values.items() for arg in (f"--{name.replace('_', '-')}", value)]
            return CliApp.run(settings_cls, cli_args=cli_args)
        err_msg = f"Unknown settings source: {source!r}"
        raise ValueError(err_msg)

    return _make
