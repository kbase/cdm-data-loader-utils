"""Shared fixtures and functions for cdm_data_loaders tests."""

from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, Final

import pytest
from frozendict import frozendict
from pydantic_settings import BaseSettings, CliApp

from cdm_data_loaders.core.fields import (
    BUFFER_SIZE,
    DEFAULTS,
    DLT_DEV_MODE,
    INPUT_DIR,
    LOG_CONFIG_FILE,
    LOG_INTERVAL,
    OUTPUT_DIR,
    USE_DESTINATION,
    USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
)
from tests.conftest import LOCAL_FS


def check_settings(
    settings_object: BaseSettings,
    expected: dict[str, Any] | frozendict[str, Any],
) -> None:
    """Check that the settings object has the expected values."""
    assert settings_object.model_dump() == expected


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


def local_io_expected(input_dir: Path, output_dir: Path) -> dict[str, Any]:
    """Expected model_dump of LocalInputOutputSettings with default logging."""
    return {
        "log_config_file": DEFAULTS[LOG_CONFIG_FILE],
        "input_dir": input_dir.resolve(),
        "output_dir": output_dir.resolve(),
    }


def cts_expected(input_dir: Path, output_dir: Path, **overrides: Any) -> dict[str, Any]:  # noqa: ANN401
    """Expected model_dump of CtsSettings: defaults, the given dirs, and any overrides."""
    return {
        **local_io_expected(input_dir, output_dir),
        BUFFER_SIZE: DEFAULTS[BUFFER_SIZE],
        DLT_DEV_MODE: DEFAULTS[DLT_DEV_MODE],
        LOG_INTERVAL: DEFAULTS[LOG_INTERVAL],
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: DEFAULTS[USE_OUTPUT_DIR_FOR_PIPELINE_METADATA],
        "pipeline_dir": None,
        **overrides,
    }


@pytest.fixture
def input_dir(tmp_path: Path) -> Path:
    """Existing input directory."""
    path = tmp_path / "input"
    path.mkdir(exist_ok=True)
    return path.resolve()


@pytest.fixture
def output_dir(tmp_path: Path) -> Path:
    """Output directory that does not exist yet."""
    return (tmp_path / "output").resolve()


@pytest.fixture
def test_defaults(input_dir: Path, output_dir: Path) -> dict[str, Any]:
    """Defaults for required fields, used for testing."""
    return {INPUT_DIR: input_dir, OUTPUT_DIR: output_dir, USE_DESTINATION: LOCAL_FS}


@pytest.fixture
def test_cts_settings(input_dir: Path, output_dir: Path) -> frozendict[str, Any]:
    """Settings dictionary with all CtsSettings fields specified."""
    return frozendict(
        {
            DLT_DEV_MODE: "false",
            LOG_CONFIG_FILE: "some/path",
            INPUT_DIR: input_dir,
            OUTPUT_DIR: output_dir,
            USE_DESTINATION: LOCAL_FS,
            USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: "true",
            LOG_INTERVAL: "5000",
            BUFFER_SIZE: 25,
        }
    )


@pytest.fixture
def test_cts_settings_reconciled(test_cts_settings: frozendict[str, Any]) -> frozendict[str, Any]:
    """Reconciled CtsSettings dictionary matching model_dump()."""
    return frozendict(
        {
            **test_cts_settings,
            DLT_DEV_MODE: False,
            USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: True,
            LOG_INTERVAL: 5000,
            "pipeline_dir": test_cts_settings[OUTPUT_DIR] / ".dlt_conf",
        }
    )
