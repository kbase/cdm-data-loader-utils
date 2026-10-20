"""Tests for pipelines.ncbi_ftp_download — settings, batch orchestration, CLI."""

from pathlib import Path

import pytest
from pydantic import ValidationError

from cdm_data_loaders.core.fields import INPUT_MOUNT, OUTPUT_MOUNT
from cdm_data_loaders.ncbi_ftp.assembly import FTP_HOST
from cdm_data_loaders.pipelines.ncbi_ftp_download import (
    DownloadSettings,
)

_DEFAULT_THREADS = 4
_CUSTOM_THREADS = 8
_ALIAS_THREADS = 16
_BOUNDARY_MIN = 1
_BOUNDARY_MAX = 32
_OVER_MAX = 64
_CUSTOM_LIMIT = 100
_ALIAS_LIMIT = 50

# Settings defaults / all params / aliases

_DEFAULT_SETTINGS = {
    "manifest": Path(INPUT_MOUNT) / "transfer_manifest.txt",
    "output_dir": Path(OUTPUT_MOUNT),
    "threads": _DEFAULT_THREADS,
    "ftp_host": FTP_HOST,
    "limit": None,
}

_ALL_PARAMS_SETTINGS = {
    "manifest": Path("/") / "data" / "my_manifest.txt",
    "output_dir": Path("/") / "data" / "output",
    "threads": _CUSTOM_THREADS,
    "ftp_host": "ftp.example.com",
    "limit": _CUSTOM_LIMIT,
}

_ALIAS_SETTINGS = {
    "manifest": Path("/") / "data" / "m.txt",
    "output_dir": Path("/") / "data" / "o",
    "threads": _ALIAS_THREADS,
    "limit": _ALIAS_LIMIT,
}

# (kwargs, expected resolved values) — human-readable id per case
_SETTINGS_CASES = [
    pytest.param(({}, _DEFAULT_SETTINGS), id="defaults"),
    pytest.param((dict(_ALL_PARAMS_SETTINGS), _ALL_PARAMS_SETTINGS), id="all_params"),
    pytest.param(
        (
            {
                "m": _ALIAS_SETTINGS["manifest"],
                "output_dir": _ALIAS_SETTINGS["output_dir"],
                "t": _ALIAS_SETTINGS["threads"],
                "l": _ALIAS_SETTINGS["limit"],
            },
            _ALIAS_SETTINGS,
        ),
        id="aliases",
    ),
]

invalid_settings = [
    pytest.param({"threads": 0}, id="threads_too_low"),
    pytest.param({"threads": _OVER_MAX}, id="threads_too_high"),
    pytest.param({"limit": 0}, id="limit_must_be_positive"),
]


@pytest.mark.parametrize("case", _SETTINGS_CASES)
def test_settings_resolve_expected_values(case: tuple) -> None:
    """Defaults, explicit params, and CLI aliases all resolve to the expected settings values."""
    kwargs, expected = case
    s = DownloadSettings(**kwargs)
    for field, expected_value in expected.items():
        assert getattr(s, field) == expected_value


@pytest.mark.parametrize("kwargs", invalid_settings)
def test_settings_invalid_constraints_raise(kwargs: dict) -> None:
    """Out-of-range threads/limit values raise ValidationError."""
    with pytest.raises(ValidationError):
        DownloadSettings(**kwargs)


@pytest.mark.parametrize(
    "threads",
    [
        pytest.param(_BOUNDARY_MIN, id="threads_boundary_1"),
        pytest.param(_BOUNDARY_MAX, id="threads_boundary_32"),
    ],
)
def test_settings_threads_boundaries_accepted(threads: int) -> None:
    """Boundary thread counts (1 and 32) are accepted."""
    s = DownloadSettings(threads=threads)  # pyright: ignore[reportCallIssue]
    assert s.threads == threads
