"""Tests for pipelines.ncbi_ftp_promote — settings, pipeline orchestration, CLI."""

from pathlib import Path, PurePosixPath
from typing import Any, cast

import pytest
from pydantic import ValidationError

from cdm_data_loaders.pipelines.ncbi_ftp_promote import (
    DEFAULT_DESTINATION_BUCKET,
    DEFAULT_DESTINATION_PREFIX,
    DEFAULT_STAGING_BUCKET,
    DEFAULT_TRANSFER_MANIFEST_FILE,
    PromoteSettings,
)
from tests.conftest import _generate_dlt_config

_DEFAULT_STAGING_PATH: PurePosixPath = PurePosixPath("staging") / "run1"


def make_settings(**kwargs: str | int | bool | Path | PurePosixPath | None) -> PromoteSettings:
    """Generate a validated PromoteSettings object with a required staging_path default."""
    kwargs.setdefault("staging_path", _DEFAULT_STAGING_PATH)
    return PromoteSettings(**kwargs)  # pyright: ignore[reportArgumentType]


# Settings defaults / all params — (id, kwargs, expected resolved values)
_DEFAULT_RESOLVED = {
    "staging_bucket": DEFAULT_STAGING_BUCKET,
    "destination_bucket": DEFAULT_DESTINATION_BUCKET,
    "destination_path": DEFAULT_DESTINATION_PREFIX,
    "removed_manifest_path": None,
    "updated_manifest_path": None,
    "transfer_manifest_path": _DEFAULT_STAGING_PATH / DEFAULT_TRANSFER_MANIFEST_FILE,
    "dry_run": False,
}

_ALL_PARAMS_RESOLVED: dict[str, Any] = {
    "staging_bucket": PurePosixPath("my-staging-bucket"),
    "destination_bucket": PurePosixPath("my-dest-bucket"),
    "staging_path": PurePosixPath("staging") / "run42",
    "destination_path": PurePosixPath("warehouse") / "ncbi",
    "dry_run": True,
}


@pytest.mark.parametrize("case", [pytest.param(_DEFAULT_RESOLVED, id="defaults")])
def test_promote_settings_resolve_expected_values(case: dict) -> None:
    """Default PromoteSettings fields resolve to the expected values."""
    s = make_settings()
    for field, expected_value in case.items():
        assert getattr(s, field) == expected_value


def test_promote_settings_all_params_resolve_expected_values(tmp_path: Path) -> None:
    """PromoteSettings fields are correctly set when all params are provided."""
    removed = tmp_path / "removed.txt"
    updated = tmp_path / "updated.txt"
    transfer = PurePosixPath("staging") / "run42" / "transfer_manifest.txt"
    kwargs = {
        **_ALL_PARAMS_RESOLVED,
        "removed_manifest": removed,
        "updated_manifest": updated,
        "transfer_manifest": transfer,
    }

    s = make_settings(**kwargs)

    expected = {
        **_ALL_PARAMS_RESOLVED,
        "removed_manifest_path": removed,
        "updated_manifest_path": updated,
        "transfer_manifest_path": transfer,
    }
    for field, expected_value in expected.items():
        assert getattr(s, field) == expected_value


# Settings aliases — (alias kwargs, expected resolved values)
_ALIAS_CASES = [
    pytest.param(({"s": PurePosixPath("staging") / "runX"}, "staging_path"), id="staging_path_alias_s"),
    pytest.param(
        ({"destination_path": PurePosixPath("warehouse") / "custom"}, "destination_path"), id="destination_path_alias"
    ),
    pytest.param(({"r": Path("removed.txt")}, "removed_manifest_path"), id="removed_manifest_alias_r"),
    pytest.param(({"u": Path("updated.txt")}, "updated_manifest_path"), id="updated_manifest_alias_u"),
    pytest.param(
        ({"t": PurePosixPath("staging") / "run1" / "manifest.txt"}, "transfer_manifest_path"),
        id="transfer_manifest_alias_t",
    ),
    pytest.param(({"staging_bucket": PurePosixPath("alt-staging")}, "staging_bucket"), id="staging_bucket_alias"),
    pytest.param(
        ({"destination_bucket": PurePosixPath("alt-dest")}, "destination_bucket"), id="destination_bucket_alias"
    ),
]


@pytest.mark.parametrize("case", _ALIAS_CASES)
def test_promote_settings_alias_resolves_to_expected_field(case: tuple) -> None:
    """Each CLI alias resolves to the expected settings field."""
    kwargs, field = case
    s = make_settings(**kwargs)
    assert getattr(s, field) == next(iter(kwargs.values()))


# Settings validation
def test_promote_settings_staging_path_required() -> None:
    """Omitting staging_path raises ValidationError."""
    settings_ctor = cast("Any", PromoteSettings)
    with pytest.raises((ValidationError, Exception)):
        settings_ctor(_cli_parse_args=[], dlt_config=_generate_dlt_config())


def test_promote_settings_transfer_manifest_path_can_be_none() -> None:
    """transfer_manifest_path can be explicitly set to None."""
    s = make_settings(transfer_manifest=None)
    assert s.transfer_manifest_path is None
