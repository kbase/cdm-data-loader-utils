"""Tests for pipelines.ncbi_ftp_promote — settings, pipeline orchestration, CLI."""

from pathlib import Path, PurePosixPath
from typing import Any, cast
from unittest.mock import patch

import pytest
from pydantic import ValidationError

from cdm_data_loaders.pipelines.ncbi_ftp_promote import (
    DEFAULT_DESTINATION_BUCKET,
    DEFAULT_DESTINATION_PREFIX,
    DEFAULT_STAGING_BUCKET,
    DEFAULT_TRANSFER_MANIFEST_FILE,
    PromoteSettings,
    run_promote,
)
from tests.conftest import _generate_dlt_config

_DEFAULT_STAGING_PATH: PurePosixPath = PurePosixPath("staging") / "run1"


def make_settings(**kwargs: str | int | bool | Path | PurePosixPath | None) -> PromoteSettings:
    """Generate a validated PromoteSettings object with a required staging_path default."""
    settings_ctor = cast("Any", PromoteSettings)
    kwargs.setdefault("staging_path", _DEFAULT_STAGING_PATH)
    return settings_ctor(_cli_parse_args=[], **kwargs)


# Settings defaults


class TestPromoteSettingsDefaults:
    """Test default settings."""

    def test_staging_bucket_default(self) -> None:
        """Verify default staging_bucket matches DEFAULT_STAGING_BUCKET constant."""
        s = make_settings()
        assert s.staging_bucket == DEFAULT_STAGING_BUCKET

    def test_destination_bucket_default(self) -> None:
        """Verify default destination_bucket matches DEFAULT_DESTINATION_BUCKET constant."""
        s = make_settings()
        assert s.destination_bucket == DEFAULT_DESTINATION_BUCKET

    def test_destination_path_default(self) -> None:
        """Verify default destination_path matches DEFAULT_DESTINATION_PREFIX constant."""
        s = make_settings()
        assert s.destination_path == DEFAULT_DESTINATION_PREFIX

    def test_removed_manifest_default_none(self) -> None:
        """Verify default removed_manifest_path is None."""
        s = make_settings()
        assert s.removed_manifest_path is None

    def test_updated_manifest_default_none(self) -> None:
        """Verify default updated_manifest_path is None."""
        s = make_settings()
        assert s.updated_manifest_path is None

    def test_transfer_manifest_default(self) -> None:
        """Verify default transfer_manifest_path is derived from staging_path."""
        s = make_settings()
        assert s.transfer_manifest_path == _DEFAULT_STAGING_PATH / DEFAULT_TRANSFER_MANIFEST_FILE

    def test_dry_run_default_false(self) -> None:
        """Verify default dry_run is False."""
        s = make_settings()
        assert s.dry_run is False


# Settings all params


class TestPromoteSettingsAllParams:
    """Test with all params explicitly set."""

    def test_all_params(self, tmp_path: Path) -> None:
        """Verify all parameters are correctly set when provided."""
        staging = PurePosixPath("my-staging-bucket")
        dest = PurePosixPath("my-dest-bucket")
        staging_path = PurePosixPath("staging") / "run42"
        destination_path = PurePosixPath("warehouse") / "ncbi"
        removed = tmp_path / "removed.txt"
        updated = tmp_path / "updated.txt"
        transfer = PurePosixPath("staging") / "run42" / "transfer_manifest.txt"

        s = make_settings(
            staging_bucket=staging,
            destination_bucket=dest,
            staging_path=staging_path,
            destination_path=destination_path,
            removed_manifest=removed,
            updated_manifest=updated,
            transfer_manifest=transfer,
            dry_run=True,
        )

        assert s.staging_bucket == staging
        assert s.destination_bucket == dest
        assert s.staging_path == staging_path
        assert s.destination_path == destination_path
        assert s.removed_manifest_path == removed
        assert s.updated_manifest_path == updated
        assert s.transfer_manifest_path == transfer
        assert s.dry_run is True


# Settings aliases


class TestPromoteSettingsAliases:
    """Test CLI alias resolution."""

    def test_staging_path_alias_s(self) -> None:
        """Verify 's' alias resolves to staging_path."""
        path = PurePosixPath("staging") / "runX"
        s = make_settings(s=path)
        assert s.staging_path == path

    def test_destination_path_alias_destination_path(self) -> None:
        """Verify 'destination_path' alias resolves to destination_path."""
        path = PurePosixPath("warehouse") / "custom"
        s = make_settings(destination_path=path)
        assert s.destination_path == path

    def test_removed_manifest_alias_r(self, tmp_path: Path) -> None:
        """Verify 'r' alias resolves to removed_manifest_path."""
        p = tmp_path / "removed.txt"
        s = make_settings(r=p)
        assert s.removed_manifest_path == p

    def test_updated_manifest_alias_u(self, tmp_path: Path) -> None:
        """Verify 'u' alias resolves to updated_manifest_path."""
        p = tmp_path / "updated.txt"
        s = make_settings(u=p)
        assert s.updated_manifest_path == p

    def test_transfer_manifest_alias_t(self) -> None:
        """Verify 't' alias resolves to transfer_manifest_path."""
        p = PurePosixPath("staging") / "run1" / "manifest.txt"
        s = make_settings(t=p)
        assert s.transfer_manifest_path == p

    def test_staging_bucket_alias(self) -> None:
        """Verify 'staging_bucket' alias resolves to staging_bucket."""
        bucket = PurePosixPath("alt-staging")
        s = make_settings(staging_bucket=bucket)
        assert s.staging_bucket == bucket

    def test_destination_bucket_alias(self) -> None:
        """Verify 'destination_bucket' alias resolves to destination_bucket."""
        bucket = PurePosixPath("alt-dest")
        s = make_settings(destination_bucket=bucket)
        assert s.destination_bucket == bucket


# Settings validation


class TestPromoteSettingsValidation:
    """Test validation constraints."""

    def test_staging_path_required(self) -> None:
        """Verify omitting staging_path raises ValidationError."""
        settings_ctor = cast("Any", PromoteSettings)
        with pytest.raises((ValidationError, Exception)):
            settings_ctor(_cli_parse_args=[], dlt_config=_generate_dlt_config())

    def test_transfer_manifest_path_can_be_none(self) -> None:
        """Verify transfer_manifest_path can be explicitly set to None."""
        s = make_settings(transfer_manifest=None)
        assert s.transfer_manifest_path is None


# run_promote


_MOCK_REPORT_SUCCESS: dict[str, Any] = {
    "timestamp": "2026-01-01T00:00:00+00:00",
    "promoted": 5,
    "archived": 0,
    "failed": 0,
    "dry_run": False,
}

_MOCK_REPORT_WITH_FAILURES: dict[str, Any] = {
    **_MOCK_REPORT_SUCCESS,
    "promoted": 3,
    "failed": 2,
}


class TestRunPromote:
    """Test run_promote orchestration."""

    def test_calls_promote_from_s3_with_correct_args(self, tmp_path: Path) -> None:
        """Verify run_promote passes all PromoteSettings fields to promote_from_s3."""
        staging_path = PurePosixPath("staging") / "run1"
        dest_path = PurePosixPath("warehouse") / "ncbi"
        removed = tmp_path / "removed.txt"
        updated = tmp_path / "updated.txt"
        transfer = PurePosixPath("staging") / "run1" / "transfer_manifest.txt"
        config = make_settings(
            staging_bucket=PurePosixPath("my-staging"),
            destination_bucket=PurePosixPath("my-dest"),
            staging_path=staging_path,
            destination_path=dest_path,
            removed_manifest=removed,
            updated_manifest=updated,
            transfer_manifest=transfer,
            dry_run=True,
        )

        with patch(
            "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
            return_value=_MOCK_REPORT_SUCCESS,
        ) as mock_promote:
            run_promote(config)

        mock_promote.assert_called_once_with(
            staging_bucket=PurePosixPath("my-staging"),
            staging_key_prefix=staging_path,
            lakehouse_bucket=PurePosixPath("my-dest"),
            lakehouse_key_prefix=dest_path,
            removed_manifest_path=removed,
            updated_manifest_path=updated,
            manifest_s3_key=transfer,
            dry_run=True,
        )

    def test_no_error_on_zero_failures(self) -> None:
        """Verify run_promote does not raise when promote_from_s3 reports zero failures."""
        config = make_settings()
        with patch(
            "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
            return_value=_MOCK_REPORT_SUCCESS,
        ):
            run_promote(config)  # should not raise

    def test_raises_runtime_error_on_failures(self) -> None:
        """Verify run_promote raises RuntimeError when promote_from_s3 reports failures."""
        config = make_settings()
        with (
            patch(
                "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
                return_value=_MOCK_REPORT_WITH_FAILURES,
            ),
            pytest.raises(RuntimeError, match="2 failures"),
        ):
            run_promote(config)

    def test_dry_run_forwarded(self) -> None:
        """Verify dry_run=True is forwarded to promote_from_s3."""
        config = make_settings(dry_run=True)
        with patch(
            "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
            return_value=_MOCK_REPORT_SUCCESS,
        ) as mock_promote:
            run_promote(config)

        _, kwargs = mock_promote.call_args
        assert kwargs["dry_run"] is True

    def test_transfer_manifest_none_forwarded(self) -> None:
        """Verify transfer_manifest_path=None is forwarded to promote_from_s3 as manifest_s3_key=None."""
        config = make_settings(transfer_manifest=None)
        with patch(
            "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
            return_value=_MOCK_REPORT_SUCCESS,
        ) as mock_promote:
            run_promote(config)

        _, kwargs = mock_promote.call_args
        assert kwargs["manifest_s3_key"] is None
