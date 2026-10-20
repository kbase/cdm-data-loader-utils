"""Tests for pipelines.ncbi_ftp_promote — settings, pipeline orchestration, CLI."""

from pathlib import Path, PurePosixPath
from typing import Any
from unittest.mock import patch

import pytest

from cdm_data_loaders.pipelines.ncbi_ftp_promote import (
    run_promote,
)
from tests.cdm_data_loaders.pipelines.ncbi_ftp_promote.test_ncbi_ftp_promote_settings import make_settings

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


def test_run_promote_calls_promote_from_s3_with_correct_args(tmp_path: Path) -> None:
    """run_promote passes all PromoteSettings fields to promote_from_s3."""
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


def test_run_promote_no_error_on_zero_failures() -> None:
    """run_promote does not raise when promote_from_s3 reports zero failures."""
    config = make_settings()
    with patch(
        "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
        return_value=_MOCK_REPORT_SUCCESS,
    ):
        run_promote(config)  # should not raise


def test_run_promote_raises_runtime_error_on_failures() -> None:
    """run_promote raises RuntimeError when promote_from_s3 reports failures."""
    config = make_settings()
    with (
        patch(
            "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
            return_value=_MOCK_REPORT_WITH_FAILURES,
        ),
        pytest.raises(RuntimeError, match="2 failures"),
    ):
        run_promote(config)


def test_run_promote_dry_run_forwarded() -> None:
    """dry_run=True is forwarded to promote_from_s3."""
    config = make_settings(dry_run=True)
    with patch(
        "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
        return_value=_MOCK_REPORT_SUCCESS,
    ) as mock_promote:
        run_promote(config)

    _, kwargs = mock_promote.call_args
    assert kwargs["dry_run"] is True


def test_run_promote_transfer_manifest_none_forwarded() -> None:
    """transfer_manifest_path=None is forwarded to promote_from_s3 as manifest_s3_key=None."""
    config = make_settings(transfer_manifest=None)
    with patch(
        "cdm_data_loaders.pipelines.ncbi_ftp_promote.promote_from_s3",
        return_value=_MOCK_REPORT_SUCCESS,
    ) as mock_promote:
        run_promote(config)

    _, kwargs = mock_promote.call_args
    assert kwargs["manifest_s3_key"] is None
