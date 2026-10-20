"""Unit tests for transfer_config.py."""

import logging
import math

import pytest
from boto3.s3.transfer import TransferConfig

from cdm_data_loaders.utils.file_transfer.s3.transfer_config import (
    DEFAULT_MULTIPART_CHUNKSIZE,
    DEFAULT_MULTIPART_THRESHOLD,
    MB,
    S3_MAX_PARTS,
    S3_MIN_PART_SIZE,
    UNKNOWN_SIZE_CHUNKSIZE,
    build_transfer_config,
    compute_multipart_chunksize,
)

"""Pointless tests for module-level constants to get coverage up to 100%."""


def test_constants_pass_mb_is_1048576() -> None:
    """MB is defined as 1024**2 bytes."""
    assert MB == 1024**2


def test_constants_pass_s3_max_parts_matches_aws_limit() -> None:
    """S3_MAX_PARTS matches AWS's documented multipart upload part limit."""
    assert S3_MAX_PARTS == 10_000


def test_constants_pass_s3_min_part_size_is_5mb() -> None:
    """S3_MIN_PART_SIZE matches AWS's documented minimum multipart part size."""
    assert S3_MIN_PART_SIZE == 5 * MB


def test_constants_pass_defaults_match_s3transfer_defaults() -> None:
    """DEFAULT_MULTIPART_THRESHOLD and DEFAULT_MULTIPART_CHUNKSIZE both start at 8 MB, matching s3transfer's own defaults."""
    assert DEFAULT_MULTIPART_THRESHOLD == 8 * MB
    assert DEFAULT_MULTIPART_CHUNKSIZE == 8 * MB


def test_constants_pass_unknown_size_chunksize_is_64mb() -> None:
    """UNKNOWN_SIZE_CHUNKSIZE is 64 MB, used when file size cannot be determined ahead of time."""
    assert UNKNOWN_SIZE_CHUNKSIZE == 64 * MB


"""Tests for compute_multipart_chunksize"""


def test_compute_multipart_chunksize_pass_no_args_returns_default() -> None:
    """With no file_size and no requested_chunksize, the default chunksize is returned unchanged."""
    assert compute_multipart_chunksize(None) == DEFAULT_MULTIPART_CHUNKSIZE


def test_compute_multipart_chunksize_pass_file_size_none_uses_requested_chunksize() -> None:
    """When file_size is None, no scaling is attempted and the requested chunksize is used as-is."""
    assert compute_multipart_chunksize(None, 16 * MB) == 16 * MB


@pytest.mark.parametrize("file_size", [0, -1, -1000])
def test_compute_multipart_chunksize_pass_non_positive_file_size_skips_scaling(file_size: int) -> None:
    """A zero or negative file_size is treated the same as "no scaling needed", same as file_size=None."""
    assert compute_multipart_chunksize(file_size, 16 * MB) == 16 * MB


def test_compute_multipart_chunksize_pass_zero_requested_chunksize_falls_back_to_default() -> None:
    """A requested_chunksize of 0 is falsy and is treated as "not provided", falling back to the default.

    Documents current behavior: 0 is indistinguishable from None due to `requested_chunksize or DEFAULT...`.
    """
    assert compute_multipart_chunksize(None, 0) == DEFAULT_MULTIPART_CHUNKSIZE


def test_compute_multipart_chunksize_pass_small_file_uses_requested_chunksize_unchanged() -> None:
    """A file well within the part limit at the requested chunksize is not scaled up."""
    assert compute_multipart_chunksize(1000, 16 * MB) == 16 * MB


def test_compute_multipart_chunksize_pass_result_never_below_min_part_size() -> None:
    """A tiny requested_chunksize is floored at S3_MIN_PART_SIZE, since S3 rejects smaller multipart parts."""
    assert compute_multipart_chunksize(None, 1024) == S3_MIN_PART_SIZE


def test_compute_multipart_chunksize_pass_negative_requested_chunksize_floored_at_min() -> None:
    """A negative requested_chunksize (truthy, so not replaced by the default) is still floored at S3_MIN_PART_SIZE."""
    assert compute_multipart_chunksize(None, -100) == S3_MIN_PART_SIZE


def test_compute_multipart_chunksize_pass_exact_boundary_does_not_trigger_scaling() -> None:
    """A file size that exactly requires S3_MAX_PARTS parts at the given chunksize is NOT scaled up.

    The scaling condition is a strict `>`, so landing exactly on the limit is fine.
    """
    chunksize = 16 * MB
    file_size = chunksize * S3_MAX_PARTS
    assert compute_multipart_chunksize(file_size, chunksize) == chunksize


def test_compute_multipart_chunksize_pass_one_byte_over_boundary_triggers_scaling() -> None:
    """A file size just one byte over the exact-limit boundary forces the chunksize to scale up."""
    chunksize = 16 * MB
    file_size = chunksize * S3_MAX_PARTS + 1
    expected = math.ceil(file_size / S3_MAX_PARTS)
    result = compute_multipart_chunksize(file_size, chunksize)
    assert result == expected
    assert result > chunksize


def test_compute_multipart_chunksize_pass_scaling_matches_ceil_division() -> None:
    """When scaling is required, the resulting chunksize is exactly ceil(file_size / S3_MAX_PARTS)."""
    file_size = 6 * 1024**4  # 6 TB
    expected = math.ceil(file_size / S3_MAX_PARTS)
    assert compute_multipart_chunksize(file_size) == expected


def test_compute_multipart_chunksize_pass_huge_file_result_still_respects_min_part_size() -> None:
    """Even after scaling for a huge file, the result never drops below S3_MIN_PART_SIZE."""
    file_size = 6 * 1024**4  # 6 TB
    result = compute_multipart_chunksize(file_size)
    assert result >= S3_MIN_PART_SIZE


def test_compute_multipart_chunksize_pass_warns_when_scaling_occurs(caplog: pytest.LogCaptureFixture) -> None:
    """A warning is logged, including the file size and both chunksize values, when scaling is triggered."""
    chunksize = 16 * MB
    file_size = chunksize * S3_MAX_PARTS + 1
    with caplog.at_level(logging.WARNING):
        result = compute_multipart_chunksize(file_size, chunksize)

    assert len(caplog.records) == 1
    message = caplog.records[0].getMessage()
    assert str(file_size) in message
    assert str(S3_MAX_PARTS) in message
    assert str(chunksize) in message
    assert str(result) in message
    assert caplog.records[0].levelno == logging.WARNING


def test_compute_multipart_chunksize_pass_no_warning_when_scaling_not_needed(caplog: pytest.LogCaptureFixture) -> None:
    """No warning is logged for a file size that does not require scaling."""
    with caplog.at_level(logging.WARNING):
        compute_multipart_chunksize(1000, 16 * MB)
    assert len(caplog.records) == 0


def test_compute_multipart_chunksize_pass_no_warning_for_unknown_file_size(caplog: pytest.LogCaptureFixture) -> None:
    """No warning is logged when file_size is None, since scaling is never attempted."""
    with caplog.at_level(logging.WARNING):
        compute_multipart_chunksize(None, 16 * MB)
    assert len(caplog.records) == 0


def test_compute_multipart_chunksize_pass_return_type_is_int() -> None:
    """The returned chunksize is always a plain int, even after ceil-division scaling."""
    result = compute_multipart_chunksize(6 * 1024**4)
    assert isinstance(result, int)


"""Tests for build_transfer_config"""


def test_build_transfer_config_pass_returns_transferconfig_instance() -> None:
    """build_transfer_config always returns a real TransferConfig instance."""
    config = build_transfer_config(1000)
    assert isinstance(config, TransferConfig)


def test_build_transfer_config_pass_default_multipart_threshold() -> None:
    """multipart_threshold defaults to DEFAULT_MULTIPART_THRESHOLD when not overridden."""
    config = build_transfer_config(1000)
    assert config.multipart_threshold == DEFAULT_MULTIPART_THRESHOLD


def test_build_transfer_config_pass_multipart_threshold_override_respected() -> None:
    """An explicit multipart_threshold override is used instead of the default."""
    config = build_transfer_config(1000, multipart_threshold=100 * MB)
    assert config.multipart_threshold == 100 * MB


def test_build_transfer_config_pass_known_file_size_no_override_uses_computed_chunksize() -> None:
    """For a known file_size with no explicit chunksize, multipart_chunksize matches compute_multipart_chunksize(file_size, None)."""
    file_size = 500 * MB
    config = build_transfer_config(file_size)
    assert config.multipart_chunksize == compute_multipart_chunksize(file_size, None)
    assert config.multipart_chunksize == DEFAULT_MULTIPART_CHUNKSIZE


def test_build_transfer_config_pass_known_file_size_with_override_used_as_starting_point() -> None:
    """For a known file_size, an explicit multipart_chunksize override is used as the starting point for scaling."""
    file_size = 1000
    config = build_transfer_config(file_size, multipart_chunksize=32 * MB)
    assert config.multipart_chunksize == compute_multipart_chunksize(file_size, 32 * MB)
    assert config.multipart_chunksize == 32 * MB


def test_build_transfer_config_pass_known_file_size_override_below_min_is_floored() -> None:
    """A too-small multipart_chunksize override is still floored at S3_MIN_PART_SIZE for a known file_size."""
    config = build_transfer_config(1000, multipart_chunksize=1024)
    assert config.multipart_chunksize == S3_MIN_PART_SIZE


def test_build_transfer_config_pass_known_huge_file_size_scales_override() -> None:
    """A known, very large file_size scales up even an explicitly-requested chunksize override."""
    chunksize = 16 * MB
    file_size = chunksize * S3_MAX_PARTS + 1
    config = build_transfer_config(file_size, multipart_chunksize=chunksize)
    assert config.multipart_chunksize == compute_multipart_chunksize(file_size, chunksize)
    assert config.multipart_chunksize > chunksize


def test_build_transfer_config_pass_unknown_file_size_no_override_uses_unknown_size_chunksize() -> None:
    """With file_size=None and no override, multipart_chunksize is set to UNKNOWN_SIZE_CHUNKSIZE."""
    config = build_transfer_config(None)
    assert config.multipart_chunksize == UNKNOWN_SIZE_CHUNKSIZE


def test_build_transfer_config_fail_unknown_file_size_override_is_silently_discarded() -> None:
    """When file_size is None, an explicit multipart_chunksize override is discarded in favor of UNKNOWN_SIZE_CHUNKSIZE.

    Documents a real quirk: `base_chunksize = requested_chunksize if file_size is not None else UNKNOWN_SIZE_CHUNKSIZE`
    ignores the caller's override entirely once file_size is unknown, rather than treating it as a
    starting point the way the known-file_size branch does.
    """
    config = build_transfer_config(None, multipart_chunksize=999 * MB)
    assert config.multipart_chunksize == UNKNOWN_SIZE_CHUNKSIZE
    assert config.multipart_chunksize != 999 * MB


def test_build_transfer_config_pass_other_overrides_passed_through() -> None:
    """Non-chunksize/threshold overrides (e.g. max_concurrency, use_threads) are passed straight through to TransferConfig."""
    max_concurrency = 3
    max_download_attempts = 2
    config = build_transfer_config(1000, max_concurrency=3, use_threads=False, num_download_attempts=2)
    assert config.max_concurrency == max_concurrency
    assert config.use_threads is False
    assert config.num_download_attempts == max_download_attempts


def test_build_transfer_config_pass_default_use_threads_and_concurrency_when_not_overridden() -> None:
    """Non-overridden TransferConfig fields fall back to boto3's own built-in defaults."""
    config = build_transfer_config(1000)
    assert config.use_threads is TransferConfig.DEFAULTS["use_threads"]
    assert config.max_concurrency == TransferConfig.DEFAULTS["max_concurrency"]


def test_build_transfer_config_fail_unknown_kwarg_raises_type_error() -> None:
    """An override key that TransferConfig itself doesn't recognize propagates as a TypeError."""
    with pytest.raises(TypeError):
        build_transfer_config(1000, not_a_real_transferconfig_field=True)


def test_build_transfer_config_pass_independent_calls_produce_independent_instances() -> None:
    """Successive calls to build_transfer_config return distinct TransferConfig objects."""
    config_a = build_transfer_config(1000)
    config_b = build_transfer_config(1000)
    assert config_a is not config_b
    assert config_a.multipart_chunksize == config_b.multipart_chunksize


def test_build_transfer_config_pass_multipart_threshold_and_chunksize_can_both_be_overridden_together() -> None:
    """multipart_threshold and multipart_chunksize overrides can be combined without interfering with each other."""
    config = build_transfer_config(1000, multipart_threshold=50 * MB, multipart_chunksize=20 * MB)
    assert config.multipart_threshold == 50 * MB
    assert config.multipart_chunksize == 20 * MB
