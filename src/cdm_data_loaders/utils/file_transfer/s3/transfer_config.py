"""Compute a safe boto3 TransferConfig for a given (or unknown) file size.

Kept free of any S3/HTTP I/O so the chunksize arithmetic can be unit tested
directly against plain integers.
"""

import math
from logging import Logger, getLogger
from typing import Any, Final

from boto3.s3.transfer import TransferConfig

logger: Logger = getLogger(__name__)

MB: Final[int] = 1024**2

# S3 hard limits (see AWS multipart upload documentation / s3transfer constants)
S3_MAX_PARTS: Final[int] = 10_000
S3_MIN_PART_SIZE: Final[int] = 5 * MB

# Matches s3transfer's own defaults, used as the starting point before scaling up.
DEFAULT_MULTIPART_THRESHOLD: Final[int] = 8 * MB
DEFAULT_MULTIPART_CHUNKSIZE: Final[int] = 8 * MB

# Chunk size used when a file's size cannot be determined ahead of the
# transfer (e.g. chunked transfer-encoding with no Content-Length). Large
# enough to safely cover files up to ~625 GB without hitting the part limit,
# at the cost of slightly less parallelism for small unknown-size transfers.
UNKNOWN_SIZE_CHUNKSIZE: Final[int] = 64 * MB


def compute_multipart_chunksize(file_size: int | None, requested_chunksize: int | None = None) -> int:
    """Compute a multipart chunksize that will not exceed S3's 10,000-part limit.

    If `file_size` would require more than `S3_MAX_PARTS` parts at
    `requested_chunksize`, the chunksize is scaled up just enough to bring
    the part count back within S3's limit; a warning is logged in that case,
    since it usually indicates an unusually large file. The result is always
    at least `S3_MIN_PART_SIZE`, since S3 rejects smaller multipart parts.

    :param file_size: size of the file to be uploaded, in bytes, or None if unknown
    :type file_size: int | None
    :param requested_chunksize: chunksize to use as a starting point, in
        bytes, defaults to `DEFAULT_MULTIPART_CHUNKSIZE` if not given
    :type requested_chunksize: int | None, optional
    :return: a chunksize, in bytes, safe to use as `TransferConfig.multipart_chunksize`
        for a file of `file_size` bytes
    :rtype: int
    """
    chunksize = requested_chunksize or DEFAULT_MULTIPART_CHUNKSIZE

    if file_size is not None and file_size > 0:
        required_chunksize = math.ceil(file_size / S3_MAX_PARTS)
        if required_chunksize > chunksize:
            logger.warning(
                "File size (%d bytes) would require more than %d parts at chunksize=%d bytes; "
                "increasing multipart_chunksize to %d bytes to stay within S3's part limit",
                file_size,
                S3_MAX_PARTS,
                chunksize,
                required_chunksize,
            )
            chunksize = required_chunksize

    return max(chunksize, S3_MIN_PART_SIZE)


def build_transfer_config(file_size: int | None, **overrides: Any) -> TransferConfig:
    """Build a `boto3.s3.transfer.TransferConfig` sized appropriately for `file_size`.

    Always computes a safe `multipart_chunksize` via `compute_multipart_chunksize`,
    using any explicitly requested `multipart_chunksize` in `overrides` as the
    starting point (and scaling it up if necessary), rather than silently
    ignoring the caller's preference.

    :param file_size: size of the file to be uploaded, in bytes, or None if
        unknown ahead of time (e.g. no Content-Length header available)
    :type file_size: int | None
    :param overrides: any other `TransferConfig` keyword arguments to pass
        through (e.g. `multipart_threshold`, `max_concurrency`); if
        `multipart_chunksize` is included, it is treated as a minimum rather
        than an absolute value
    :type overrides: Any
    :return: a TransferConfig with a chunksize guaranteed not to exceed
        S3's 10,000-part limit for a file of `file_size` bytes
    :rtype: TransferConfig
    """
    requested_chunksize = overrides.pop("multipart_chunksize", None)
    base_chunksize = requested_chunksize if file_size is not None else UNKNOWN_SIZE_CHUNKSIZE
    chunksize = compute_multipart_chunksize(file_size, base_chunksize)

    overrides.setdefault("multipart_threshold", DEFAULT_MULTIPART_THRESHOLD)

    return TransferConfig(multipart_chunksize=chunksize, **overrides)
