"""File downloader with tenacity for retry support."""

import hashlib
import logging
from collections.abc import AsyncIterable, Iterable
from pathlib import Path
from typing import Any

import httpx

from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger, log_and_die

logger: logging.Logger = get_cdm_logger()


class DownloadError(Exception):
    """Retryable download error."""


class NonRetryableDownloadError(Exception):
    """Retryable download error."""


class ChecksumMismatchError(NonRetryableDownloadError):
    """Checksum validation failed."""


class DownloadCore:
    """
    HTTP file downloader with retries and checksum testing.
    """

    RETRYABLE_EXCEPTIONS = (
        httpx.TimeoutException,
        httpx.TransportError,
        DownloadError,
    )

    @staticmethod
    def validate_args(
        destination: str | Path,
        checksum_fn: str | None = None,
        extra_headers: dict[str, Any] | None = None,
    ) -> tuple[Path, str, dict[str, Any]]:
        """Validate download arguments, apply defaults, and convert them to the appropriate format.

        If no checksum function is supplied, it will be set to "sha256".

        :param destination: location to save file to
        :type destination: str | Path
        :param checksum_fn: function for calculating file checksum
        :type checksum_fn: str | None, optional
        :param extra_headers: extra headers to add to the httpx request, defaults to None
        :type extra_headers: dict[str, Any] | None, optional
        :return: arguments spruced up and in shape!
        :rtype: tuple[Path, str, dict[str, Any]]
        """
        if not checksum_fn:
            checksum_fn = "sha256"

        # shake hashing algorithms require a length for the hexdigest => too much faff, error out
        if checksum_fn not in hashlib.algorithms_available or checksum_fn.startswith("shake"):
            log_and_die(f"Hashing algorithm {checksum_fn} not supported.", ValueError)

        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)

        if not extra_headers:
            extra_headers = {}

        return (destination, checksum_fn, extra_headers)

    @staticmethod
    def validate_response(response: httpx.Response) -> bool:
        """Check the response from the server and act accordingly."""
        status = response.status_code

        if status == 304:  # noqa: PLR2004
            return False

        # successful response
        if 200 <= status < 300:  # noqa: PLR2004
            return True

        # client error
        if 400 <= status < 500:  # noqa: PLR2004
            msg = f"Client error: {status} {response.reason_phrase}"
            raise NonRetryableDownloadError(msg)

        # server errors
        msg = f"Server error: {status} {response.reason_phrase}"
        raise DownloadError(msg)

    @staticmethod
    def check_hash(
        url: str,
        expected_checksum: str,
        checksum_fn: str,
        hash_obj,  #: hashlib._Hash,
    ) -> None:
        """Compare the expected to actual checksum.

        :param url: download URL
        :type url: str
        :param expected_checksum: expected checksum
        :type expected_checksum: str
        :param checksum_fn: checksum algorithm
        :type checksum_fn: str
        :param hash_obj: hash object for calculating checksum
        :type hash_obj: bool
        :raises ChecksumMismatchError: if the checksum is incorrect
        """
        if hash_obj:
            # check the checksum matches expectations
            actual = hash_obj.hexdigest()
            if actual.lower() != expected_checksum.lower():
                msg = f"{url}: Checksum mismatch: expected={expected_checksum}, actual={actual}"
                raise ChecksumMismatchError(msg)
            logger.info("%s: %s checksum matches", url, checksum_fn)

    @staticmethod
    def write_and_hash(
        url: str,
        destination: Path,
        chunks: Iterable[bytes],
        expected_checksum: str | None,
        checksum_fn: str,
    ) -> None:
        """Write downloaded data to disk and accrue data for calculating the checksum of the download.

        :param url: URL to download from
        :type url: str
        :param destination: where to save the file to
        :type destination: Path
        :param chunks: iterable chunks of downloaded data
        :type chunks: Iterable[bytes]
        :param expected_checksum: expected checksum (if known), defaults to None
        :type expected_checksum: str | None, optional
        :param checksum_fn: function used for calculating the checksum, defaults to "sha256"
        :type checksum_fn: str, optional
        :raises ChecksumMismatchError: if the checksum of the downloaded data does not match the expected checksum.
        """
        hash_obj = hashlib.new(checksum_fn) if expected_checksum else None

        with destination.open("wb") as fh:
            for chunk in chunks:
                fh.write(chunk)
                if hash_obj:
                    hash_obj.update(chunk)

        # include expected_checksum to stop the typechecker complaining
        if hash_obj and expected_checksum:
            DownloadCore.check_hash(url, expected_checksum, checksum_fn, hash_obj)

    @staticmethod
    async def write_and_hash_async(
        url: str,
        destination: Path,
        chunks: AsyncIterable[bytes],
        expected_checksum: str | None,
        checksum_fn: str,
    ) -> None:
        """Write downloaded data to disk and accrue data for calculating the checksum of the download.

        :param url: URL to download from
        :type url: str
        :param destination: where to save the file to
        :type destination: Path
        :param chunks: iterable chunks of downloaded data
        :type chunks: AsyncIterable[bytes]
        :param expected_checksum: expected checksum (if known), defaults to None
        :type expected_checksum: str | None, optional
        :param checksum_fn: function used for calculating the checksum, defaults to "sha256"
        :type checksum_fn: str, optional
        :raises ChecksumMismatchError: if the checksum of the downloaded data does not match the expected checksum.
        """
        hash_obj = hashlib.new(checksum_fn) if expected_checksum else None

        # NOTE: file I/O is sync, but chunk iteration must be async
        with destination.open("wb") as fh:
            async for chunk in chunks:
                fh.write(chunk)
                if hash_obj:
                    hash_obj.update(chunk)

        # include expected_checksum to stop the typechecker complaining
        if hash_obj and expected_checksum:
            DownloadCore.check_hash(url, expected_checksum, checksum_fn, hash_obj)
