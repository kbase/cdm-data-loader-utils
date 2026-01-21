"""Synchronous file download client.

Usage:

# initialise a download client that downloads files in chunks of 1024 bytes and retries the downlooad
# three times in case of connection/timeout errors
sync_client = FileDownloader(max_attempts=3, chunk_size=1024)

# output_path will be the path to the saved file
output_path = sync_client.download(
    "https://example.com/file.txt",
    Path("/path") / "to" / "save" / "file.txt"),
    expected_checksum="some_checksum",
    checksum_fn="sha512",
    extra_headers={
        "If-None-Match": "abc",
        "If-Modified-Since": "Wed, 21 Oct 2015 07:28:00 GMT"
    })

"""

import logging
from pathlib import Path
from typing import Any

import httpx
from tenacity import (
    before_sleep_log,
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from cdm_data_loader_utils.utils.download.core import (
    DownloadCore,
    DownloadError,
    NonRetryableDownloadError,
)
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger

logger: logging.Logger = get_cdm_logger()


def get_httpx_client() -> httpx.Client:
    """Get a basic client for executing http requests."""
    return httpx.Client(
        timeout=httpx.Timeout(30.0),
        limits=httpx.Limits(
            max_connections=20,
            max_keepalive_connections=10,
        ),
        follow_redirects=True,
    )


class FileDownloader:
    """
    Synchronous downloader interface.
    """

    RETRYABLE_EXCEPTIONS = (
        httpx.TimeoutException,
        httpx.TransportError,
        DownloadError,
    )

    def __init__(
        self,
        client: httpx.Client | None = None,
        max_attempts: int = 5,
        min_backoff: int = 1,
        max_backoff: int = 30,
        chunk_size: int = 8192,
    ) -> None:
        """Initialise a synchronous download client.

        :param client: an httpx.Client object, defaults to None
        :type client: httpx.Client | None, optional
        :param max_attempts: how many times to retry the request defaults to 5
        :type max_attempts: int, optional
        :param min_backoff: minimum backoff for retries, defaults to 1
        :type min_backoff: int, optional
        :param max_backoff: maximum backoff for retries, defaults to 30
        :type max_backoff: int, optional
        :param chunk_size: chunk size for the downloader, defaults to 8192
        :type chunk_size: int, optional
        """
        self.client = client or get_httpx_client()
        self.chunk_size = chunk_size

        self._retry = retry(
            retry=retry_if_exception_type(self.RETRYABLE_EXCEPTIONS),
            stop=stop_after_attempt(max_attempts),
            wait=wait_exponential(
                min=min_backoff,
                max=max_backoff,
            ),
            reraise=True,
            before_sleep=before_sleep_log(logger, logging.WARNING),
        )

    def download(
        self,
        url: str,
        destination: str | Path,
        expected_checksum: str | None = None,
        checksum_fn: str | None = None,
        extra_headers: dict[str, Any] | None = None,
    ) -> Path | None:
        """Download a file from ``url`` to ``destination`` on disk.

        Returns the file path if successful; None otherwise.
        304 / Not Modified responses return None.

        :param url: URL to download from
        :type url: str
        :param destination: where to save the file to
        :type destination: str | Path
        :param expected_checksum: expected checksum (if known), defaults to None
        :type expected_checksum: str | None, optional
        :param checksum_fn: function used for calculating the checksum. Set to "sha256" if not supplied.
        :type checksum_fn: str, optional
        :param extra_headers: allow extra headers to be passed, defaults to None
        :type last_modified: dict[str, Any | None, optional
        :return: either the path to the downloaded file or None if no file was downloaded
        :rtype: Path | None
        """
        destination, checksum_fn, extra_headers = DownloadCore.validate_args(destination, checksum_fn, extra_headers)

        @self._retry
        def _once() -> Path | None:
            return self._download_once(url, destination, expected_checksum, checksum_fn, extra_headers)

        return _once()

    def _download_once(
        self,
        url: str,
        destination: Path,
        expected_checksum: str | None,
        checksum_fn: str,
        extra_headers: dict[str, Any],
    ) -> Path | None:
        """Core synchronous download function.

        :param url: URL to download from
        :type url: str
        :param destination: where to save the file to
        :type destination: Path
        :param expected_checksum: expected checksum (if known), defaults to None
        :type expected_checksum: str | None, optional
        :param checksum_fn: function used for calculating the checksum, defaults to "sha256"
        :type checksum_fn: str, optional
        :param extra_headers: any extra headers to pass to the request, defaults to None
        :type extra_headers: dict[str, Any] | None, optional
        :raises DownloadError: for unrecoverable download errors
        :raises ChecksumMismatchError: if the checksum does not match expected_checksum
        :return: either the path to the downloaded file or None if no file was downloaded
        :rtype: Path | None
        """
        try:
            with self.client.stream("GET", url, headers=extra_headers) as response:
                continue_dl = DownloadCore.validate_response(response)
                if not continue_dl:
                    logger.info("%s: resource has not been modified", url, extra={"url": url})
                    return None

                chunks = response.iter_bytes(chunk_size=self.chunk_size)

                DownloadCore.write_and_hash(
                    url,
                    destination,
                    chunks,
                    expected_checksum=expected_checksum,
                    checksum_fn=checksum_fn,
                )

        except NonRetryableDownloadError as exc:
            logger.exception("%s: %s; retry not possible", url, exc.args[0], extra={"url": url})
            raise

        except Exception as exc:
            logger.exception("%s: %s; retry possible", url, exc.args[0], extra={"url": url})
            raise DownloadError(str(exc)) from exc

        logger.info("%s: download successful", url, extra={"path": str(destination)})

        return destination
