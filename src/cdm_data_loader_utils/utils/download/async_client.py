"""Asynchronous file download client.

Usage:

# initialise a download client that downloads files in chunks of 1024 bytes and retries the downlooad
# three times in case of connection/timeout errors
async_client = AsyncFileDownloader(max_attempts=3, chunk_size=1024)

# output_path will be the path to the saved file
output_path = await async_client.download(
    "https://example.com/file.txt",
    Path("/path") / "to" / "save" / "file.txt"),
    expected_checksum="some_checksum",
    checksum_fn="sha512",
    extra_headers={
        "If-None-Match": "abc",
        "If-Modified-Since": "Wed, 21 Oct 2015 07:28:00 GMT"
    })

"""

import asyncio
import logging
from pathlib import Path
from typing import Any

import httpx
from tenacity import (
    AsyncRetrying,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger
from cdm_data_loader_utils.utils.download.core import (
    DownloadCore,
    DownloadError,
    NonRetryableDownloadError,
)

logger: logging.Logger = get_cdm_logger()


def get_async_httpx_client() -> httpx.AsyncClient:
    """Get a basic client for executing http requests."""
    return httpx.AsyncClient(
        timeout=httpx.Timeout(30.0),
        limits=httpx.Limits(
            max_connections=20,
            max_keepalive_connections=10,
        ),
        follow_redirects=True,
    )


class AsyncFileDownloader:
    """
    Asynchronous downloader interface.
    """

    RETRYABLE_EXCEPTIONS = (
        httpx.TimeoutException,
        httpx.TransportError,
        DownloadError,
    )

    def __init__(  # noqa: PLR0913
        self,
        client: httpx.AsyncClient | None = None,
        max_attempts: int = 5,
        min_backoff: int = 1,
        max_backoff: int = 30,
        chunk_size: int = 8192,
        max_concurrency: int | None = None,
    ) -> None:
        """Initialise an async download client.

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
        :param max_concurrency: maximum concurrency for the downloader, defaults to None
        :type max_concurrency: int, optional
        """
        self.client = client or get_async_httpx_client()
        self.chunk_size = chunk_size
        self.semaphore = asyncio.Semaphore(max_concurrency) if max_concurrency else None

        self._retry = AsyncRetrying(
            retry=retry_if_exception_type(self.RETRYABLE_EXCEPTIONS),
            stop=stop_after_attempt(max_attempts),
            wait=wait_exponential(
                min=min_backoff,
                max=max_backoff,
            ),
            reraise=True,
        )

    async def download(
        self,
        url: str,
        destination: str | Path,
        expected_checksum: str | None = None,
        checksum_fn: str | None = "sha256",
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
        :param checksum_fn: function used for calculating the checksum, defaults to "sha256"
        :type checksum_fn: str, optional
        :param extra_headers: allow extra headers to be passed, defaults to None
        :type last_modified: dict[str, Any | None, optional
        :return: either the path to the downloaded file or None if no file was downloaded
        :rtype: Path | None
        """
        destination, checksum_fn, extra_headers = DownloadCore.validate_args(destination, checksum_fn, extra_headers)

        if self.semaphore:
            async with self.semaphore:
                return await self._download_with_retry(
                    url,
                    destination,
                    expected_checksum,
                    checksum_fn,
                    extra_headers,
                )

        return await self._download_with_retry(
            url,
            destination,
            expected_checksum,
            checksum_fn,
            extra_headers,
        )

    async def _download_with_retry(self, *args) -> Path | None:
        async for attempt in self._retry:
            with attempt:
                return await self._download_once(*args)
        # this should never be raised under normal circumstances
        # tenacity should never exhaust silently
        msg = "Iterator exhausted, unreachable code reached!"
        raise RuntimeError(msg)

    async def _download_once(
        self,
        url: str,
        destination: Path,
        expected_checksum: str | None,
        checksum_fn: str,
        extra_headers: dict[str, Any],
    ) -> Path | None:
        """Core download function.

        :param url: URL to download from
        :type url: str
        :param destination: where to save the file to
        :type destination: Path
        :param expected_checksum: expected checksum (if known), defaults to None
        :type expected_checksum: str | None, optional
        :param checksum_fn: function used for calculating the checksum
        :type checksum_fn: str
        :param extra_headers: any extra headers to pass to the request
        :type extra_headers: dict[str, Any]
        :raises DownloadError: for unrecoverable download errors
        :raises ChecksumMismatchError: if the checksum does not match expected_checksum
        :return: either the path to the downloaded file or None if no file was downloaded
        :rtype: Path | None
        """
        try:
            async with self.client.stream("GET", url, headers=extra_headers) as response:
                # check the response - should the download continue?
                if not DownloadCore.validate_response(response):
                    logger.info("%s: resource has not been modified", url, extra={"url": url})
                    return None

                chunks = response.aiter_bytes(chunk_size=self.chunk_size)

                await DownloadCore.write_and_hash_async(
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
