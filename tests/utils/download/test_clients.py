"""Tests for the FileDownloader and AsyncFileDownloader modules."""

import hashlib
import logging
from pathlib import Path

import httpx
import pytest

from cdm_data_loader_utils.utils.download.async_client import AsyncFileDownloader
from cdm_data_loader_utils.utils.download.core import (
    ChecksumMismatchError,
    DownloadError,
    NonRetryableDownloadError,
)
from cdm_data_loader_utils.utils.download.sync_client import FileDownloader
from tests.utils.download.conftest import DownloaderAdapter

DOWNLOAD_URL = "https://example.com/file.txt"


@pytest.fixture
def sample_content() -> bytes:
    """Generate some extremely interesting sample content."""
    return b"Hello world!"


@pytest.fixture
def sample_checksum(sample_content: bytes) -> str:
    """Generate a checksum for the sample content."""
    return hashlib.sha256(sample_content).hexdigest()


@pytest.mark.parametrize("save_loc", [Path("file.txt"), Path("path") / "to" / "file.txt"])
@pytest.mark.asyncio
async def test_download_to_file_system_success(
    save_loc: Path, tmp_path: Path, sample_content: bytes, sample_checksum: str, downloader_adapter: DownloaderAdapter
) -> None:
    """Simple happy path test of a successful download."""

    def handler(_: httpx.Request) -> httpx.Response:
        return httpx.Response(
            200,
            content=sample_content,
        )

    downloader = downloader_adapter.make_downloader(handler)
    destination = tmp_path / save_loc

    result = await downloader_adapter.download(
        downloader,
        DOWNLOAD_URL,
        destination,
        expected_checksum=sample_checksum,
    )
    assert result == destination
    assert destination.is_file()
    assert destination.read_bytes() == sample_content


@pytest.mark.asyncio
async def test_chunked_streaming(tmp_path: Path, downloader_adapter: DownloaderAdapter) -> None:
    """Ensure that chunked downloads work correctly."""
    content = b"a" * 10_000

    def handler(_: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=content)

    downloader = downloader_adapter.make_downloader(handler, chunk_size=1024)

    destination = tmp_path / "file.txt"
    await downloader_adapter.download(downloader, DOWNLOAD_URL, destination)

    # TODO: add in count of number of times fh.write is called by _download_once?
    assert destination.read_bytes() == content


@pytest.mark.parametrize(
    ("status", "status_text", "error", "msg"),
    [
        (200, None, None, None),
        (302, "Found", DownloadError, "Server"),
        (304, None, None, None),
        (404, "Not Found", NonRetryableDownloadError, "Client"),
        (400, "Bad Request", NonRetryableDownloadError, "Client"),
        (503, "Service Unavailable", DownloadError, "Server"),
    ],
)
@pytest.mark.asyncio
async def test_download_validate_response(  # noqa: PLR0913
    status: int,
    status_text: str | None,
    error: type[Exception] | None,
    msg: str | None,
    tmp_path: Path,
    caplog: pytest.LogCaptureFixture,
    downloader_adapter: DownloaderAdapter,
) -> None:
    """Check that responses of various types emit the correct errors."""

    def handler(_: httpx.Request) -> httpx.Response:
        return httpx.Response(status)

    downloader = downloader_adapter.make_downloader(handler, max_attempts=1)
    destination = tmp_path / "file.txt"

    if error:
        with pytest.raises(error, match=f"{msg} error: {status} {status_text}"):
            await downloader_adapter.download(downloader, DOWNLOAD_URL, destination)

        last_log_msg = caplog.records[-1]
        assert last_log_msg.levelno == logging.ERROR
        assert (
            last_log_msg.message
            == f"{DOWNLOAD_URL}: {msg} error: {status} {status_text}; retry{' not ' if msg == 'Client' else ' '}possible"
        )
        return

    result = await downloader_adapter.download(
        downloader,
        DOWNLOAD_URL,
        destination,
    )
    last_log_msg = caplog.records[-1]
    assert last_log_msg.levelno == logging.INFO

    if status == 304:  # noqa: PLR2004
        assert result is None
        assert last_log_msg.message.startswith(f"{DOWNLOAD_URL}: resource has not been modified")
    else:
        assert result == destination
        assert last_log_msg.message == f"{DOWNLOAD_URL}: download successful"


@pytest.mark.asyncio
async def test_extra_headers(
    tmp_path: Path, caplog: pytest.LogCaptureFixture, downloader_adapter: DownloaderAdapter
) -> None:
    """Test that appropriate headers are added to a request.

    Also tests the 304 response, where a resource is not downloaded as it has not been modified.
    """
    extra_headers = {
        "If-Modified-Since": "Wed, 21 Oct 2015 07:28:00 GMT",
        "X-Sender": "my-user-name",
        "Some-Other-Header": "some value",
    }

    def handler(request: httpx.Request) -> httpx.Response:
        for header, val in extra_headers.items():
            assert request.headers[header] == val
        return httpx.Response(304)

    downloader = downloader_adapter.make_downloader(handler)

    result = await downloader_adapter.download(
        downloader,
        DOWNLOAD_URL,
        tmp_path / "file.txt",
        extra_headers=extra_headers,
    )

    assert result is None
    last_log_msg = caplog.records[-1]
    assert last_log_msg.levelno == logging.INFO
    assert last_log_msg.message.startswith(f"{DOWNLOAD_URL}: resource has not been modified")


@pytest.mark.parametrize("checksum_fn", [None, *hashlib.algorithms_available])
@pytest.mark.parametrize("content", [b"Hello, world!", b"a" * 256])
@pytest.mark.asyncio
async def test_checksum_success(
    tmp_path: Path,
    content: bytes,
    checksum_fn: str | None,
    caplog: pytest.LogCaptureFixture,
    downloader_adapter: DownloaderAdapter,
) -> None:
    """Test that a correct checksum does not throw an error."""
    checksum_fn_used = checksum_fn or "sha256"
    if checksum_fn_used.startswith("shake"):
        pytest.skip("shake algorithms not supported")

    def handler(_: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=content)

    downloader = downloader_adapter.make_downloader(handler, chunk_size=16)
    destination = tmp_path / "file.txt"

    hl = hashlib.new(checksum_fn_used)
    hl.update(content)
    expected_checksum = hl.hexdigest()

    result = await downloader_adapter.download(
        downloader, DOWNLOAD_URL, destination, expected_checksum=expected_checksum, checksum_fn=checksum_fn
    )
    assert result == destination
    assert destination.read_bytes() == content
    assert f"{DOWNLOAD_URL}: {checksum_fn_used} checksum matches" in [m.message for m in caplog.records]


@pytest.mark.parametrize("checksum_fn", ["shake128", " shake256", "shake-n-vac"])
def test_checksum_invalid_sync_downloader(checksum_fn: str, tmp_path: Path) -> None:
    """Ensure that invalid checksums throw an error."""
    downloader = FileDownloader()
    with pytest.raises(ValueError, match=f"Hashing algorithm {checksum_fn} not supported"):
        downloader.download(DOWNLOAD_URL, tmp_path, checksum_fn=checksum_fn)


@pytest.mark.parametrize("checksum_fn", ["shake128", " shake256", "shake-n-vac"])
@pytest.mark.asyncio
async def test_checksum_invalid_async_downloader(checksum_fn: str, tmp_path: Path) -> None:
    """Ensure that invalid checksums throw an error."""
    downloader = AsyncFileDownloader()
    with pytest.raises(ValueError, match=f"Hashing algorithm {checksum_fn} not supported"):
        await downloader.download(DOWNLOAD_URL, tmp_path, checksum_fn=checksum_fn)


@pytest.mark.asyncio
async def test_checksum_mismatch(tmp_path: Path, sample_content: bytes, downloader_adapter: DownloaderAdapter) -> None:
    """Test checksum checking!"""

    def handler(_: httpx.Request) -> httpx.Response:
        return httpx.Response(200, content=sample_content)

    downloader = downloader_adapter.make_downloader(handler)
    destination = tmp_path / "file.txt"

    with pytest.raises(ChecksumMismatchError):
        await downloader_adapter.download(
            downloader,
            DOWNLOAD_URL,
            destination,
            expected_checksum="some_checksum_or_other",
        )
    # N.b. file still exists, even though checksum validation failed.
    assert destination.read_bytes() == sample_content


@pytest.mark.parametrize("max_attempts", range(1, 5))
@pytest.mark.parametrize("response_type", ["error", "timeout"])
@pytest.mark.asyncio
async def test_timeout_and_server_error_retries(  # noqa: PLR0913
    tmp_path: Path,
    sample_content: bytes,
    max_attempts: int,
    response_type: str,
    caplog: pytest.LogCaptureFixture,
    downloader_adapter: DownloaderAdapter,
) -> None:
    """Test retry functionality with different numbers of retries. Server errors are retried."""
    calls = {"n": 0}
    # the attempt where the download will succeed
    successful_attempt = 3

    def handler_503(_: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        if calls["n"] < successful_attempt:
            return httpx.Response(503)
        return httpx.Response(200, content=sample_content)

    def handler_timeout(_: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        if calls["n"] < successful_attempt:
            msg = "Connection has timed out"
            raise httpx.ConnectTimeout(msg)
        return httpx.Response(200, content=sample_content)

    downloader = downloader_adapter.make_downloader(
        handler_503 if response_type == "error" else handler_timeout, max_attempts=max_attempts
    )
    destination = tmp_path / "file.txt"

    if max_attempts < successful_attempt:
        expected_match = (
            "Server error: 503 Service Unavailable" if response_type == "error" else "Connection has timed out"
        )
        with pytest.raises(DownloadError, match=expected_match):
            await downloader_adapter.download(
                downloader,
                DOWNLOAD_URL,
                tmp_path / "file.txt",
            )
        assert calls["n"] == max_attempts
    else:
        # successful download
        result = await downloader_adapter.download(downloader, DOWNLOAD_URL, destination)
        assert result == destination
        assert calls["n"] == successful_attempt

    # check the logs: we should have messages for each error/timeout until the successful download
    retry_possible_msgs = [r for r in caplog.records if r.message.endswith("retry possible")]
    if max_attempts < successful_attempt:
        assert len(retry_possible_msgs) == max_attempts
    else:
        assert len(retry_possible_msgs) == successful_attempt - 1


@pytest.mark.parametrize("max_attempts", range(1, 5))
@pytest.mark.asyncio
async def test_client_error_retries(
    tmp_path: Path, sample_content: bytes, max_attempts: int, downloader_adapter: DownloaderAdapter
) -> None:
    """Test retry functionality with different numbers of retries. Client errors are not retried."""
    calls = {"n": 0}
    # the attempt where the download will succeed
    successful_attempt = 3

    def handler(_: httpx.Request) -> httpx.Response:
        calls["n"] += 1
        if calls["n"] < successful_attempt:
            return httpx.Response(404)
        return httpx.Response(200, content=sample_content)

    downloader = downloader_adapter.make_downloader(handler, max_attempts=max_attempts)

    with pytest.raises(NonRetryableDownloadError, match="Client error: 404 Not Found"):
        await downloader_adapter.download(
            downloader,
            DOWNLOAD_URL,
            tmp_path / "file.txt",
        )
    # there will only ever be one retry
    assert calls["n"] == 1
