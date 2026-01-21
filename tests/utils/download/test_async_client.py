"""Tests for the AsyncFileDownloader client that have no counterpart in the sync version."""

import asyncio
from pathlib import Path
from typing import NoReturn, Self

import httpx
import pytest

from cdm_data_loader_utils.utils.download.async_client import AsyncFileDownloader

DOWNLOAD_URL = "https://example.com/file.txt"


@pytest.mark.asyncio
async def test_max_concurrency_is_enforced(tmp_path: Path) -> None:
    """
    No more than max_concurrency downloads may run at once.
    """
    max_concurrency = 2
    active = 0
    peak = 0
    lock = asyncio.Lock()
    block = asyncio.Event()

    async def handler(_: httpx.Request) -> httpx.Response:
        nonlocal active, peak
        async with lock:
            active += 1
            peak = max(peak, active)

        # Block all requests here until released
        await block.wait()
        async with lock:
            active -= 1

        return httpx.Response(200, content=b"ok")

    downloader = AsyncFileDownloader(
        httpx.AsyncClient(transport=httpx.MockTransport(handler)),
        max_concurrency=max_concurrency,
    )

    async def start_download(i: int) -> None:
        dest = tmp_path / f"file_{i}.txt"
        await downloader.download(
            DOWNLOAD_URL,
            dest,
        )

    tasks = [asyncio.create_task(start_download(i)) for i in range(5)]

    # Give tasks time to start and block
    await asyncio.sleep(0.1)
    # Release blocked requests
    block.set()
    await asyncio.gather(*tasks)
    # peak concurrency should be limited by the max_concurrency value
    assert peak == max_concurrency


@pytest.mark.asyncio
async def test_no_concurrency_limit_when_none(tmp_path: Path) -> None:
    """
    When max_concurrency is None, all downloads may run concurrently.
    """
    active = 0
    peak = 0
    lock = asyncio.Lock()
    block = asyncio.Event()

    async def handler(_: httpx.Request) -> httpx.Response:
        nonlocal active, peak
        async with lock:
            active += 1
            peak = max(peak, active)

        await block.wait()
        async with lock:
            active -= 1

        return httpx.Response(200, content=b"ok")

    client = httpx.AsyncClient(transport=httpx.MockTransport(handler))

    downloader = AsyncFileDownloader(
        client,
        max_concurrency=None,
    )

    async def start_download(i: int) -> None:
        dest = tmp_path / f"file_{i}.txt"
        await downloader.download(
            DOWNLOAD_URL,
            dest,
        )

    n_tasks = 5
    tasks = [asyncio.create_task(start_download(i)) for i in range(n_tasks)]

    await asyncio.sleep(0.1)
    block.set()
    await asyncio.gather(*tasks)

    # All tasks should have entered concurrently
    assert peak == n_tasks


class EmptyAsyncRetry:
    """Async iterator that yields no retry attempts."""

    def __aiter__(self) -> Self:
        """Iterator."""
        return self

    async def __anext__(self) -> NoReturn:
        """Next."""
        raise StopAsyncIteration


@pytest.mark.asyncio
async def test_async_client_unreachable_code(monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that the unreachable code is reached if the retry iterator yields no attempts.

    Defensive test, should never occur IRL.
    """
    downloader = AsyncFileDownloader(client=None)

    # Replace the retry mechanism with an empty async iterator
    monkeypatch.setattr(downloader, "_retry", EmptyAsyncRetry())

    with pytest.raises(RuntimeError, match="Iterator exhausted, unreachable code reached!"):
        await downloader._download_with_retry(DOWNLOAD_URL, "file.txt")  # noqa: SLF001
