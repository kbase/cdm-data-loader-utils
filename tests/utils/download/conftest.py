"""Fixtures for testing downloads."""

from collections.abc import Awaitable, Callable
from dataclasses import dataclass
from pathlib import Path

import httpx
import pytest

from cdm_data_loader_utils.utils.download.async_client import AsyncFileDownloader
from cdm_data_loader_utils.utils.download.sync_client import FileDownloader


@dataclass
class DownloaderAdapter:
    """Class providing a uniform API to the sync and async file downloader."""

    name: str
    make_downloader: Callable
    download: Callable[..., Awaitable[Path | None]]


@pytest.fixture(params=["sync", "async"])
def downloader_adapter(request: pytest.FixtureRequest) -> DownloaderAdapter:
    """Fixture to generate sync or async file download clients."""

    def make_sync(handler: Callable, **kwargs) -> FileDownloader:  # noqa: ANN003
        client = httpx.Client(transport=httpx.MockTransport(handler))
        return FileDownloader(client, **kwargs)

    def make_async(handler: Callable, **kwargs) -> AsyncFileDownloader:  # noqa: ANN003
        client = httpx.AsyncClient(transport=httpx.MockTransport(handler))
        return AsyncFileDownloader(client, **kwargs)

    async def download_sync(downloader: FileDownloader, *args, **kwargs) -> Path | None:  # noqa: ANN002, ANN003
        # Call sync method inside async test
        return downloader.download(*args, **kwargs)

    async def download_async(downloader: AsyncFileDownloader, *args, **kwargs) -> Path | None:  # noqa: ANN002, ANN003
        return await downloader.download(*args, **kwargs)

    if request.param == "sync":
        return DownloaderAdapter(
            name="sync",
            make_downloader=make_sync,
            download=download_sync,
        )

    return DownloaderAdapter(
        name="async",
        make_downloader=make_async,
        download=download_async,
    )
