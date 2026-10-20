"""Shared fixtures for the XML reader tests."""

import gzip
from collections.abc import Callable
from pathlib import Path

import pytest


@pytest.fixture
def xml_path_factory(tmp_path: Path) -> Callable[..., Path]:
    """Return a factory that writes XML content to a real file, optionally gzip-compressed."""

    def _make(content: str, *, gzip_compress: bool = False, filename: str = "data.xml") -> Path:
        if gzip_compress:
            path = tmp_path / f"{filename}.gz"
            with gzip.open(path, "wb") as f:
                f.write(content.encode("utf-8"))
        else:
            path = tmp_path / filename
            path.write_text(content, encoding="utf-8")
        return path

    return _make
