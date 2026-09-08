"""Tests for utils.file_transfer.ftp_client module — mock ftplib for keepalive, retry, thread-local."""

import socket
import time
from collections.abc import Callable
from ftplib import FTP, error_temp
from pathlib import Path, PurePosixPath
from unittest.mock import MagicMock, patch

import pytest

from cdm_data_loaders.utils.file_transfer.ftp_client import (
    ThreadLocalFTP,
    _set_keepalive,
    connect_ftp,
    ftp_download_file,
    ftp_list_dir,
    ftp_noop_keepalive,
    ftp_retrieve_text,
)

_IDLE_SECONDS = 30
_KEEPIDLE_VALUE = 60
_KEEPALIVE_INTERVAL = 25
_EXPECTED_RETRY_COUNT = 2
_FTP_TIMEOUT = 30
_ERR_421 = "421 timeout"


class TestSetKeepalive:
    """Test TCP keepalive socket options."""

    def test_sets_so_keepalive(self) -> None:
        """Verify SO_KEEPALIVE is set on the socket."""
        mock_ftp = MagicMock(spec=FTP)
        mock_sock = MagicMock()
        mock_ftp.sock = mock_sock
        _set_keepalive(mock_ftp)
        mock_sock.setsockopt.assert_any_call(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)

    def test_sets_tcp_keepidle(self) -> None:
        """Verify TCP_KEEPIDLE is set when available."""
        mock_ftp = MagicMock(spec=FTP)
        mock_sock = MagicMock()
        mock_ftp.sock = mock_sock
        _set_keepalive(mock_ftp, idle=_KEEPIDLE_VALUE)
        if hasattr(socket, "TCP_KEEPIDLE"):
            mock_sock.setsockopt.assert_any_call(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, _KEEPIDLE_VALUE)


class TestConnectFtp:
    """Test connect_ftp creates and configures an FTP connection."""

    @patch("cdm_data_loaders.utils.file_transfer.ftp_client.FTP")
    def test_connect_and_login(self, mock_ftp_cls: MagicMock) -> None:
        """Verify FTP object is created, login called, and returned."""
        mock_ftp = MagicMock()
        mock_ftp.sock = MagicMock()
        mock_ftp_cls.return_value = mock_ftp
        result = connect_ftp("ftp.example.com", timeout=_FTP_TIMEOUT)
        mock_ftp_cls.assert_called_once_with("ftp.example.com", timeout=_FTP_TIMEOUT)
        mock_ftp.login.assert_called_once()
        assert result is mock_ftp


class TestFtpNoopKeepalive:
    """Test NOOP keepalive logic."""

    def test_sends_noop_when_idle(self) -> None:
        """Verify NOOP is sent when idle exceeds interval."""
        mock_ftp = MagicMock(spec=FTP)
        old_time = time.monotonic() - _IDLE_SECONDS
        new_time = ftp_noop_keepalive(mock_ftp, old_time, interval=_KEEPALIVE_INTERVAL)
        mock_ftp.sendcmd.assert_called_once_with("NOOP")
        assert new_time > old_time

    def test_no_noop_when_recent(self) -> None:
        """Verify no NOOP is sent when activity is recent."""
        mock_ftp = MagicMock(spec=FTP)
        recent = time.monotonic()
        result = ftp_noop_keepalive(mock_ftp, recent, interval=_KEEPALIVE_INTERVAL)
        mock_ftp.sendcmd.assert_not_called()
        assert result == recent


class TestFtpListDir:
    """Test ftp_list_dir with retry."""

    def test_returns_file_list(self) -> None:
        """Verify file listing is returned correctly."""
        mock_ftp = MagicMock(spec=FTP)

        def fake_retrlines(_cmd: str, callback: Callable[[str], None]) -> None:
            for name in ["file1.txt", "file2.gz"]:
                callback(name)

        mock_ftp.retrlines.side_effect = fake_retrlines
        result = ftp_list_dir(mock_ftp, PurePosixPath("/") / "some" / "path")
        assert result == [PurePosixPath("file1.txt"), PurePosixPath("file2.gz")]
        mock_ftp.cwd.assert_called_once_with("/some/path")

    def test_retries_on_error_temp(self) -> None:
        """Verify retry logic on FTP temporary errors."""
        mock_ftp = MagicMock(spec=FTP)
        call_count = 0

        def fake_retrlines(_cmd: str, callback: Callable[[str], None]) -> None:
            nonlocal call_count
            call_count += 1
            if call_count < _EXPECTED_RETRY_COUNT:
                raise error_temp(_ERR_421)  # noqa: S321
            callback("file.txt")

        mock_ftp.retrlines.side_effect = fake_retrlines
        result = ftp_list_dir(mock_ftp, PurePosixPath("/path"), retries=3)
        assert result == [PurePosixPath("file.txt")]
        assert call_count == _EXPECTED_RETRY_COUNT

    def test_raises_after_exhausted_retries(self) -> None:
        """Verify error is raised after all retries are exhausted."""
        mock_ftp = MagicMock(spec=FTP)
        mock_ftp.retrlines.side_effect = error_temp(_ERR_421)  # noqa: S321
        with pytest.raises(error_temp):
            ftp_list_dir(mock_ftp, PurePosixPath("/") / "path", retries=_EXPECTED_RETRY_COUNT)


class TestFtpDownloadFile:
    """Test ftp_download_file with retry."""

    def test_downloads_file(self, tmp_path: Path) -> None:
        """Verify file is downloaded and written to disk."""
        mock_ftp = MagicMock(spec=FTP)

        def fake_retrbinary(_cmd: str, callback: Callable[[bytes], None]) -> None:
            callback(b"file data")

        mock_ftp.retrbinary.side_effect = fake_retrbinary
        local = tmp_path / "out.bin"
        ftp_download_file(mock_ftp, PurePosixPath("remote.bin"), local)
        assert local.read_bytes() == b"file data"

    def test_retries_on_error_temp(self, tmp_path: Path) -> None:
        """Verify download retries on FTP temporary errors."""
        mock_ftp = MagicMock(spec=FTP)
        call_count = 0

        def fake_retrbinary(_cmd: str, callback: Callable[[bytes], None]) -> None:
            nonlocal call_count
            call_count += 1
            if call_count < _EXPECTED_RETRY_COUNT:
                msg = "421"
                raise error_temp(msg)  # noqa: S321
            callback(b"ok")

        mock_ftp.retrbinary.side_effect = fake_retrbinary
        local = tmp_path / "out.bin"
        ftp_download_file(mock_ftp, PurePosixPath("remote.bin"), local, retries=3)
        assert call_count == _EXPECTED_RETRY_COUNT


class TestFtpRetrieveText:
    """Test ftp_retrieve_text."""

    def test_returns_content(self) -> None:
        """Verify text content is retrieved and joined with newlines."""
        mock_ftp = MagicMock(spec=FTP)

        def fake_retrlines(_cmd: str, callback: Callable[[str], None]) -> None:
            for line in ["line1", "line2"]:
                callback(line)

        mock_ftp.retrlines.side_effect = fake_retrlines
        result = ftp_retrieve_text(mock_ftp, PurePosixPath("remote.txt"))
        assert result == "line1\nline2"


class TestThreadLocalFTP:
    """Test thread-local FTP connection management."""

    @patch("cdm_data_loaders.utils.file_transfer.ftp_client.connect_ftp")
    def test_get_returns_same_connection(self, mock_connect: MagicMock) -> None:
        """Verify get() returns the same FTP connection on repeated calls."""
        mock_ftp = MagicMock()
        mock_connect.return_value = mock_ftp
        pool = ThreadLocalFTP("ftp.example.com")
        ftp1 = pool.get()
        ftp2 = pool.get()
        assert ftp1 is ftp2
        mock_connect.assert_called_once()

    @patch("cdm_data_loaders.utils.file_transfer.ftp_client.connect_ftp")
    def test_close_all(self, mock_connect: MagicMock) -> None:
        """Verify close_all() quits the FTP connection."""
        mock_ftp = MagicMock()
        mock_connect.return_value = mock_ftp
        pool = ThreadLocalFTP("ftp.example.com")
        pool.get()
        pool.close_all()
        mock_ftp.quit.assert_called_once()
