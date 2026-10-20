"""General-purpose FTP client utilities.

Provides resilient FTP connections with TCP keepalive, NOOP pings, retry
on transient errors, and thread-local connection management for parallel
downloads.  Protocol-agnostic — callers supply the FTP hostname.
"""

import contextlib
import logging
import socket
import threading
import time
from ftplib import FTP, error_temp
from logging import Logger, getLogger
from pathlib import Path, PurePosixPath

from tenacity import (
    before_sleep_log,
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_fixed,
)

logger: Logger = getLogger(__name__)

DEFAULT_TIMEOUT = 60


def connect_ftp(host: str, timeout: int = DEFAULT_TIMEOUT) -> FTP:
    """Connect and log in to an FTP server with TCP keepalive enabled.

    :param host: FTP hostname
    :param timeout: connection timeout in seconds
    :return: logged-in FTP connection
    """
    ftp = FTP(host, timeout=timeout)  # noqa: S321
    ftp.login()
    _set_keepalive(ftp)
    return ftp


def _set_keepalive(ftp: FTP, idle: int = 30, interval: int = 10, count: int = 3) -> None:
    """Enable TCP keepalive on the FTP control socket.

    Prevents idle-timeout disconnects (e.g. '421 No transfer timeout') when
    the control connection sits idle during data transfers or checksum
    verification.
    """
    sock = ftp.sock
    if sock is None:
        return
    sock.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE, 1)
    if hasattr(socket, "TCP_KEEPIDLE"):
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPIDLE, idle)
    if hasattr(socket, "TCP_KEEPINTVL"):
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPINTVL, interval)
    if hasattr(socket, "TCP_KEEPCNT"):
        sock.setsockopt(socket.IPPROTO_TCP, socket.TCP_KEEPCNT, count)


def ftp_noop_keepalive(ftp: FTP, last_activity: float, interval: int = 25) -> float:
    """Send NOOP if the connection has been idle longer than *interval* seconds.

    :param ftp: active FTP connection
    :param last_activity: monotonic timestamp of last FTP activity
    :param interval: seconds of idle time before sending NOOP
    :return: updated last-activity timestamp
    """
    if time.monotonic() - last_activity > interval:
        with contextlib.suppress(Exception):
            ftp.sendcmd("NOOP")
        return time.monotonic()
    return last_activity


def ftp_list_dir(ftp: FTP, path: PurePosixPath, retries: int = 3) -> list[PurePosixPath]:
    """List files in an FTP directory with retry on transient errors.

    :param ftp: active FTP connection
    :param path: remote directory path
    :param retries: number of retry attempts
    :return: list of filenames
    """
    ftp.cwd(str(path))

    @retry(
        retry=retry_if_exception_type(error_temp),
        stop=stop_after_attempt(retries),
        wait=wait_fixed(2),
        reraise=True,
        before_sleep=before_sleep_log(logger, logging.WARNING),
    )
    def _list() -> list[PurePosixPath]:
        file_strings: list[str] = []
        ftp.retrlines("NLST", file_strings.append)
        return [PurePosixPath(f) for f in file_strings]

    return _list()


def ftp_download_file(ftp: FTP, remote_path: PurePosixPath, local_path: Path, retries: int = 3) -> None:
    """Download a single file from FTP with retry on transient errors.

    :param ftp: active FTP connection
    :param remote_path: full remote file path
    :param local_path: local destination path
    :param retries: number of retry attempts
    """

    @retry(
        retry=retry_if_exception_type(error_temp),
        stop=stop_after_attempt(retries),
        wait=wait_fixed(2),
        reraise=True,
        before_sleep=before_sleep_log(logger, logging.WARNING),
    )
    def _download() -> None:
        with local_path.open("wb") as f:
            ftp.retrbinary(f"RETR {remote_path}", f.write)

    _download()


def ftp_retrieve_text(ftp: FTP, remote_path: PurePosixPath) -> str:
    """Retrieve a text file from FTP, returning its content as a string.

    :param ftp: active FTP connection
    :param remote_path: full remote file path
    :return: file content
    """
    lines: list[str] = []
    ftp.retrlines(f"RETR {remote_path}", lines.append)
    return "\n".join(lines)


class ThreadLocalFTP:
    """Manage thread-local FTP connections for parallel downloads.

    Each thread gets its own FTP connection, created on first access.
    Call :meth:`close_all` when done to cleanly shut down all connections.
    """

    def __init__(self, host: str, timeout: int = DEFAULT_TIMEOUT) -> None:
        """Initialise with FTP host and timeout.

        :param host: FTP hostname (required — no default)
        :param timeout: connection timeout in seconds
        """
        self._host = host
        self._timeout = timeout
        self._local = threading.local()
        self._lock = threading.Lock()
        self._connections: list[FTP] = []

    def get(self) -> FTP:
        """Return the FTP connection for the current thread, reconnecting if stale.

        Sends a NOOP to verify the connection is still alive before returning it.
        If the server has closed the connection (e.g. after an idle timeout or
        session limit), the dead socket is discarded and a fresh connection is
        established transparently.
        """
        ftp = getattr(self._local, "ftp", None)
        if ftp is not None:
            try:
                ftp.voidcmd("NOOP")
            except Exception:  # noqa: BLE001
                # Connection is stale — discard it and reconnect below
                with contextlib.suppress(Exception):
                    ftp.quit()
                with self._lock, contextlib.suppress(ValueError):
                    self._connections.remove(ftp)
                ftp = None
                self._local.ftp = None
        if ftp is None:
            ftp = connect_ftp(self._host, self._timeout)
            self._local.ftp = ftp
            with self._lock:
                self._connections.append(ftp)
        return ftp

    def close_all(self) -> None:
        """Close all thread-local FTP connections."""
        with self._lock:
            for ftp in self._connections:
                with contextlib.suppress(Exception):
                    ftp.quit()
            self._connections.clear()
