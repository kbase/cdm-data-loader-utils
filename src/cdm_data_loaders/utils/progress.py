"""Progress bar helpers built on tqdm, usable across download and upload paths."""

import threading
from collections.abc import Callable
from typing import Any, Final

from tqdm import tqdm

DEFAULT_UNIT: Final[str] = "B"
DEFAULT_UNIT_DIVISOR: Final[int] = 1024


class SynchronizedCallback:
    """Wraps a callable so it can be safely invoked from multiple threads.

    boto3's `Callback` parameter (used by `upload_fileobj`/`download_fileobj`)
    may be invoked concurrently from multiple transfer threads during
    multipart transfers. tqdm's `update()` is not guaranteed to be
    thread-safe, so any callback driving a tqdm bar should be wrapped with this.
    """

    def __init__(self, callback: Callable[[int], Any]) -> None:
        """Initialise a thread-safe wrapper around a single-argument callback.

        :param callback: callable that accepts the number of bytes
            transferred since the previous invocation
        :type callback: Callable[[int], Any]
        """
        self._callback = callback
        self._lock = threading.Lock()

    def __call__(self, bytes_amount: int) -> None:
        """Invoke the wrapped callback under a lock.

        :param bytes_amount: number of bytes transferred since the previous call
        :type bytes_amount: int
        """
        with self._lock:
            self._callback(bytes_amount)


def make_progress_bar(  # noqa: PLR0913
    *,
    total: int | None,
    desc: str,
    disable: bool = False,
    unit: str = DEFAULT_UNIT,
    unit_scale: bool = True,
    unit_divisor: int = DEFAULT_UNIT_DIVISOR,
) -> tqdm:
    """Create a tqdm progress bar, always returning a usable instance.

    When `disable` is True, a real tqdm instance is still returned (with all
    output suppressed) rather than None, so callers don't need branching
    logic to support the "no progress bar" case — it can always be used as
    a context manager and always be updated via `Callback`/`.update()`.

    :param total: total number of bytes expected, or None if unknown (e.g.
        the server did not return a Content-Length header)
    :type total: int | None
    :param desc: label shown alongside the progress bar
    :type desc: str
    :param disable: if True, suppress all output from the progress bar, defaults to False
    :type disable: bool, optional
    :param unit: unit label for the progress bar, defaults to "B"
    :type unit: str, optional
    :param unit_scale: whether to auto-scale the unit (e.g. KB/MB/GB), defaults to True
    :type unit_scale: bool, optional
    :param unit_divisor: divisor used when auto-scaling units, defaults to 1024
    :type unit_divisor: int, optional
    :return: a tqdm progress bar instance; intended to be used as a context manager
    :rtype: tqdm
    """
    return tqdm(
        total=total,
        desc=desc,
        disable=disable,
        unit=unit,
        unit_scale=unit_scale,
        unit_divisor=unit_divisor,
    )
