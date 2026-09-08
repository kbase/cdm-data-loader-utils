"""Unit tests for progress.py."""

import re
import threading
import time

import pytest
from tqdm import tqdm

from cdm_data_loaders.utils.progress import (
    DEFAULT_UNIT,
    DEFAULT_UNIT_DIVISOR,
    SynchronizedCallback,
    make_progress_bar,
)

"""Tests for SynchronizedCallback"""


def test_synchronizedcallback_pass_forwards_single_call() -> None:
    """Calling the wrapper invokes the wrapped callback exactly once with the given argument."""
    calls: list[int] = []
    wrapped = SynchronizedCallback(calls.append)
    wrapped(42)
    assert calls == [42]


def test_synchronizedcallback_pass_forwards_multiple_calls_in_order() -> None:
    """Sequential calls are forwarded to the wrapped callback in the order they were made."""
    calls: list[int] = []
    wrapped = SynchronizedCallback(calls.append)
    for amount in (1, 2, 3, 4):
        wrapped(amount)
    assert calls == [1, 2, 3, 4]


def test_synchronizedcallback_pass_works_with_lambda() -> None:
    """SynchronizedCallback works with a lambda, not just named functions."""
    n1 = 5
    n2 = 7
    total = {"value": 0}
    wrapped = SynchronizedCallback(lambda n: total.__setitem__("value", total["value"] + n))
    wrapped(n1)
    wrapped(n2)
    assert total["value"] == n1 + n2


def test_synchronizedcallback_pass_works_with_bound_method() -> None:
    """SynchronizedCallback works with a bound method as the wrapped callable."""
    n1 = 5
    n2 = 7

    class Accumulator:
        def __init__(self) -> None:
            self.total = 0

        def add(self, amount: int) -> None:
            self.total += amount

    acc = Accumulator()
    wrapped = SynchronizedCallback(acc.add)
    wrapped(n1)
    wrapped(n2)
    assert acc.total == n1 + n2


def test_synchronizedcallback_pass_thread_safe_concurrent_calls_all_counted() -> None:
    """Many concurrent calls from multiple threads all reach the wrapped callback exactly once each.

    Uses a plain (non-atomic) int increment, which is prone to lost updates under a
    race condition, to prove the lock actually serialises access.
    """
    counter = {"value": 0}

    def unsafe_increment(amount: int) -> None:
        # read-modify-write with no atomicity of its own; relies entirely on the
        # SynchronizedCallback's lock for correctness under concurrency
        current = counter["value"]
        current += amount
        counter["value"] = current

    wrapped = SynchronizedCallback(unsafe_increment)

    n_threads = 20
    calls_per_thread = 50

    def worker() -> None:
        for _ in range(calls_per_thread):
            wrapped(1)

    threads = [threading.Thread(target=worker) for _ in range(n_threads)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    assert counter["value"] == n_threads * calls_per_thread


def test_synchronizedcallback_pass_lock_serialises_overlapping_calls() -> None:
    """A slow callback call blocks a concurrent second call until the first completes.

    Demonstrates that the internal lock provides real mutual exclusion, not just
    "correct by luck" behaviour.
    """
    events: list[str] = []
    entered = threading.Event()
    release = threading.Event()

    def slow_callback(n: int) -> None:
        events.append(f"enter-{n}")
        entered.set()
        # hold the lock until the main thread says we can proceed
        release.wait(timeout=5)
        events.append(f"exit-{n}")

    wrapped = SynchronizedCallback(slow_callback)

    thread = threading.Thread(target=wrapped, args=(1,))
    thread.start()
    assert entered.wait(timeout=5)

    # second call, made from the main thread, must block until `thread` releases the lock
    second_call_done = threading.Event()

    def second_call() -> None:
        wrapped(2)
        second_call_done.set()

    second_thread = threading.Thread(target=second_call)
    second_thread.start()

    # give the second thread a moment to attempt (and be blocked on) acquiring the lock
    time.sleep(0.1)
    assert not second_call_done.is_set()

    release.set()
    thread.join(timeout=5)
    second_thread.join(timeout=5)

    assert second_call_done.is_set()
    assert events == ["enter-1", "exit-1", "enter-2", "exit-2"]


def test_synchronizedcallback_fail_exception_propagates_and_lock_is_released() -> None:
    """An exception raised inside the wrapped callback propagates to the caller.

    Also verifies the lock is correctly released afterward (no deadlock on the next call).
    """

    def failing_callback(amount: int) -> None:
        msg = f"boom: {amount}"
        raise ValueError(msg)

    wrapped = SynchronizedCallback(failing_callback)

    with pytest.raises(ValueError, match="boom: 1"):
        wrapped(1)

    # the lock must have been released; a subsequent call should not hang
    with pytest.raises(ValueError, match="boom: 2"):
        wrapped(2)


"""Tests for make_progress_bar"""


def test_make_progress_bar_pass_returns_tqdm_instance() -> None:
    """make_progress_bar always returns a real tqdm instance."""
    bar = make_progress_bar(total=100, desc="test")
    try:
        assert isinstance(bar, tqdm)
    finally:
        bar.close()


def test_make_progress_bar_pass_disable_false_by_default() -> None:
    """When disable is not specified, the returned bar is enabled (not suppressed)."""
    bar = make_progress_bar(total=100, desc="test")
    try:
        assert bar.disable is False
    finally:
        bar.close()


def test_make_progress_bar_pass_disable_true_still_returns_usable_instance() -> None:
    """disable=True returns a real, usable tqdm instance rather than None."""
    bar = make_progress_bar(total=100, desc="test", disable=True)
    try:
        assert isinstance(bar, tqdm)
        assert bar.disable is True
    finally:
        bar.close()


def test_make_progress_bar_pass_total_none_is_accepted() -> None:
    """total=None is accepted for cases where the size is unknown (e.g. missing Content-Length)."""
    bar = make_progress_bar(total=None, desc="test")
    try:
        assert bar.total is None
    finally:
        bar.close()


def test_make_progress_bar_pass_total_int_is_stored() -> None:
    """A numeric total is stored on the returned bar."""
    total = 12345
    bar = make_progress_bar(total=total, desc="test")
    try:
        assert bar.total == total
    finally:
        bar.close()


def test_make_progress_bar_pass_desc_is_set() -> None:
    """The desc argument is applied to the returned bar."""
    bar = make_progress_bar(total=10, desc="my-transfer")
    try:
        assert bar.desc == "my-transfer"
    finally:
        bar.close()


def test_make_progress_bar_pass_default_unit_settings() -> None:
    """Default unit, unit_scale, and unit_divisor match the documented defaults."""
    bar = make_progress_bar(total=10, desc="test")
    try:
        assert bar.unit == DEFAULT_UNIT
        assert bar.unit_scale is True
        assert bar.unit_divisor == DEFAULT_UNIT_DIVISOR
    finally:
        bar.close()


def test_make_progress_bar_pass_custom_unit_settings_respected() -> None:
    """Custom overrides for unit, unit_scale, and unit_divisor are applied instead of the defaults."""
    divisor = 1000
    bar = make_progress_bar(
        total=10,
        desc="test",
        unit="it",
        unit_scale=False,
        unit_divisor=divisor,
    )
    try:
        assert bar.unit == "it"
        assert bar.unit_scale is False
        assert bar.unit_divisor == divisor
    finally:
        bar.close()


def test_make_progress_bar_pass_usable_as_context_manager() -> None:
    """The returned bar can be used as a context manager, as required by its documented usage."""
    n = 5
    with make_progress_bar(total=10, desc="test") as bar:
        assert isinstance(bar, tqdm)
        bar.update(n)
        assert bar.n == n


@pytest.mark.parametrize("disable", [True, False])
def test_make_progress_bar_pass_update_advances_n_regardless_of_disable(disable: bool) -> None:
    """update() only advances the internal counter (.n) if disable=False."""
    bar = make_progress_bar(total=10, desc="test", disable=disable)
    try:
        bar.update(3)
        assert bar.n == 0 if disable else 3
        bar.update(4)
        assert bar.n == 0 if disable else 7
    finally:
        bar.close()


def test_make_progress_bar_fail_requires_keyword_arguments() -> None:
    """All parameters are keyword-only; calling with positional arguments raises TypeError."""
    with pytest.raises(TypeError, match=r"make_progress_bar\(\) takes 0 positional arguments but 2 were given"):
        make_progress_bar(100, "test")  # type: ignore[misc,call-arg]


@pytest.mark.parametrize("args", [{}, {"total": 10}, {"desc": "test"}])
def test_make_progress_bar_fail_missing_required_arguments(args: dict[str, int | str]) -> None:
    """Keywords 'total' and 'desc' are required; omitting them raises TypeError."""
    regex = re.compile(r"make_progress_bar\(\) missing [12] required keyword-only arguments?")
    with pytest.raises(TypeError, match=regex):
        make_progress_bar(**args)  # type: ignore[call-arg]


def test_make_progress_bar_pass_independent_instances() -> None:
    """Each call to make_progress_bar returns a distinct, independently-tracked instance."""
    n = 5
    bar_a = make_progress_bar(total=10, desc="a")
    bar_b = make_progress_bar(total=10, desc="b")
    try:
        bar_a.update(n)
        assert bar_a.n == n
        assert bar_b.n == 0
        assert bar_a is not bar_b
    finally:
        bar_a.close()
        bar_b.close()
