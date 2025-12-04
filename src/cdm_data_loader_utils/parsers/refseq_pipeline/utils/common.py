from typing import Iterable, TypeVar

T = TypeVar("T")


def chunks(iterable: list[T], size: int) -> Iterable[list[T]]:
    """
    Yield successive chunks of a list.
    """
    for i in range(0, len(iterable), size):
        yield iterable[i : i + size]


def ceildiv(a: int, b: int) -> int:
    """
    Ceiling division. Returns ⌈a / b⌉
    """
    return -(-a // b)
