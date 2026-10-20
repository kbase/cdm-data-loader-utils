"""Generic nested dictionary path operations."""

from collections.abc import Sequence
from typing import Any

from cdm_data_loaders.converters.core.errors import ConversionError


class NestedPathError(ConversionError):
    """Raised for an invalid nested path or a non-dictionary intermediate value."""


def set_nested(target: dict[str, Any], path: Sequence[str], value: object) -> None:
    """Set a leaf through a nonempty sequence of string keys, creating missing dicts.

    Existing leaves are replaced, including dicts and lists. Empty string keys are
    valid. Intermediate values must be dicts; collisions fail without mutation.
    The assigned value is not copied. A string alone is not a path sequence.
    """
    if not path or isinstance(path, (str, bytes)) or not all(isinstance(fragment, str) for fragment in path):
        message = "Path must be a nonempty sequence of string keys."
        raise NestedPathError(message)
    current = target
    for fragment in path[:-1]:
        if fragment not in current:
            current[fragment] = {}
        child = current[fragment]
        if not isinstance(child, dict):
            message = f"Path collision at {fragment!r}: intermediate value must be a dict."
            raise NestedPathError(message)
        current = child
    current[path[-1]] = value
