"""Tests for generic dictionary path assignment."""

from collections.abc import Sequence
from copy import deepcopy
from typing import Any

import pytest

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.paths import NestedPathError, set_nested


@pytest.mark.parametrize(
    ("target", "path", "expected"),
    [
        ({}, ["parent", "child"], {"parent": {"child": ["new"]}}),
        ({"parent": {"keep": 1}}, ("parent", "child"), {"parent": {"keep": 1, "child": ["new"]}}),
        ({"leaf": "old"}, ["leaf"], {"leaf": ["new"]}),
        ({"leaf": {"old": True}}, ["leaf"], {"leaf": ["new"]}),
        ({"leaf": ["old"]}, ["leaf"], {"leaf": ["new"]}),
        ({}, [""], {"": ["new"]}),
        ({}, ["branch", "branch", "leaf"], {"branch": {"branch": {"leaf": ["new"]}}}),
    ],
    ids=[
        "create-dicts",
        "preserve-siblings-tuple-path",
        "replace-scalar",
        "replace-dict",
        "replace-list",
        "empty-key",
        "repeated-key",
    ],
)
def test_set_nested_pass_assignment(target: dict[str, Any], path: Sequence[str], expected: dict[str, Any]) -> None:
    """Preserve existing siblings and replace leaves without copying the assigned value."""
    value = ["new"]
    assert set_nested(target, path, value) is None
    assert target == expected
    current = target
    for fragment in path[:-1]:
        current = current[fragment]
    assert current[path[-1]] is value


@pytest.mark.parametrize(
    ("target", "path", "message"),
    [
        ({}, [], "nonempty"),
        ({}, "leaf", "nonempty sequence"),
        ({}, b"leaf", "nonempty sequence"),
        ({}, ["new", 1], "string keys"),
        ({"branch": None}, ["branch", "leaf"], "collision"),
        ({"branch": []}, ["branch", "leaf"], "collision"),
        ({"branch": {"nested": 1}}, ["branch", "nested", "leaf"], "collision"),
    ],
    ids=[
        "empty-path",
        "bare-string",
        "bare-bytes",
        "invalid-key-before-mutation",
        "null-intermediate",
        "list-intermediate",
        "deep-scalar",
    ],
)
def test_set_nested_fail_invalid_path(target: dict[str, Any], path: Any, message: str) -> None:  # noqa: ANN401
    """Report path errors with a ConversionError subtype and preserve target state."""
    original = deepcopy(target)
    with pytest.raises(NestedPathError, match=message) as raised:
        set_nested(target, path, "new")
    assert isinstance(raised.value, ConversionError)
    assert target == original
