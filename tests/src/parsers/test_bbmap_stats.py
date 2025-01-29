"""Tests for the bbmap stats parser."""

import re
from pathlib import Path
from typing import Any

import pytest

from src.parsers.bbmap_stats import get_bbmap_stats
from tests.conftest import RESULTS


def test_get_bbmap_stats_empty(tmp_path: Path) -> None:
    """Test that an error is thrown by a non-existent file."""
    full_path = tmp_path / "some-file"
    with pytest.raises(
        RuntimeError,
        match=re.escape("error parsing stats_file: [Errno 2] No such file or directory"),
    ):
        get_bbmap_stats(full_path)


def test_get_bbmap_stats_empty_file(tmp_path: Path) -> None:
    """Test that an empty file throws an error."""
    full_path = tmp_path / "some-file"
    full_path.touch()

    with pytest.raises(
        RuntimeError,
        match=re.escape("no valid data found in stats_file"),
    ):
        get_bbmap_stats(full_path)


err_msgs = {
    "wrong_format": "stats_file is not in the correct format",
    "not_dicts": "stats_file is not in the correct format: stats must be in dictionary form",
    "invalid_entry": "invalid entry format",
}

format_errors = {
    # top level data struct is invalid JSON
    "unclosed_str": "error parsing stats_file: Unterminated string starting at: line 1 column 10 (char 9)",
    "null_key": "error parsing stats_file: Expecting property name enclosed in double quotes: line 1 column 11 (char 10)",
    "quoted_ws": "Invalid control character at: line 1 column 3 (char 2)",
    "str": "error parsing stats_file: Expecting value: line 1 column 2 (char 1)",
    # not an array of dicts
    "null": err_msgs["not_dicts"],
    "quoted_str": err_msgs["not_dicts"],
    "empty_str": err_msgs["not_dicts"],
    "empty_array": err_msgs["not_dicts"],
    "array_null": err_msgs["not_dicts"],
    "array_of_str": err_msgs["not_dicts"],
    "array_mixed": err_msgs["not_dicts"],
    "array_of_arrays": err_msgs["not_dicts"],
    # errors with the dict content
    "object": "invalid entry format",
    # no data
    "ws": "no valid data found in stats_file",
    "empty": "no valid data found in stats_file",
}

err_types = {
    err_id: RuntimeError
    for err_id in format_errors
    if format_errors[err_id] != err_msgs["not_dicts"]
}


@pytest.mark.parametrize(
    "params",
    [
        pytest.param(
            {
                "err_type": err_types.get(err_id, TypeError),
                "err_msg": format_errors.get(err_id, err_msgs["not_dicts"]),
                "input": err_id,
            },
            id=err_id,
        )
        for err_id in format_errors
    ],
)
def test_get_bbmap_stats_no_valid_data(
    params: dict[str, str],
    test_data_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    json_test_strings: dict[str, Any],
) -> None:
    """Test that invalid JSON structures throw an error."""
    full_path = test_data_dir / "some-file"

    def mockreturn(_) -> dict[str, Any] | list[str | Any]:
        return json_test_strings[params["input"]]

    # patch Path.read_text to return the data structure in params
    monkeypatch.setattr(Path, "read_text", mockreturn)

    with pytest.raises(
        params["err_type"],
        match=re.escape(params["err_msg"]),
    ):
        get_bbmap_stats(full_path)


@pytest.mark.parametrize("key", list(RESULTS))
def test_bbmap_stats(key: str, test_data_dir: Path, checkm2_stats_results: dict[str, Any]) -> None:
    """Test that stats are correctly parsed."""
    results_dir = test_data_dir / f"results_{key}" / "stats.json"
    assert get_bbmap_stats(results_dir) == checkm2_stats_results[key]["stats"]
