"""Tests for the genome paths file parser."""

import json
import re
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loader_utils.parsers.genome_paths import get_genome_paths

GPF_DIR = "genome_paths_file"


def test_get_genome_paths_empty(tmp_path: Path) -> None:
    """Test that an error is thrown by a non-existent file."""
    full_path = tmp_path / "some-file"
    with pytest.raises(
        RuntimeError,
        match=re.escape("error parsing genome_paths_file: [Errno 2] No such file or directory"),
    ):
        get_genome_paths(full_path)


def test_get_genome_paths_empty_file(tmp_path: Path) -> None:
    """Test that an empty file throws an error."""
    full_path = tmp_path / "some-file"
    full_path.touch()

    with pytest.raises(
        RuntimeError,
        match=re.escape("error parsing genome_paths_file: Expecting value: line 1 column 1 (char 0)"),
    ):
        get_genome_paths(full_path)


format_errors = {
    # JSON parse errors
    "empty": "error parsing genome_paths_file: Expecting value: line 1 column 1 (char 0)",
    "ws": "error parsing genome_paths_file: Expecting value: line 4 column 5 (char 8)",
    "str": "error parsing genome_paths_file: Expecting value: line 1 column 1 (char 0)",
    "unclosed_str": "error parsing genome_paths_file: Unterminated string starting at: line 1 column 9 (char 8)",
    "null_key": "error parsing genome_paths_file: Expecting property name enclosed in double quotes: line 1 column 10 (char 9)",
    # wrong format (array)
    "array_of_objects": "genome_paths_file is not in the correct format",
    "empty_array": "genome_paths_file is not in the correct format",
    # no data
    "empty_object": "no valid data found in genome_paths_file",
}

err_types = {"empty_array": TypeError, "array_of_objects": TypeError}


@pytest.mark.parametrize(
    "params",
    [
        pytest.param(
            {
                "err_msg": format_errors[err_id],
                "input": err_id,
            },
            id=err_id,
        )
        for err_id in format_errors
    ],
)
def test_get_genome_paths_invalid_format(
    params: dict[str, str],
    test_data_dir: Path,
    monkeypatch: pytest.MonkeyPatch,
    json_test_strings: dict[str, Any],
) -> None:
    """Test that invalid JSON structures throw an error."""
    full_path = test_data_dir / GPF_DIR / "valid.json"

    def mockreturn(_) -> dict[str, Any] | list[str | Any]:
        return json.loads(json_test_strings[params["input"]])

    # patch json.load to return the data structure in params
    monkeypatch.setattr(json, "load", mockreturn)

    with pytest.raises(
        err_types.get(params["input"], RuntimeError),
        match=re.escape(params["err_msg"]),
    ):
        get_genome_paths(full_path)


error_list = {
    "no_entry": [{"": {"this": "that"}}, 'No ID specified for entry {"this": "that"}', ValueError],
    "invalid_entry_format_arr": [{"id": []}, "id: invalid entry format"],
    "invalid_entry_format_str": [{"id": "some string"}, "id: invalid entry format"],
    "invalid_entry_format_None": [{"id": None}, "id: invalid entry format"],
    "no_valid_paths": [{"id": {}}, "id: no valid file types or paths found"],
    "invalid_keys": [{"id": {"pap": 1, "pip": 2, "pop": 3}}, "id: invalid keys: pap, pip, pop"],
}


@pytest.mark.parametrize(
    "params",
    [
        pytest.param(
            {
                "err_msg": error_list[err_id][1],
                "input": error_list[err_id][0],
            },
            id=err_id,
        )
        for err_id in error_list
    ],
)
def test_get_genome_paths_valid_input_invalid_format(
    params: dict[str, str], test_data_dir: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Test that invalid JSON structures throw an error."""
    full_path = test_data_dir / GPF_DIR / "valid.json"

    def mockreturn(_) -> dict[str, Any] | list[str | Any]:
        return params["input"]

    # patch json.load to return the data structure in params
    monkeypatch.setattr(json, "load", mockreturn)

    with pytest.raises(
        RuntimeError,
        match=f"Please ensure that the genome_paths_file is in the correct format.\n\n{params['err_msg']}",
    ):
        get_genome_paths(full_path)


def test_get_genome_paths(test_data_dir: Path) -> None:
    """Test that the genome paths file can be correctly parsed."""
    full_path = test_data_dir / GPF_DIR / "valid.json"
    assert get_genome_paths(full_path) == {
        "FW305-3-2-15-C-TSA1.1": {
            "fna": "tests/data/FW305-3-2-15-C-TSA1/FW305-3-2-15-C-TSA1_scaffolds.fna",
            "gff": "tests/data/FW305-3-2-15-C-TSA1/FW305-3-2-15-C-TSA1_genes.gff",
            "protein": "tests/data/FW305-3-2-15-C-TSA1/FW305-3-2-15-C-TSA1_genes.faa",
        },
        "FW305-C-112.1": {
            "fna": "tests/data/FW305-C-112.1/FW305-C-112.1_scaffolds.fna",
            "gff": "tests/data/FW305-C-112.1/FW305-C-112.1_genes.gff",
            "protein": "tests/data/FW305-C-112.1/FW305-C-112.1_genes.faa",
        },
    }
