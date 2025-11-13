"""Tests for the checkm2 results parser."""

import re
from pathlib import Path

import pytest

from cdm_data_loader_utils.parsers.checkm2 import get_checkm2_data
from tests.conftest import RESULTS


def test_get_checkm2_data_empty(tmp_path: Path) -> None:
    """Test that an error is thrown by a non-existent file."""
    full_path = tmp_path / "some-file"
    with pytest.raises(
        RuntimeError,
        match=re.escape("error parsing checkm2_file: [Errno 2] No such file or directory"),
    ):
        get_checkm2_data(full_path)


def test_get_checkm2_data_empty_file(tmp_path: Path) -> None:
    """Test that an empty file throws an error."""
    full_path = tmp_path / "some-file"
    full_path.touch()

    with pytest.raises(
        RuntimeError,
        match=re.escape("error parsing checkm2_file: file is not in TSV format"),
    ):
        get_checkm2_data(full_path)


def test_get_checkm2_data_not_tsv(test_data_dir: Path) -> None:
    """Test that a file in the wrong format throws an error."""
    full_path = test_data_dir / "genome_paths_file" / "valid.json"

    with pytest.raises(
        RuntimeError,
        match=re.escape(
            "error parsing checkm2_file: checkm2 output is missing the following columns: Name, Completeness, Contamination"
        ),
    ):
        get_checkm2_data(full_path)


def test_get_checkm2_data_wrong_headers(test_data_dir: Path) -> None:
    """Test that a file with the wrong headers throws an error."""
    full_path = test_data_dir / "checkm2_file" / "wrong_headers.tsv"

    with pytest.raises(
        RuntimeError,
        match=re.escape(
            "error parsing checkm2_file: checkm2 output is missing the following columns: Name, Completeness, Contamination"
        ),
    ):
        get_checkm2_data(full_path)


def test_get_checkm2_data_headers_only(test_data_dir: Path) -> None:
    """Test that a file with only headers throws an error."""
    full_path = test_data_dir / "checkm2_file" / "headers_only.tsv"

    with pytest.raises(
        RuntimeError,
        match=re.escape("no valid data found in checkm2_file"),
    ):
        get_checkm2_data(full_path)


def test_get_checkm2_data_missing_names(test_data_dir: Path) -> None:
    """Ensure that a file with rows lacking a name throws an error."""
    full_path = test_data_dir / "checkm2_file" / "missing_names.tsv"
    with pytest.raises(
        RuntimeError,
        match=re.escape("errors found in checkm2_file:\nrow 2 has no Name value\nrow 3 has no Name value"),
    ):
        get_checkm2_data(full_path)


@pytest.mark.parametrize("key", list(RESULTS))
def test_get_checkm2_scores(key: str, test_data_dir: Path) -> None:
    """Test that checkm2 scores are correctly parsed."""
    results_dir = test_data_dir / f"results_{key}" / "quality_report.tsv"
    assert get_checkm2_data(results_dir) == RESULTS[key]["checkm2"]
