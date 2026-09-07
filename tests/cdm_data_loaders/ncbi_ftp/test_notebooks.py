"""Smoke tests for NCBI FTP notebooks — syntax and import validation."""

import ast
import json
from pathlib import Path, PurePosixPath

import pytest

from cdm_data_loaders.ncbi_ftp.assembly import FTP_HOST
from cdm_data_loaders.ncbi_ftp.manifest import (
    AssemblyRecord,
    compute_diff,
    download_assembly_summary,
    write_updated_manifest,
)
from cdm_data_loaders.pipelines.ncbi_ftp_download import (
    DEFAULT_STAGING_KEY_PREFIX,
    download_and_stage,
)

NOTEBOOKS_DIR = Path(__file__).resolve().parents[2] / "notebooks"

NCBI_NOTEBOOKS = [
    "ncbi_ftp_manifest.ipynb",
    "ncbi_ftp_download.ipynb",
]


def _extract_code_cells(notebook_path: Path) -> list[str]:
    """Extract source code from all code cells in a notebook.

    :param notebook_path: path to the .ipynb file
    :return: list of source code strings, one per code cell
    """
    with notebook_path.open() as f:
        nb = json.load(f)
    return ["".join(cell.get("source", [])) for cell in nb.get("cells", []) if cell.get("cell_type") == "code"]


@pytest.mark.parametrize("notebook", NCBI_NOTEBOOKS)
def test_notebook_syntax(notebook: str) -> None:
    """Every code cell is syntactically valid Python and non-empty."""
    path = NOTEBOOKS_DIR / notebook
    assert path.exists(), f"Notebook not found: {path}"
    cells = _extract_code_cells(path)
    assert len(cells) > 0, f"No code cells found in {notebook}"
    for i, source in enumerate(cells, 1):
        assert source.strip(), f"{notebook} cell {i} is empty"
        try:
            ast.parse(source, filename=f"{notebook}:cell{i}")
        except SyntaxError as exc:
            pytest.fail(f"{notebook} cell {i} has a syntax error: {exc}")


def test_manifest_notebook_imports() -> None:
    """All manifest notebook imports are verified at module load time above."""
    assert isinstance(FTP_HOST, str)
    assert FTP_HOST
    assert AssemblyRecord is not None
    assert callable(download_assembly_summary)
    assert callable(compute_diff)
    assert callable(write_updated_manifest)


def test_download_notebook_imports() -> None:
    """All download notebook imports resolve without error."""
    assert callable(download_and_stage)
    assert isinstance(DEFAULT_STAGING_KEY_PREFIX, PurePosixPath)
