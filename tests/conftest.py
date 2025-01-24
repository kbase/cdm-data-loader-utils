"""Global configuration settings for tests."""

from pathlib import Path
from typing import Any

import pytest

stats_output = {
    "FW305-3-2-15-C-TSA1_scaffolds.fna": {
        "scaffolds": 59,
        "contigs": 61,
        "gc_avg": 0.65985,
        "gc_std": 0.02052,
        "filename": "/home/runner/work/cdm-data-loader-utils/cdm-data-loader-utils/tests/data/FW305-3-2-15-C-TSA1/FW305-3-2-15-C-TSA1_scaffolds.fna",
    },
    "FW305-C-112.1_scaffolds.fna": {
        "scaffolds": 28,
        "contigs": 31,
        "gc_avg": 0.69368,
        "gc_std": 0.03501,
        "filename": "/home/runner/work/cdm-data-loader-utils/cdm-data-loader-utils/tests/data/FW305-C-112.1/FW305-C-112.1_scaffolds.fna",
    },
}

checkm2_output = {
    "FW305-3-2-15-C-TSA1_scaffolds": {
        "checkm2_completeness": 99.99,
        "checkm2_contamination": 1.09,
    },
    "FW305-C-112.1_scaffolds": {"checkm2_completeness": None, "checkm2_contamination": 12.27},
}

RESULTS = {
    "no_file": {
        "stats": {},
        "checkm2": {},
    },
    "empty": {
        "stats": {},
        "checkm2": {},
    },
    "single": {
        "stats": {
            "FW305-3-2-15-C-TSA1_scaffolds.fna": stats_output["FW305-3-2-15-C-TSA1_scaffolds.fna"]
        },
        "checkm2": {
            "FW305-3-2-15-C-TSA1_scaffolds": checkm2_output["FW305-3-2-15-C-TSA1_scaffolds"]
        },
    },
    "multi": {"stats": stats_output, "checkm2": checkm2_output},
}


@pytest.fixture(scope="session")
def checkm2_stats_results() -> dict[str, Any]:
    """Expected results from parsing checkm2 and bbmap stats files."""
    return RESULTS


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Test data directory."""
    return Path("tests") / "data"


@pytest.fixture(scope="session")
def genome_paths_file() -> str:
    """Input file for tests."""
    return "tests/data/genome_paths.json"
