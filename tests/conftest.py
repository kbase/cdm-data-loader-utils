"""Global configuration settings for tests."""

import pytest


@pytest.fixture(scope="session")
def genome_paths_file() -> str:
    """Input file for tests."""
    return "tests/data/genome_paths.json"
