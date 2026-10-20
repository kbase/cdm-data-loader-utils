"""Shared fixtures and functions for cdm_data_loaders.readers."""

from collections.abc import Callable
from pathlib import Path

import pytest
import xmlschema


@pytest.fixture
def load_schema(test_data_dir: Path) -> Callable[[str], xmlschema.XMLSchema]:
    """Return a function that builds an XMLSchema from a file name in the data directory."""

    def _load(filename: str) -> xmlschema.XMLSchema:
        return xmlschema.XMLSchema(str(test_data_dir / "xsd" / filename))

    return _load
