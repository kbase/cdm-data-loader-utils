"""Shared fixtures for xmltopyarrow pipeline tests."""

import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.core.settings import CtsSettings
from cdm_data_loaders.pipelines.xmltopyarrow.settings import XmlToPyArrowSettings
from tests.cdm_data_loaders.pipelines.helpers import make_cts_settings
from tests.xml_samples import LIBRARY_XSD


@pytest.fixture
def library_xsd(tmp_path: Path) -> Path:
    """Write LIBRARY_XSD to tmp_path and return its path."""
    xsd_path = tmp_path / "library.xsd"
    xsd_path.write_text(LIBRARY_XSD, encoding="utf-8")
    return xsd_path


@pytest.fixture
def settings_factory(tmp_path: Path, library_xsd: Path) -> Callable[..., Any]:
    """Return a factory that builds a valid XmlToPyArrowSettings (or CtsSettings), with fields open to override."""

    def _factory(
        settings_cls: type[CtsSettings] = XmlToPyArrowSettings,
        **overrides: Any,  # noqa: ANN401
    ) -> Any:  # noqa: ANN401
        defaults: dict[str, Any] = {
            "dlt_dev_mode": True,
        }
        if issubclass(settings_cls, XmlToPyArrowSettings):
            defaults["file_glob"] = "*.xml*"
        if settings_cls is XmlToPyArrowSettings:
            defaults.update(
                {
                    "dataset_name": "xml_test_dataset",
                    "table_name": "book",
                    "xml_tag": "book",
                    "xsd_file": str(library_xsd),
                }
            )
        defaults.update(overrides)
        return make_cts_settings(settings_cls, tmp_path, **defaults)

    return _factory


@pytest.fixture
def copy_xsd_fixture(tmp_path: Path, test_data_dir: Path) -> Callable[[str], str]:
    """Return a function that copies an XSD from tests/data/xsd into the input dir and returns its path."""

    def _copy(filename: str) -> str:
        input_dir = tmp_path / "input"
        input_dir.mkdir(exist_ok=True)
        dest = input_dir / filename
        shutil.copy(test_data_dir / "xsd" / filename, dest)
        return str(dest)

    return _copy
