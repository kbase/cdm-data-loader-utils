"""Shared fixtures for xml pipeline tests."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings
from tests.cdm_data_loaders.pipelines.helpers import make_cts_settings


@pytest.fixture
def settings_factory(tmp_path: Path) -> Callable[..., XmlToDictSettings]:
    """Return a factory that builds a valid XmlToDictSettings, with fields open to override."""

    def _factory(**overrides: Any) -> XmlToDictSettings:  # noqa: ANN401
        defaults: dict[str, Any] = {
            "dlt_dev_mode": True,
            "dataset_name": "xml_test_dataset",
            "table_name": "book",
            "xml_tag": "book",
            "file_glob": "*.xml*",
        }
        defaults.update(overrides)
        return make_cts_settings(XmlToDictSettings, tmp_path, **defaults)

    return _factory

