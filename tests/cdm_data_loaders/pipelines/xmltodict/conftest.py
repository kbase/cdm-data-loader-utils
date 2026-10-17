"""Shared fixtures for xml pipeline tests."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings
from tests.conftest import LOCAL_FS


@pytest.fixture
def settings_factory(tmp_path: Path) -> Callable[..., XmlToDictSettings]:
    """Return a factory that builds a valid XmlToDictSettings, with fields open to override."""

    def _factory(**overrides: Any) -> XmlToDictSettings:  # noqa: ANN401
        log_config_file = tmp_path / "logging.json"
        log_config_file.write_text('{"version": 1}')
        input_dir = tmp_path / "input"
        input_dir.mkdir(exist_ok=True)
        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        kwargs: dict[str, Any] = {
            "log_config_file": str(log_config_file),
            "input_dir": str(input_dir),
            "output_dir": str(output_dir),
            "dlt_dev_mode": True,
            "use_destination": LOCAL_FS,
            "use_output_dir_for_pipeline_metadata": False,
            "buffer_size": 100,
            "log_interval": 1000,
            "dataset_name": "xml_test_dataset",
            "table_name": "book",
            "xml_tag": "book",
            "file_glob": "*.xml*",
        }
        kwargs.update(overrides)
        return XmlToDictSettings(**kwargs)

    return _factory
