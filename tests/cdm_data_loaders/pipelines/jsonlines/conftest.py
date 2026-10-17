"""Shared fixtures for jsonlines pipeline tests."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.pipelines.jsonlines.settings import JsonlExtractSettings, JsonlPydanticIngestSettings


@pytest.fixture
def settings_factory(
    tmp_path: Path, dlt_destination_config: str, entity_models_module: str
) -> Callable[..., JsonlPydanticIngestSettings]:
    """Return a factory that builds a valid JsonlPydanticIngestSettings, with fields open to override."""

    def _factory(**overrides: Any) -> JsonlPydanticIngestSettings:  # noqa: ANN401
        log_config_file = tmp_path / "logging.conf"
        log_config_file.touch()
        input_dir = tmp_path / "input"
        input_dir.mkdir(exist_ok=True)
        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        kwargs: dict[str, Any] = {
            "dataset_name": "test_dataset",
            "dlt_dev_mode": True,
            "entity_models_module": entity_models_module,
            "input_dir": str(input_dir),
            "log_config_file": str(log_config_file),
            "output_dir": str(output_dir),
            "use_destination": dlt_destination_config,
            "use_output_dir_for_pipeline_metadata": False,
        }
        kwargs.update(overrides)
        return JsonlPydanticIngestSettings(**kwargs)

    return _factory


@pytest.fixture
def extract_settings_factory(tmp_path: Path, dlt_destination_config: str) -> Callable[..., JsonlExtractSettings]:
    """Return a factory that builds a valid JsonlExtractSettings, with fields open to override."""

    def _factory(**overrides: Any) -> JsonlExtractSettings:  # noqa: ANN401
        log_config_file = tmp_path / "logging.conf"
        log_config_file.touch()
        input_dir = tmp_path / "input"
        input_dir.mkdir(exist_ok=True)
        output_dir = tmp_path / "output"
        output_dir.mkdir(exist_ok=True)

        kwargs: dict[str, Any] = {
            "dataset_name": "test_dataset",
            "dlt_dev_mode": True,
            "input_dir": str(input_dir),
            "log_config_file": str(log_config_file),
            "output_dir": str(output_dir),
            "table_name": "widget",
            "use_destination": dlt_destination_config,
            "use_output_dir_for_pipeline_metadata": False,
        }
        kwargs.update(overrides)
        return JsonlExtractSettings(**kwargs)

    return _factory
