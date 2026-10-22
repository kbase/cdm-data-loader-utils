"""Shared fixtures for jsonlines pipeline tests."""

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.pipelines.jsonlines.settings import JsonlExtractSettings, JsonlPydanticIngestSettings
from tests.helpers import make_cts_settings


def make_item(
    record: object = None,
    line_no: int = 1,
    raw_record: str | None = None,
    parse_error: str | None = None,
    source_file: str = "file.jsonl",
) -> dict[str, Any]:
    """Build a line record matching stream_jsonl_lines' output shape."""
    if raw_record is None:
        raw_record = "" if parse_error is not None and record is None else json.dumps(record, separators=(",", ":"))
    return {
        "record": record,
        "raw_record": raw_record,
        "parse_error": parse_error,
        "source_file": source_file,
        "line_no": line_no,
    }


@pytest.fixture
def settings_factory(
    tmp_path: Path, dlt_destination_config: str, entity_models_module: str
) -> Callable[..., JsonlPydanticIngestSettings]:
    """Return a factory that builds a valid JsonlPydanticIngestSettings, with fields open to override."""

    def _factory(**overrides: Any) -> JsonlPydanticIngestSettings:  # noqa: ANN401
        defaults: dict[str, Any] = {
            "dataset_name": "test_dataset",
            "dlt_dev_mode": True,
            "entity_models_module": entity_models_module,
            "use_destination": dlt_destination_config,
        }
        defaults.update(overrides)
        return make_cts_settings(JsonlPydanticIngestSettings, tmp_path, **defaults)

    return _factory


@pytest.fixture
def extract_settings_factory(tmp_path: Path, dlt_destination_config: str) -> Callable[..., JsonlExtractSettings]:
    """Return a factory that builds a valid JsonlExtractSettings, with fields open to override."""

    def _factory(**overrides: Any) -> JsonlExtractSettings:  # noqa: ANN401
        defaults: dict[str, Any] = {
            "dataset_name": "test_dataset",
            "dlt_dev_mode": True,
            "table_name": "widget",
            "use_destination": dlt_destination_config,
        }
        defaults.update(overrides)
        return make_cts_settings(JsonlExtractSettings, tmp_path, **defaults)

    return _factory
