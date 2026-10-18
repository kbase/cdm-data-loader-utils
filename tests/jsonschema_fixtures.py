"""Isolated registry and settings fixtures for JSON Schema pipeline tests."""

import sys
from collections.abc import Callable
from itertools import count
from pathlib import Path
from types import ModuleType
from typing import Any
from uuid import uuid4

import dlt
import pytest

from cdm_data_loaders.pipelines.jsonlines.settings import JsonlJsonschemaIngestSettings


@pytest.fixture
def widget_schema() -> dict[str, Any]:
    """Return an object schema with required, numeric, and format constraints."""
    return {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "required": ["widget_id", "count"],
        "properties": {
            "widget_id": {"type": "string", "minLength": 1},
            "count": {"type": "integer", "minimum": 0},
            "collected": {"type": "string", "format": "date"},
            "label": {"type": "string", "default": "not inserted"},
        },
    }


@pytest.fixture
def schema_module_factory(monkeypatch: pytest.MonkeyPatch) -> Callable[[object], str]:
    """Register uniquely named schema modules and remove them after each test."""

    def _register(registry: object) -> str:
        module_name = f"_test_entity_schemas_{uuid4().hex}"
        module = ModuleType(module_name)
        module.__dict__["ENTITY_SCHEMAS"] = registry
        monkeypatch.setitem(sys.modules, module_name, module)
        return module_name

    return _register


@pytest.fixture
def schema_module(schema_module_factory: Callable[[object], str], widget_schema: dict[str, Any]) -> str:
    """Provide a registry with two entities and an intentionally absent directory."""
    return schema_module_factory({"widget": widget_schema, "gadget": widget_schema, "missing": widget_schema})


@pytest.fixture
def schema_settings_factory(
    tmp_path: Path, dlt_destination_config: str, schema_module: str
) -> Callable[..., JsonlJsonschemaIngestSettings]:
    """Create settings with fresh output and metadata directories for each run."""
    run_counter = count()

    def _factory(**overrides: object) -> JsonlJsonschemaIngestSettings:
        run_index = next(run_counter)
        input_dir = tmp_path / "input"
        input_dir.mkdir(exist_ok=True)
        output_dir = tmp_path / f"output_{run_index}"
        output_dir.mkdir()
        log_config_file = tmp_path / "logging.json"
        log_config_file.write_text('{"version": 1, "disable_existing_loggers": false}', encoding="utf-8")
        kwargs = {
            "input_dir": str(input_dir),
            "output_dir": str(output_dir),
            "log_config_file": str(log_config_file),
            "use_destination": dlt_destination_config,
            "use_output_dir_for_pipeline_metadata": True,
            "dlt_dev_mode": False,
            "dataset_name": "schema_test",
            "file_glob": "*.jsonl*",
            "loader_file_format": "parquet",
            "buffer_size": 2,
            "schema_files_module": schema_module,
        }
        kwargs.update(overrides)
        settings = JsonlJsonschemaIngestSettings(**kwargs)
        dlt.config["destination.local_fs.bucket_url"] = str(settings.output_dir)
        return settings

    return _factory
