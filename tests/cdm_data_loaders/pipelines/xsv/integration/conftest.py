"""Shared fixtures for xsv pipeline integration tests."""

import json
import os
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, Final

import pytest

from cdm_data_loaders.pipelines.xsv.settings import XsvIngestSettings

VALID_SCHEMA_URI: Final[str] = "https://json-schema.org/draft/2020-12/schema"
COLUMNS: Final[list[str]] = ["number", "date", "float", "boolean", "string"]

# The `float` column is deliberately typed `["string", "null"]`, not `"number"`, matching the
# convention used throughout the xsv_validator test fixtures: it lets a field hold either a
# numeric-looking string or a null placeholder without needing float-parsing semantics.
DEFAULT_XSV_SCHEMA: Final[dict[str, Any]] = {
    "$schema": VALID_SCHEMA_URI,
    "title": "XSV integration test schema",
    "type": "object",
    "required": COLUMNS,
    "properties": {
        "number": {"type": "integer"},
        "date": {"type": "string"},
        "float": {"type": ["string", "null"]},
        "boolean": {"type": "boolean"},
        "string": {"type": "string"},
    },
    "x-xsv-config": {"x-delimiter": ",", "x-has-header": True},
}

VALID_ROWS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "3.14", "true", "key:value1"],
    ["3", "2023-02-20", "1.11", "false", "key:value2"],
]

# One ragged row (missing the trailing `string` field) among otherwise-valid rows: recoverable.
PARTIAL_RAGGED_ROWS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "3.14", "true", "key:value1"],
    ["3", "2023-02-20", "1.11", "false"],
    ["4", "2023-03-25", "2.22", "true", "key:value3"],
]

# Every row here is ragged (too few or too many fields): unrecoverable.
ALL_RAGGED_ROWS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "3.14", "true"],
    ["3", "2023-02-20", "1.11", "false", "key:value2", "extra"],
]


def rows_to_csv(rows: list[list[str]], header: list[str] = COLUMNS, delimiter: str = ",") -> str:
    """Render rows (plus header) as delimited text content, one row per line."""
    lines = [delimiter.join(header), *(delimiter.join(row) for row in rows)]
    return "\n".join(lines) + "\n"


@pytest.fixture(scope="session")
def qsv_cmd() -> str:
    """Provide the path/name of the qsv binary; override via the QSV_BIN env var if not on PATH.

    Tests using this fixture are skipped when the qsv binary is not available.
    """
    cmd = shutil.which("qsv") or os.environ.get("QSV_BIN")
    if not cmd or not Path(cmd).is_file():
        pytest.skip(f"qsv binary not found at {cmd!r}; set QSV_BIN to point at a real binary")
    return cmd


def make_settings(
    tmp_path: Path,
    dlt_destination_config: str,
    schema: dict[str, Any],
    **overrides: Any,
) -> dict[str, Any]:
    """Build a valid xsv ingest settings dict, with fields open to override."""
    log_config_file = overrides.pop("log_config_file", tmp_path / "logging.conf")
    Path(log_config_file).touch(exist_ok=True)
    output_dir = overrides.pop("output_dir", tmp_path / "output")
    Path(output_dir).mkdir(exist_ok=True)
    input_dir = overrides.pop("input_dir", tmp_path / "input")
    Path(input_dir).mkdir(exist_ok=True)
    (Path(input_dir) / "schema.json").write_text(json.dumps(schema))

    kwargs: dict[str, Any] = {
        "buffer_size": 10,
        "dataset_name": "xsv_test_dataset",
        "dlt_dev_mode": False,
        "file_glob": "*",
        "input_dir": str(input_dir),
        "loader_file_format": "jsonl",
        "log_config_file": str(log_config_file),
        "log_interval": 1000,
        "output_dir": str(output_dir),
        "schema_file": str(Path(input_dir) / "schema.json"),
        "table_name": "my_table",
        "use_destination": dlt_destination_config,
        "use_output_dir_for_pipeline_metadata": False,
    }
    kwargs.update(overrides)
    return kwargs


@pytest.fixture
def settings_factory(tmp_path: Path, dlt_destination_config: str) -> Callable[..., XsvIngestSettings]:
    """Return a factory that builds a valid XsvIngestSettings, with fields open to override."""

    def _factory(schema: dict[str, Any] = DEFAULT_XSV_SCHEMA, **overrides: Any) -> XsvIngestSettings:
        settings_dict = make_settings(tmp_path, dlt_destination_config, schema, **overrides)
        return XsvIngestSettings(**settings_dict)

    return _factory


@pytest.fixture
def write_xsv_file() -> Callable[[XsvIngestSettings, str, str], Path]:
    """Return a function that writes text content to a new file inside settings.input_dir."""

    def _write(settings: XsvIngestSettings, file_name: str, content: str) -> Path:
        path = Path(settings.input_dir) / file_name
        path.write_text(content)
        return path

    return _write
