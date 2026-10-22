"""Fixtures shared by the XSV ingestion pipeline tests."""

import json
import os
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any, Final

import pytest

from cdm_data_loaders.pipelines.xsv.pipeline import XsvWorkPaths
from cdm_data_loaders.pipelines.xsv.settings import XsvIngestSettings
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import generate_header
from tests.cdm_data_loaders.pipelines.helpers import make_cts_settings

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


@pytest.fixture
def settings_factory(tmp_path: Path, dlt_destination_config: str) -> Callable[..., XsvIngestSettings]:
    """Return a factory that builds a valid XsvIngestSettings, with fields open to override."""

    def _factory(
        schema: dict[str, Any] | None = DEFAULT_XSV_SCHEMA,
        data_files: dict[str, str] | None = None,
        **overrides: Any,  # noqa: ANN401
    ) -> XsvIngestSettings:
        input_dir = Path(overrides.get("input_dir", tmp_path / "input"))
        input_dir.mkdir(parents=True, exist_ok=True)

        if schema is not None:
            (input_dir / "schema.json").write_text(json.dumps(schema), encoding="utf-8")
        if data_files:
            for file_name, content in data_files.items():
                (input_dir / file_name).write_text(content, encoding="utf-8")

        defaults: dict[str, Any] = {
            "dataset_name": "xsv_test_dataset",
            "file_glob": "*",
            "loader_file_format": "jsonl",
            "schema_file": str(input_dir / "schema.json"),
            "table_name": "my_table",
            "use_destination": dlt_destination_config,
        }
        defaults.update(overrides)
        return make_cts_settings(XsvIngestSettings, tmp_path, **defaults)

    return _factory


@pytest.fixture
def work_paths_factory(tmp_path: Path, qsv_cmd: str) -> Callable[[XsvIngestSettings], XsvWorkPaths]:
    """Return a factory that builds XsvWorkPaths (scratch dirs + header file) for one settings object."""

    def _factory(settings: XsvIngestSettings) -> XsvWorkPaths:
        tmp_dir = tmp_path / "tmp"
        qsv_output_dir = tmp_path / "qsv_output"
        validated_dir = tmp_path / "validated"
        for directory in (tmp_dir, qsv_output_dir, validated_dir):
            directory.mkdir(exist_ok=True)
        return XsvWorkPaths(
            qsv_cmd=qsv_cmd,
            header_file_path=generate_header(settings.validated_schema, tmp_dir),
            tmp_dir=tmp_dir,
            qsv_output_dir=qsv_output_dir,
            validated_dir=validated_dir,
        )

    return _factory


@pytest.fixture
def write_xsv_file() -> Callable[[XsvIngestSettings, str, str], Path]:
    """Return a function that writes text content to a new file inside settings.input_dir."""

    def _write(settings: XsvIngestSettings, file_name: str, content: str) -> Path:
        path = Path(settings.input_dir) / file_name
        path.write_text(content)
        return path

    return _write
