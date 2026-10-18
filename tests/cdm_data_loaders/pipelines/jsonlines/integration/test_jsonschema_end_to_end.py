"""End-to-end JSON Schema ingestion tests using real filesystem destinations."""

import gzip
import json
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from cdm_data_loaders.pipelines.jsonlines.extract_jsonschema_validate_pipeline import (
    cli,
    run_jsonlines_pipeline,
)
from cdm_data_loaders.pipelines.jsonlines.settings import JSONSCHEMA_PIPELINE_NAME, JsonlJsonschemaIngestSettings
from tests.cdm_data_loaders.pipelines.helpers import LOADER_FILE_FORMATS, parse_json_string, read_data_rows

pytestmark = pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)

WIDGET_RECORD = {"widget_id": "a", "count": 0, "nested": {"tags": ["one", "two"]}}
GADGET_RECORD = {"widget_id": "g", "count": 3}
TYPE_ERROR_RAW = '{"widget_id":"bad", "count":"2"}'
FORMAT_ERROR_RAW = '{"widget_id":"bad","count":1,"collected":"2023-02-29"}'
MALFORMED_RAW = "{broken"
PARSE_ERROR = "Expecting property name enclosed in double quotes: line 1 column 2 (char 1)"


@pytest.fixture
def schema_input_dir(tmp_path: Path) -> Path:
    """Create two entities containing plain, compressed, valid, and invalid inputs."""
    input_dir = tmp_path / "input"
    widget_dir = input_dir / "widget"
    widget_dir.mkdir(parents=True)
    widget_line = json.dumps(WIDGET_RECORD)
    (widget_dir / "plain.jsonl").write_text(
        f"\n{widget_line}\n{MALFORMED_RAW}\n{TYPE_ERROR_RAW}\n{widget_line}", encoding="utf-8"
    )
    with gzip.open(widget_dir / "compressed.jsonl.gz", "wt", encoding="utf-8") as stream:
        stream.write(f"{widget_line}\n{FORMAT_ERROR_RAW}\nnull")
    (widget_dir / "ignored.txt").write_text(MALFORMED_RAW, encoding="utf-8")
    gadget_dir = input_dir / "gadget"
    gadget_dir.mkdir()
    (gadget_dir / "records.jsonl").write_text(json.dumps(GADGET_RECORD), encoding="utf-8")
    return input_dir


def read_data_tables(settings: JsonlJsonschemaIngestSettings, dataset_name: str) -> dict[str, list[dict[str, Any]]]:
    """Read persisted rows without dlt metadata, preserving rejection text."""
    return read_data_rows(Path(settings.output_dir) / dataset_name, str(settings.loader_file_format))


@pytest.mark.parametrize("buffer_size", [1, 2, 10], ids=["single-row", "partial-pages", "single-page"])
def test_run_jsonlines_pipeline_pass_persists_valid_and_rejected_rows(
    schema_input_dir: Path,
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings],
    loader_file_format: str,
    buffer_size: int,
) -> None:
    """Mixed files load every row exactly once with intact nested data and rejection provenance."""
    settings = schema_settings_factory(
        input_dir=str(schema_input_dir), loader_file_format=loader_file_format, buffer_size=buffer_size
    )

    load_info = run_jsonlines_pipeline(settings)

    assert load_info is not None
    assert not load_info.has_failed_jobs
    assert load_info.pipeline.pipeline_name == JSONSCHEMA_PIPELINE_NAME
    assert load_info.dataset_name == "schema_test"
    tables = read_data_tables(settings, load_info.dataset_name)
    assert set(tables) == {"widget", "widget_rejected", "gadget"}
    assert [{key: parse_json_string(value) for key, value in row.items()} for row in tables["widget"]] == [
        WIDGET_RECORD
    ] * 3
    assert tables["gadget"] == [GADGET_RECORD]
    assert sorted(tables["widget_rejected"], key=lambda row: (row["source_file"], row["line_no"])) == [
        {
            "source_file": "compressed.jsonl.gz",
            "line_no": 2,
            "raw_record": FORMAT_ERROR_RAW,
            "error_detail": json.dumps(["'2023-02-29' is not a 'date'"]),
        },
        {
            "source_file": "compressed.jsonl.gz",
            "line_no": 3,
            "raw_record": "null",
            "error_detail": json.dumps(["None is not of type 'object'"]),
        },
        {"source_file": "plain.jsonl", "line_no": 3, "raw_record": MALFORMED_RAW, "error_detail": PARSE_ERROR},
        {
            "source_file": "plain.jsonl",
            "line_no": 4,
            "raw_record": TYPE_ERROR_RAW,
            "error_detail": json.dumps(["'2' is not of type 'integer'"]),
        },
    ]


def test_run_jsonlines_pipeline_pass_all_invalid(
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings], loader_file_format: str
) -> None:
    """An entirely invalid entity still produces its complete rejection table."""
    settings = schema_settings_factory(loader_file_format=loader_file_format, buffer_size=1)
    table_dir = Path(settings.input_dir) / "widget"
    table_dir.mkdir()
    (table_dir / "invalid.jsonl").write_text(f"{MALFORMED_RAW}\n{TYPE_ERROR_RAW}", encoding="utf-8")

    load_info = run_jsonlines_pipeline(settings)

    assert load_info is not None
    assert not load_info.has_failed_jobs
    assert read_data_tables(settings, load_info.dataset_name) == {
        "widget_rejected": [
            {"source_file": "invalid.jsonl", "line_no": 1, "raw_record": MALFORMED_RAW, "error_detail": PARSE_ERROR},
            {
                "source_file": "invalid.jsonl",
                "line_no": 2,
                "raw_record": TYPE_ERROR_RAW,
                "error_detail": json.dumps(["'2' is not of type 'integer'"]),
            },
        ]
    }


def test_cli_pass_loads_selected_entity(
    schema_input_dir: Path,
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings],
    loader_file_format: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The real CTS CLI accepts a registry shortcut, table selection, and the requested file format."""
    settings = schema_settings_factory(input_dir=str(schema_input_dir), loader_file_format=loader_file_format)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            JSONSCHEMA_PIPELINE_NAME,
            "--input-dir",
            str(settings.input_dir),
            "--output-dir",
            str(settings.output_dir),
            "--log-config-file",
            str(settings.log_config_file),
            "--use-destination",
            str(settings.use_destination),
            "--use-output-dir-for-pipeline-metadata",
            "true",
            "--loader-file-format",
            loader_file_format,
            "--dataset-name",
            "schema_cli",
            "-m",
            settings.schema_files_module,
            "--table-names",
            "gadget",
        ],
    )

    load_info = cli.run()

    assert load_info is not None
    assert not load_info.has_failed_jobs
    assert load_info.dataset_name == "schema_cli"
    assert read_data_tables(settings, load_info.dataset_name) == {"gadget": [GADGET_RECORD]}
