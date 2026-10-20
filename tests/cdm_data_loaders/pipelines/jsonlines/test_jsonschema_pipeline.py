"""Tests for JSON Schema validation and rejected-record routing."""

import gzip
import json
import sys
from collections.abc import Callable
from pathlib import Path
from types import ModuleType
from typing import Any
from uuid import uuid4

import pytest
from pydantic import ValidationError

from cdm_data_loaders.pipelines.jsonlines.base import make_page_router
from cdm_data_loaders.pipelines.jsonlines.extract_jsonschema_validate_pipeline import (
    _make_resolve_error,
    build_entity_resource,
    load_entity_schemas,
    run_jsonlines_pipeline,
)
from cdm_data_loaders.pipelines.jsonlines.settings import JsonlJsonschemaIngestSettings


def make_item(record: object, line_no: int = 1) -> dict[str, Any]:
    """Build the shared reader's parsed-line envelope."""
    return {
        "record": record,
        "raw_record": json.dumps(record, separators=(",", ":")),
        "parse_error": None,
        "source_file": "data.jsonl",
        "line_no": line_no,
    }


def make_router(schema: dict[str, Any], buffer_size: int = 2) -> Callable[[list[dict[str, Any]]], Any]:
    """Build the router the jsonschema pipeline uses for the widget table."""
    return make_page_router("widget", buffer_size, _make_resolve_error(schema))


def test_make_page_router_pass_parse_error_preserves_raw_record() -> None:
    """Malformed JSON bypasses validation and preserves its source text."""
    validate = make_router({"$schema": "https://json-schema.org/draft/2020-12/schema"})
    item = {
        "record": None,
        "raw_record": "{bad json",
        "parse_error": "bad json",
        "source_file": "data.jsonl",
        "line_no": 3,
    }

    results = list(validate([item]))

    assert [(result.meta.table_name, result.data) for result in results] == [
        (
            "widget_rejected",
            [{"source_file": "data.jsonl", "line_no": 3, "raw_record": "{bad json", "error_detail": "bad json"}],
        )
    ]


def test_make_page_router_fail_uses_declared_draft() -> None:
    """Draft 4 boolean exclusiveMinimum rejects the boundary value."""
    schema: dict[str, Any] = {
        "$schema": "http://json-schema.org/draft-04/schema#",
        "type": "object",
        "properties": {"count": {"type": "number", "minimum": 5, "exclusiveMinimum": True}},
    }
    validate = make_router(schema)
    item = {
        "record": {"count": 5},
        "raw_record": '{"count": 5}',
        "parse_error": None,
        "source_file": "data.jsonl",
        "line_no": 1,
    }

    results = list(validate([item]))

    assert [result.meta.table_name for result in results] == ["widget_rejected"]


def test_load_entity_schemas_pass_returns_registry(
    schema_module_factory: Callable[[object], str], widget_schema: dict[str, Any]
) -> None:
    """Valid schemas from multiple supported drafts are returned unchanged."""
    registry = {
        "widget": widget_schema,
        "older": {"$schema": "http://json-schema.org/draft-04/schema#", "type": "object"},
    }
    assert load_entity_schemas(schema_module_factory(registry)) is registry


@pytest.mark.parametrize(
    ("registry", "message"),
    [
        pytest.param(None, "non-empty ENTITY_SCHEMAS", id="none"),
        pytest.param({}, "non-empty ENTITY_SCHEMAS", id="empty"),
        pytest.param([], "non-empty ENTITY_SCHEMAS", id="list"),
        pytest.param({"": {}}, "non-empty table names", id="empty-table"),
        pytest.param({1: {}}, "non-empty table names", id="non-string-table"),
        pytest.param({"widget": {}}, r"missing the \$schema keyword", id="missing-draft"),
        pytest.param({"widget": {"$schema": "https://example.invalid/schema"}}, "Unsupported", id="unknown-draft"),
    ],
)
def test_load_entity_schemas_fail_malformed_registry(
    schema_module_factory: Callable[[object], str], registry: object, message: str
) -> None:
    """Malformed registries and unsupported drafts fail before extraction."""
    with pytest.raises(ValueError, match=message):
        load_entity_schemas(schema_module_factory(registry))


@pytest.mark.parametrize(
    "schema",
    [pytest.param(True, id="boolean-schema"), pytest.param({"$schema": []}, id="non-string-draft")],
)
def test_load_entity_schemas_fail_invalid_types(schema_module_factory: Callable[[object], str], schema: object) -> None:
    """Schema dictionaries and draft declarations must have the documented types."""
    with pytest.raises(TypeError, match="JSON Schema for widget must"):
        load_entity_schemas(schema_module_factory({"widget": schema}))


def test_load_entity_schemas_fail_missing_registry(monkeypatch: pytest.MonkeyPatch) -> None:
    """A module must explicitly provide ENTITY_SCHEMAS."""
    module_name = f"_test_missing_registry_{uuid4().hex}"
    monkeypatch.setitem(sys.modules, module_name, ModuleType(module_name))
    with pytest.raises(ValueError, match="non-empty ENTITY_SCHEMAS"):
        load_entity_schemas(module_name)


def test_load_entity_schemas_fail_missing_module() -> None:
    """Import failures retain their original exception."""
    with pytest.raises(ModuleNotFoundError):
        load_entity_schemas(f"_missing_schema_module_{uuid4().hex}")


def test_load_entity_schemas_fail_reports_all_invalid_schemas(
    schema_module_factory: Callable[[object], str], widget_schema: dict[str, Any]
) -> None:
    """Schema errors identify every failing table, not just the first."""
    module_name = schema_module_factory(
        {"first": {**widget_schema, "type": "invalid"}, "second": {**widget_schema, "required": 7}}
    )
    with pytest.raises(RuntimeError, match="errors were found") as error:
        load_entity_schemas(module_name)
    assert all(text in str(error.value) for text in [module_name, "first:", "second:"])


def test_make_page_router_pass_preserves_records_and_empty_pages(widget_schema: dict[str, Any]) -> None:
    """Valid data retains extra fields, dates, and nested values without inserted defaults."""
    record = {"widget_id": "a", "count": 0, "collected": "2024-02-29", "nested": {"items": [1, 2]}}
    validate = make_router(widget_schema)
    results = list(validate([make_item(record)]))
    assert [(result.meta.table_name, result.data) for result in results] == [("widget", [record])]
    assert results[0].data[0] is record
    assert "label" not in record
    assert list(validate([])) == []


@pytest.mark.parametrize(
    ("record", "messages"),
    [
        pytest.param({"count": 1}, ["'widget_id' is a required property"], id="missing-required"),
        pytest.param({"widget_id": "a", "count": "1"}, ["'1' is not of type 'integer'"], id="no-coercion"),
        pytest.param({"widget_id": "a", "count": True}, ["True is not of type 'integer'"], id="boolean-not-integer"),
        pytest.param({"widget_id": "a", "count": -1}, ["-1 is less than the minimum of 0"], id="minimum"),
        pytest.param(
            {"widget_id": "a", "count": 1, "collected": "2023-02-29"},
            ["'2023-02-29' is not a 'date'"],
            id="invalid-format",
        ),
        pytest.param({}, ["'count' is a required property", "'widget_id' is a required property"], id="all-errors"),
        pytest.param(None, ["None is not of type 'object'"], id="null"),
        pytest.param([], ["[] is not of type 'object'"], id="array"),
        pytest.param(42, ["42 is not of type 'object'"], id="scalar"),
    ],
)
def test_make_page_router_fail_reports_exact_errors(
    widget_schema: dict[str, Any], record: object, messages: list[str]
) -> None:
    """Invalid records retain their raw text and report every validation error."""
    item = make_item(record, line_no=7)
    results = list(make_router(widget_schema)([item]))
    assert [(result.meta.table_name, result.data) for result in results] == [
        (
            "widget_rejected",
            [
                {
                    "source_file": "data.jsonl",
                    "line_no": 7,
                    "raw_record": item["raw_record"],
                    "error_detail": json.dumps(messages),
                }
            ],
        )
    ]


def test_make_page_router_pass_flushes_independent_buffers(widget_schema: dict[str, Any]) -> None:
    """Full and partial valid/rejected pages flush without leaking between calls."""
    valid = {"widget_id": "a", "count": 0}
    items = [make_item(record, index) for index, record in enumerate([valid, {}, valid, {}, valid, {}], start=1)]
    validate = make_router(widget_schema)
    for _iteration in range(2):
        pages = list(validate(items))
        assert [(page.meta.table_name, len(page.data)) for page in pages] == [
            ("widget", 2),
            ("widget_rejected", 2),
            ("widget", 1),
            ("widget_rejected", 1),
        ]
        assert [row["line_no"] for page in pages if page.meta.table_name == "widget_rejected" for row in page.data] == [
            2,
            4,
            6,
        ]


@pytest.mark.parametrize("compressed", [False, True], ids=["plain", "gzip"])
def test_build_entity_resource_pass_reads_selected_files(
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings],
    widget_schema: dict[str, Any],
    compressed: bool,
) -> None:
    """The real resource reads plain/gzip files, skips blanks, and honors the glob."""
    settings = schema_settings_factory()
    entity_dir = Path(settings.input_dir) / "widget"
    entity_dir.mkdir()
    record = {"widget_id": "a", "count": 0}
    contents = f"\n{json.dumps(record)}\n\n{{broken\n"
    filename = "data.jsonl.gz" if compressed else "data.jsonl"
    if compressed:
        with gzip.open(entity_dir / filename, "wt", encoding="utf-8") as stream:
            stream.write(contents)
    else:
        (entity_dir / filename).write_text(contents, encoding="utf-8")
    (entity_dir / "ignored.txt").write_text("not JSON", encoding="utf-8")

    resource = build_entity_resource("widget", widget_schema, settings)
    rows = list(resource)

    assert record in rows
    rejected = next(row for row in rows if "error_detail" in row)
    expected_error = "Expecting property name enclosed in double quotes: line 1 column 2 (char 1)"
    assert rejected == {"source_file": filename, "line_no": 4, "raw_record": "{broken", "error_detail": expected_error}
    assert rows == [record, rejected]


@pytest.mark.parametrize("create_directory", [False, True], ids=["missing", "empty"])
def test_build_entity_resource_pass_no_files(
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings],
    widget_schema: dict[str, Any],
    create_directory: bool,
) -> None:
    """Missing and empty entity directories yield no records."""
    settings = schema_settings_factory()
    if create_directory:
        (Path(settings.input_dir) / "widget").mkdir()
    assert list(build_entity_resource("widget", widget_schema, settings)) == []


def test_run_jsonlines_pipeline_fail_unknown_tables(
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings],
) -> None:
    """Unknown tables fail before the shared runner creates pipeline state."""
    settings = schema_settings_factory(table_names=["unknown"])
    with pytest.raises(ValueError, match=r"Unknown table name\(s\) requested: \['unknown'\]"):
        run_jsonlines_pipeline(settings)
    assert list(Path(settings.output_dir).iterdir()) == []


@pytest.mark.parametrize(
    ("table_names", "expected"),
    [
        pytest.param(None, None, id="all-tables"),
        pytest.param("widget, gadget", ["widget", "gadget"], id="comma-separated"),
        pytest.param(["widget"], ["widget"], id="list"),
        pytest.param(" , ", None, id="blank"),
    ],
)
def test_settings_pass_table_selection(
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings],
    table_names: str | list[str] | None,
    expected: list[str] | None,
) -> None:
    """JSON Schema settings retain the shared table selection contract."""
    settings = schema_settings_factory(table_names=table_names)
    assert settings.table_names == expected
    assert settings.model_config["cli_prog_name"] == "jsonlines_jsonschema_validator"


@pytest.mark.parametrize("buffer_size", [0, -1], ids=["zero", "negative"])
def test_settings_fail_invalid_buffer_size(
    schema_settings_factory: Callable[..., JsonlJsonschemaIngestSettings], buffer_size: int
) -> None:
    """Settings reject nonpositive buffers before resource construction."""
    with pytest.raises(ValidationError, match="buffer_size"):
        schema_settings_factory(buffer_size=buffer_size)
