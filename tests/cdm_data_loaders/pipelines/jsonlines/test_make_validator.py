"""Tests for make_page_router: the valid/invalid routing logic."""

import json
from collections.abc import Callable
from typing import Any

import pytest
from dlt.extract.items import DataItemWithMeta, TableNameMeta
from pydantic import ValidationError

from cdm_data_loaders.pipelines.jsonlines.base import make_page_router


def make_item(
    record: dict[str, Any] | None = None,
    raw_record: str = "",
    parse_error: str | None = None,
    source_file: str = "file.jsonl",
    line_no: int = 1,
) -> dict[str, Any]:
    """Build a line record matching stream_jsonl_lines' output shape."""
    return {
        "record": record,
        "raw_record": raw_record,
        "parse_error": parse_error,
        "source_file": source_file,
        "line_no": line_no,
    }


def run_validator(validator: Callable[[list[dict[str, Any]]], Any], item: dict[str, Any]) -> list[DataItemWithMeta]:
    """Run the validator on a single-item page. Return the tagged results."""
    return list(validator([item]))


def make_router(model: Any, **kwargs: Any) -> Callable[[list[dict[str, Any]]], Any]:
    """Build a router for the widget table using a model, mirroring the pydantic pipeline."""

    def _resolve_error(record: dict[str, Any]) -> str | None:
        try:
            model.model_validate(record)
        except ValidationError as error:
            return json.dumps(error.errors(), default=str)
        return None

    return make_page_router(
        "widget",
        kwargs.get("buffer_size", 10),
        _resolve_error,
        **{key: value for key, value in kwargs.items() if key != "buffer_size"},
    )


def test_make_page_router_pass_valid_record_routed_to_main_table(widget_model: Any) -> None:
    """A record that passes validation goes to the main table, dumped from the validated model."""
    validator = make_router(widget_model)
    item = make_item(record={"widget_id": "a", "count": 5}, raw_record='{"widget_id": "a", "count": 5}')

    results = run_validator(validator, item)

    assert len(results) == 1
    result = results[0]
    assert isinstance(result.meta, TableNameMeta)
    assert result.meta.table_name == "widget"
    assert result.data == [{"widget_id": "a", "count": 5}]


def test_make_page_router_fail_json_parse_error_routed_to_rejected_without_model_validation(widget_model: Any) -> None:
    """A record with a parse_error goes to <table>_rejected. The model is not called."""
    validator = make_router(widget_model, include_record=True)
    item = make_item(record=None, raw_record="{not json", parse_error="Expecting value: line 1 column 1 (char 0)")

    results = run_validator(validator, item)

    assert len(results) == 1
    result = results[0]
    assert result.meta.table_name == "widget_rejected"
    assert result.data == [
        {
            "error_detail": "Expecting value: line 1 column 1 (char 0)",
            "line_no": 1,
            "raw_record": "{not json",
            "record": None,
            "source_file": "file.jsonl",
        }
    ]


@pytest.mark.parametrize(
    "record",
    [
        pytest.param({"count": 5}, id="missing_required_widget_id"),
        pytest.param({"widget_id": "a", "count": -1}, id="count_below_minimum_boundary"),
        pytest.param({"widget_id": "a", "count": "not-a-number"}, id="count_not_coercible"),
        pytest.param({"widget_id": "", "count": 1}, id="widget_id_below_min_length_boundary"),
    ],
)
def test_make_page_router_fail_pydantic_validation_error_routed_to_rejected(
    widget_model: Any, record: dict[str, Any]
) -> None:
    """A record that parses as JSON but fails model validation goes to <table>_rejected, with error detail."""
    validator = make_router(widget_model)
    item = make_item(record=record, raw_record=json.dumps(record))

    results = run_validator(validator, item)

    assert len(results) == 1
    result = results[0]
    assert result.meta.table_name == "widget_rejected"
    rejected_rows = result.data
    assert len(rejected_rows) == 1
    assert rejected_rows[0]["raw_record"] == json.dumps(record)
    error_detail = json.loads(rejected_rows[0]["error_detail"])
    assert isinstance(error_detail, list)
    assert len(error_detail) >= 1
    assert all({"loc", "msg", "type"} <= set(error.keys()) for error in error_detail)


def test_make_page_router_pass_mixed_batch_routes_each_record_independently(widget_model: Any) -> None:
    """A mixed batch routes each record on its own. Valid and invalid records do not affect each other."""
    validator = make_router(widget_model)
    items = [
        make_item(record={"widget_id": "a", "count": 1}, raw_record='{"widget_id": "a", "count": 1}', line_no=1),
        make_item(record={"widget_id": "b", "count": -1}, raw_record='{"widget_id": "b", "count": -1}', line_no=2),
        make_item(record={"widget_id": "c", "count": 3}, raw_record='{"widget_id": "c", "count": 3}', line_no=3),
        make_item(record=None, raw_record="{bad", parse_error="boom", line_no=4),
    ]

    results = list(validator(items))

    valid_rows = [r for r in results if r.meta.table_name == "widget"]
    assert len(valid_rows) == 1
    assert [row["widget_id"] for row in valid_rows[0].data] == ["a", "c"]

    rejected_rows = [r for r in results if r.meta.table_name == "widget_rejected"]
    assert len(rejected_rows) == 1
    assert {row["line_no"] for row in rejected_rows[0].data} == {2, 4}


def test_make_page_router_fail_all_invalid_none_reach_main_table(widget_model: Any) -> None:
    """When every record in a batch is invalid, none reach the main table."""
    validator = make_router(widget_model)
    items = [
        make_item(record={"count": 1}, raw_record='{"count": 1}', line_no=1),
        make_item(record={"widget_id": "b", "count": -5}, raw_record='{"widget_id": "b", "count": -5}', line_no=2),
        make_item(record=None, raw_record="not json at all", parse_error="bad json", line_no=3),
    ]

    results = list(validator(items))

    assert len(results) == 1
    assert results[0].meta.table_name == "widget_rejected"
    assert len(results[0].data) == 3


def test_make_page_router_pass_buffer_size_controls_batching(widget_model: Any) -> None:
    """buffer_size controls how many rows land in each yielded page per table."""
    validator = make_router(widget_model, buffer_size=2)
    items = [
        make_item(record={"widget_id": f"w{i}", "count": 1}, raw_record=f'{{"count": {i}}}', line_no=i + 1)
        for i in range(5)
    ]

    pages = tagged_pages(list(validator(items)))

    assert pages == [("widget", 2), ("widget", 2), ("widget", 1)]


def tagged_pages(results: list[DataItemWithMeta]) -> list[tuple[str, int]]:
    """Return (table_name, row_count) pairs for tagged pages."""
    return [(r.meta.table_name, len(r.data)) for r in results]
