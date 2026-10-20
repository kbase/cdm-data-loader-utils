"""Tests for the XSV reader (post-qsv row casting and parsing-config resolution)."""

from pathlib import Path
from typing import Any, Final

import pytest

from cdm_data_loaders.readers.jsonschema_xsv.xsv_reader import (
    cast_xsv_row,
    cast_xsv_value,
    read_validated_rows,
    resolve_xsv_parsing_config,
)
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import ValidatedSchema

VALID_SCHEMA_URI: Final[str] = "https://json-schema.org/draft/2020-12/schema"
COLUMNS: Final[list[str]] = ["number", "date", "float", "boolean", "string"]


def _schema(properties: dict[str, Any] | None = None, x_xsv_config: dict[str, Any] | None = None) -> ValidatedSchema:
    """Build a minimal ValidatedSchema, optionally with properties and/or an x-xsv-config block."""
    jsonschema: dict[str, Any] = {"$schema": VALID_SCHEMA_URI, "required": COLUMNS}
    if properties is not None:
        jsonschema["properties"] = properties
    if x_xsv_config is not None:
        jsonschema["x-xsv-config"] = x_xsv_config
    return ValidatedSchema(jsonschema=jsonschema)


"""cast_xsv_value"""


@pytest.mark.parametrize(
    ("value", "type_spec", "expected"),
    [
        pytest.param("5", "integer", 5, id="integer"),
        pytest.param("-3", ["integer"], -3, id="integer-in-list"),
        pytest.param("3.14", "number", 3.14, id="number"),
        pytest.param("3.14", ["number", "null"], 3.14, id="number-nullable-non-empty"),
        pytest.param("true", "boolean", True, id="boolean-true"),
        pytest.param("True", "boolean", True, id="boolean-true-mixed-case"),
        pytest.param("1", "boolean", True, id="boolean-true-numeric"),
        pytest.param("false", "boolean", False, id="boolean-false"),
        pytest.param("no", "boolean", False, id="boolean-false-unrecognised-token"),
        pytest.param("hello", "string", "hello", id="string-unchanged"),
        pytest.param("hello", None, "hello", id="no-type-non-empty-passthrough"),
        pytest.param("5", ["string", "integer"], 5, id="first-type-without-caster-skipped"),
        pytest.param("", ["string", "null"], None, id="empty-nullable-becomes-none"),
        pytest.param("", None, None, id="empty-no-type-becomes-none"),
        pytest.param("", "string", "", id="empty-non-nullable-string-stays-empty"),
    ],
)
def test_cast_xsv_value_pass_conversions(value: str, type_spec: str | list[str] | None, expected: object) -> None:
    """cast_xsv_value converts a raw string to the type(s) declared for it."""
    assert cast_xsv_value(value, type_spec) == expected


def test_cast_xsv_row_pass_casts_every_column() -> None:
    """cast_xsv_row casts each field independently using its own property's declared type."""
    properties = {
        "number": {"type": "integer"},
        "date": {"type": "string"},
        "float": {"type": ["string", "null"]},
        "boolean": {"type": "boolean"},
        "string": {"type": "string"},
    }
    row = {"number": "2", "date": "2023-01-15", "float": "", "boolean": "true", "string": "key:value1"}

    assert cast_xsv_row(row, properties) == {
        "number": 2,
        "date": "2023-01-15",
        "float": None,
        "boolean": True,
        "string": "key:value1",
    }


def test_cast_xsv_row_pass_unknown_column_defaults_to_no_type() -> None:
    """A column absent from properties is treated as having no declared type."""
    assert cast_xsv_row({"unknown": ""}, {}) == {"unknown": None}
    assert cast_xsv_row({"unknown": "text"}, {}) == {"unknown": "text"}


"""resolve_xsv_parsing_config"""


def test_resolve_xsv_parsing_config_pass_no_config_block_returns_empty_dict() -> None:
    """A schema without an x-xsv-config block resolves to an empty dict."""
    assert resolve_xsv_parsing_config(_schema()) == {}


def test_resolve_xsv_parsing_config_pass_maps_all_keys() -> None:
    """Every supported x-xsv-config key maps to its CleanerValidatorArgs kwarg name."""
    schema = _schema(
        x_xsv_config={
            "x-delimiter": ",",
            "x-comment-char": ";",
            "x-quote": "'",
            "x-escape": "\\",
            "x-null-regex": "^NA$",
            "x-null-cols": ["float"],
            "x-has-header": False,
        }
    )

    assert resolve_xsv_parsing_config(schema) == {
        "delimiter": ",",
        "comment_char": ";",
        "quote": "'",
        "escape": "\\",
        "null_regex": "^NA$",
        "null_regex_cols": ["float"],
        "missing_header": True,
    }


def test_resolve_xsv_parsing_config_pass_has_header_true_inverts_to_missing_header_false() -> None:
    """x-has-header: true resolves to missing_header: false."""
    schema = _schema(x_xsv_config={"x-has-header": True})
    assert resolve_xsv_parsing_config(schema) == {"missing_header": False}


def test_resolve_xsv_parsing_config_pass_partial_config_only_sets_present_keys() -> None:
    """Only keys present in x-xsv-config appear in the resolved config."""
    schema = _schema(x_xsv_config={"x-delimiter": "\t"})
    assert resolve_xsv_parsing_config(schema) == {"delimiter": "\t"}


def test_resolve_xsv_parsing_config_pass_quoting_policy_not_propagated() -> None:
    """x-quoting-policy configures XSV output, not qsv cleaning, and is never returned."""
    schema = _schema(x_xsv_config={"x-delimiter": ",", "x-quoting-policy": "necessary"})
    config = resolve_xsv_parsing_config(schema)
    assert "quoting_policy" not in config
    assert config == {"delimiter": ","}


"""read_validated_rows"""


def _write_xsv(tmp_path: Path, content: str, file_name: str = "data.csv") -> Path:
    path = tmp_path / file_name
    path.write_text(content)
    return path


def test_read_validated_rows_pass_casts_rows_by_property_type(tmp_path: Path) -> None:
    """Rows are cast according to the declared property types, keyed by header column name."""
    properties = {
        "number": {"type": "integer"},
        "date": {"type": "string"},
        "float": {"type": ["string", "null"]},
        "boolean": {"type": "boolean"},
        "string": {"type": "string"},
    }
    content = "number,date,float,boolean,string\n2,2023-01-15,3.14,true,key:value1\n3,2023-02-20,,false,key:value2\n"
    xsv_path = _write_xsv(tmp_path, content)

    rows = list(read_validated_rows(xsv_path, properties, delimiter=","))

    assert rows == [
        {"number": 2, "date": "2023-01-15", "float": "3.14", "boolean": True, "string": "key:value1"},
        {"number": 3, "date": "2023-02-20", "float": None, "boolean": False, "string": "key:value2"},
    ]


def test_read_validated_rows_pass_tab_delimited(tmp_path: Path) -> None:
    r"""A tab-delimited file is read correctly when delimiter="\t"."""
    properties = {"a": {"type": "integer"}, "b": {"type": "string"}}
    xsv_path = _write_xsv(tmp_path, "a\tb\n1\tfoo\n2\tbar\n", file_name="data.tsv")

    rows = list(read_validated_rows(xsv_path, properties, delimiter="\t"))

    assert rows == [{"a": 1, "b": "foo"}, {"a": 2, "b": "bar"}]


def test_read_validated_rows_pass_custom_quote_char(tmp_path: Path) -> None:
    """A field containing the delimiter is read as one field when wrapped in the custom quote char."""
    properties = {"a": {"type": "string"}, "b": {"type": "string"}}
    xsv_path = _write_xsv(tmp_path, "a,b\n'contains, a comma',plain\n")

    rows = list(read_validated_rows(xsv_path, properties, delimiter=",", quote="'"))

    assert rows == [{"a": "contains, a comma", "b": "plain"}]


def test_read_validated_rows_pass_custom_escape_char_disables_doublequote(tmp_path: Path) -> None:
    """A backslash-escaped quote inside a field is read literally when escape is set."""
    properties = {"a": {"type": "string"}}
    xsv_path = _write_xsv(tmp_path, 'a\n"contains \\"a quote\\""\n')

    rows = list(read_validated_rows(xsv_path, properties, delimiter=",", escape="\\"))

    assert rows == [{"a": 'contains "a quote"'}]


def test_read_validated_rows_pass_header_only_file_yields_no_rows(tmp_path: Path) -> None:
    """A file with only a header row yields no rows."""
    xsv_path = _write_xsv(tmp_path, "number,date,float,boolean,string\n")

    assert list(read_validated_rows(xsv_path, {}, delimiter=",")) == []
