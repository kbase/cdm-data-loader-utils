"""Tests for the constrained field types in cdm_data_loaders.core.fields that the core settings classes do not use.

Fields used by the classes in cdm_data_loaders.core.settings are exercised via those classes in test_settings.py.
BatchSize is the same PositiveInt construct as BufferSize/LogInterval. FileGlob, PreserveTableNesting and TableName
are unconstrained str/bool fields, and the module constants are plain values, so none of those are tested here.
"""

import pytest
from pydantic import BaseModel, ValidationError

from cdm_data_loaders.core.fields import DatasetName, LoaderFileFormat, LoaderFileFormatEnum, MaxTableNesting

FORMAT_ERROR = "Input should be 'jsonl' or 'parquet'"


class PipelineFieldsModel(BaseModel):
    """Model holding the constrained field types that pipeline-specific settings classes add."""

    dataset_name: DatasetName
    loader_file_format: LoaderFileFormat
    max_table_nesting: MaxTableNesting


@pytest.mark.parametrize(
    ("kwargs", "expected"),
    [
        pytest.param(
            {"dataset_name": "ds"},
            {"dataset_name": "ds", "loader_file_format": "parquet", "max_table_nesting": 0},
            id="defaults",
        ),
        pytest.param(
            {"dataset_name": "ds", "loader_file_format": "jsonl", "max_table_nesting": 3},
            {"dataset_name": "ds", "loader_file_format": "jsonl", "max_table_nesting": 3},
            id="jsonl_string-nesting_3",
        ),
        pytest.param(
            {"dataset_name": "ds", "loader_file_format": LoaderFileFormatEnum.JSONL, "max_table_nesting": 0},
            {"dataset_name": "ds", "loader_file_format": "jsonl", "max_table_nesting": 0},
            id="jsonl_enum_member-nesting_0",
        ),
    ],
)
def test_pipeline_fields_pass_valid_values(kwargs: dict[str, object], expected: dict[str, object]) -> None:
    """Valid values are accepted; the format is an enum member whose str() is the plain name passed to dlt."""
    model = PipelineFieldsModel.model_validate(kwargs)
    assert model.model_dump() == expected
    assert isinstance(model.loader_file_format, LoaderFileFormatEnum)
    assert str(model.loader_file_format) == expected["loader_file_format"]


@pytest.mark.parametrize(
    ("kwargs", "expected_errors"),
    [
        pytest.param(
            {"dataset_name": "ds", "max_table_nesting": -1},
            [(("max_table_nesting",), "Input should be greater than or equal to 0")],
            id="negative_nesting",
        ),
        pytest.param(
            {"dataset_name": "ds", "loader_file_format": "csv"},
            [(("loader_file_format",), FORMAT_ERROR)],
            id="unknown_format",
        ),
        pytest.param(
            {"dataset_name": "ds", "loader_file_format": "PARQUET"},
            [(("loader_file_format",), FORMAT_ERROR)],
            id="format_wrong_case",
        ),
        pytest.param(
            {"dataset_name": ""},
            [(("dataset_name",), "String should have at least 1 character")],
            id="empty_dataset_name",
        ),
        pytest.param({}, [(("dataset_name",), "Field required")], id="missing_dataset_name"),
    ],
)
def test_pipeline_fields_fail_invalid_values(
    kwargs: dict[str, object], expected_errors: list[tuple[tuple[str, ...], str]]
) -> None:
    """Out-of-range nesting, unknown formats, and empty or missing dataset names are rejected."""
    with pytest.raises(ValidationError) as exc_info:
        PipelineFieldsModel.model_validate(kwargs)
    assert [(err["loc"], err["msg"]) for err in exc_info.value.errors()] == expected_errors
