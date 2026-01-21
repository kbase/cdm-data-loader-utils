"""Tests for parser error handling, schema compliance, and so on."""

from typing import Any
from unittest.mock import MagicMock

import pytest
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import StructField, StructType

from cdm_data_loader_utils.audit.schema import METRICS, REJECTS, ROW_ERRORS
from cdm_data_loader_utils.core.constants import INVALID_DATA_FIELD_NAME
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.validation.dataframe_validator import DataFrameValidator, Validator
from tests.audit.conftest import create_table


@pytest.mark.requires_spark
def test_validate_dataframe_empty_df(pipeline_run: PipelineRun, empty_df: DataFrame) -> None:
    """Assert that an empty dataframe does not perform any validation."""
    dfv = DataFrameValidator(MagicMock())
    validation_fn = MagicMock()
    assert empty_df.count() == 0

    validator = Validator(validation_fn, {"invalid_col": INVALID_DATA_FIELD_NAME})
    output = dfv.validate_dataframe(
        data_to_validate=empty_df,
        schema=empty_df.schema.fields,
        run=pipeline_run,
        validator=validator,
        invalid_col=INVALID_DATA_FIELD_NAME,
    )
    assert output.records_read == 0
    assert output.records_invalid == 0
    assert output.records_valid == 0
    assert output.validation_errors == []
    assert output.valid_df.count() == 0
    validation_fn.assert_not_called()


@pytest.mark.requires_spark
def test_validate_dataframe_no_validation(  # noqa: PLR0913
    spark: SparkSession,
    csv_schema: list[StructField],
    annotated_df_data: list[dict[str, Any]],
    annotated_df_schema: StructType,
    annotated_df_errors: set[str],
    pipeline_run: PipelineRun,
) -> None:
    """Test the dataframe validator writes to the appropriate tables.

    The validator is mocked for test purposes.
    """
    for t in [METRICS, REJECTS]:
        create_table(spark, t, add_default_data=False)

    dfv = DataFrameValidator(spark)
    annotated_df = spark.createDataFrame(annotated_df_data, schema=annotated_df_schema)
    validation_fn = MagicMock(return_value=annotated_df)

    output = dfv.validate_dataframe(
        data_to_validate=annotated_df,
        schema=csv_schema,
        run=pipeline_run,
        validator=Validator(validation_fn, {}),
        invalid_col=INVALID_DATA_FIELD_NAME,
    )
    n_rows = 4  # 4 rows per csv file
    assert output.records_invalid == n_rows * 4
    assert output.records_valid == n_rows
    assert output.records_read == n_rows * 5

    assert set(output.validation_errors) == annotated_df_errors
    # the first 4 rows are valid; everything else triggers validation errors
    assert [r.asDict() for r in output.valid_df.collect()] == [
        {k: v for k, v in r.items() if k not in (INVALID_DATA_FIELD_NAME, ROW_ERRORS)} for r in annotated_df_data[:4]
    ]

    # check that metrics and rejects have been written
    metrics = spark.table(f"{pipeline_run.namespace}.{METRICS}")
    assert metrics.count() == 1
    rejects = spark.table(f"{pipeline_run.namespace}.{REJECTS}")
    assert rejects.count() == output.records_invalid
