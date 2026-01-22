"""Tests for parser error handling, schema compliance, and so on."""

from typing import Any

import pytest
from py4j.protocol import Py4JJavaError
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import StringType, StructField

from cdm_data_loader_utils.audit.schema import (
    AUDIT_SCHEMA,
    METRICS,
    N_INVALID,
    N_READ,
    N_VALID,
    VALIDATION_ERRORS,
)
from cdm_data_loader_utils.core.constants import INVALID_DATA_FIELD_NAME
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.readers.dsv import read
from cdm_data_loader_utils.validation.dataframe_validator import DataFrameValidator, Validator
from cdm_data_loader_utils.validation.df_nullable_fields import validate as nullable_fields
from tests.conftest import ALL_LINES, MISSING_REQUIRED, TEST_NS, TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH, VALID
from tests.helpers import create_empty_delta_table

PERMISSIVE = "PERMISSIVE"
DROP = "DROPMALFORMED"
FAILFAST = "FAILFAST"

INGEST_MODES = [PERMISSIVE, DROP, FAILFAST]

TEST_SCHEMA_FIELD = StructField("what", StringType(), nullable=False)


def read_with_validation(
    spark: SparkSession,
    run: PipelineRun,
    path: str,
    schema_fields: list[StructField],
    options: dict[str, Any],
) -> DataFrame:
    """Read in a delimiter-separated file, performing some minimal validation checks.

    :param spark: spark sesh
    :type spark: SparkSession
    :param run: current run info
    :type run: PipelineRun
    :param path: location of the file to parse
    :type path: str
    :param schema_fields: list of StructFields describing the expected input
    :type schema_fields: list[StructField]
    :return: just the valid rows from the dataframe
    :rtype: DataFrame
    """
    parsed_df = read(spark, path, schema_fields, options)

    # validate the output
    validator = DataFrameValidator(spark)
    result = validator.validate_dataframe(
        data_to_validate=parsed_df,
        schema=schema_fields,
        run=run,
        validator=Validator(nullable_fields, {"invalid_col": INVALID_DATA_FIELD_NAME}),
        invalid_col=INVALID_DATA_FIELD_NAME,
    )

    return result.valid_df


@pytest.mark.requires_spark
@pytest.mark.parametrize("mode", INGEST_MODES)
@pytest.mark.parametrize("csv_lines", [VALID, MISSING_REQUIRED, TYPE_MISMATCH, TOO_FEW_COLS, TOO_MANY_COLS, ALL_LINES])
def test_csv_read_with_validation_errors(  # noqa: PLR0913
    spark: SparkSession,
    mode: str,
    csv_lines: str,
    csv_schema: list[StructField],
    request: pytest.FixtureRequest,
    pipeline_run: PipelineRun,
) -> None:
    """Test ingestion of valid and invalid CSV data."""
    n_rows = 4
    csv_lines_path = request.getfixturevalue(csv_lines)

    read_options = {
        "delimiter": ",",
        "header": False,
        "comment": "#",
        "dateFormat": "yyyyMMdd",
        "ignoreLeadingWhiteSpace": True,
        "ignoreTrailingWhiteSpace": True,
        "mode": mode,
    }

    # DROPMALFORMED won't be run
    if mode in (FAILFAST, DROP):
        with pytest.raises(ValueError, match="The only permitted read mode is PERMISSIVE"):
            read_with_validation(spark, pipeline_run, str(csv_lines_path), csv_schema, options=read_options)
        return

    # prep the appropriate tables
    for t in AUDIT_SCHEMA:
        create_empty_delta_table(spark, TEST_NS, t, AUDIT_SCHEMA[t])

    test_df = read_with_validation(spark, pipeline_run, str(csv_lines_path), csv_schema, options=read_options)
    data_rows = [r.asDict() for r in test_df.collect()]
    fields = {f.name for f in csv_schema}
    assert set(test_df.schema.fieldNames()) == fields

    metrics = spark.table(f"{pipeline_run.namespace}.{METRICS}").collect()
    assert len(metrics) == 1
    result = metrics[0]
    assert len(data_rows) == result[N_VALID]

    # all modes should correctly parse the valid data
    if csv_lines == VALID:
        assert result[N_VALID] == n_rows
        assert result[N_READ] == n_rows
        assert result[N_INVALID] == 0

    if csv_lines in (TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH):
        # permissive will pull in all the data
        assert result[N_VALID] == 0
        assert result[N_READ] == n_rows
        assert result[N_INVALID] == n_rows

    # none of the modes GAF about missing required values, so all will be read
    # they will be removed by the validator
    if csv_lines == MISSING_REQUIRED:
        assert result[N_VALID] == 0
        assert result[N_READ] == n_rows
        assert result[N_INVALID] == n_rows
        for n in range(1, 6):
            assert f"missing_required: col{n}" in result[VALIDATION_ERRORS]

    if csv_lines == ALL_LINES:
        assert result[N_VALID] == n_rows
        assert result[N_READ] == n_rows * 5
        assert result[N_INVALID] == n_rows * 4
