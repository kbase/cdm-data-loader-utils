"""Tests for the nullable fields validator."""

import pytest
from pyspark.sql import SparkSession
from pyspark.sql.types import StructField, StructType

from cdm_data_loaders.audit.schema import ROW_ERRORS
from cdm_data_loaders.core.constants import INVALID_DATA_FIELD_NAME
from cdm_data_loaders.readers.dsv import INVALID_DATA_FIELD
from cdm_data_loaders.validation.df_nullable_fields import validate
from tests.conftest import ALL_LINES, MISSING_REQUIRED, TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH, VALID


@pytest.mark.requires_spark
@pytest.mark.parametrize("csv_lines", [VALID, MISSING_REQUIRED, TYPE_MISMATCH, TOO_FEW_COLS, TOO_MANY_COLS, ALL_LINES])
def test_csv_read_errors(
    spark: SparkSession,
    csv_lines: str,
    request: pytest.FixtureRequest,
    csv_schema: list[StructField],
    invalid_csv_missing_required_annots: list[list[str]],
) -> None:
    """Test ingestion of valid and invalid CSV data."""
    n_rows = 4
    csv_lines_path = request.getfixturevalue(csv_lines)

    read_options = {
        "inferSchema": False,
        "enforceSchema": True,
        "delimiter": ",",
        "header": False,
        "comment": "#",
        "dateFormat": "yyyyMMdd",
        "ignoreLeadingWhiteSpace": True,
        "ignoreTrailingWhiteSpace": True,
        "mode": "PERMISSIVE",
        "columnNameOfCorruptRecord": INVALID_DATA_FIELD_NAME,
    }
    # this mimics what the dsv reader does
    dsv_schema = StructType([*csv_schema, INVALID_DATA_FIELD])
    df = spark.read.options(**read_options).csv(str(csv_lines_path), schema=dsv_schema)

    # any data in invalid_col will show up as a parse_error in annotated_df
    original_df_rows = [r.asDict() for r in df.collect()]
    n_invalid_rows = 0
    for r in original_df_rows:
        if r[INVALID_DATA_FIELD_NAME]:
            n_invalid_rows += 1

    # missing required cols does not show up as an error (thanks, Spark schema validation)
    if csv_lines in (VALID, MISSING_REQUIRED):
        assert n_invalid_rows == 0
    elif csv_lines == ALL_LINES:
        assert n_invalid_rows == n_rows * 3
    else:
        # the other inputs all generate parse errors for each row due to incorrect number of cols or type mismatches
        assert n_invalid_rows == n_rows

    # run it through the validator
    annotated_df = validate(df, csv_schema, INVALID_DATA_FIELD_NAME)
    annotated_rows = [r.asDict() for r in annotated_df.collect()]
    row_errors = [r[ROW_ERRORS] for r in annotated_rows]
    # any data in invalid_col will show up as a parse_error in annotated_df
    assert n_invalid_rows == row_errors.count(["parse_error"])

    # get the specific errors for the missing data
    if csv_lines in (MISSING_REQUIRED, ALL_LINES):
        if csv_lines == MISSING_REQUIRED:
            assert row_errors == invalid_csv_missing_required_annots
        else:
            assert row_errors[4:8] == invalid_csv_missing_required_annots
