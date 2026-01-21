"""Tests for parser error handling, schema compliance, and so on."""

import logging
from typing import Any
from unittest.mock import MagicMock

import pytest
from py4j.protocol import Py4JJavaError
from pyspark.errors import AnalysisException
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import IntegerType, StringType, StructField, StructType

from cdm_data_loader_utils.readers.dsv import read
from tests.conftest import ALL_LINES, MISSING_REQUIRED, TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH, VALID

PERMISSIVE = "PERMISSIVE"
DROP = "DROPMALFORMED"
FAILFAST = "FAILFAST"

INGEST_MODES = [PERMISSIVE, DROP, FAILFAST]

TEST_SCHEMA_FIELD = StructField("what", StringType(), nullable=False)


@pytest.mark.parametrize(
    "schema_fields",
    [
        # valid schema but not the right format
        StructType([TEST_SCHEMA_FIELD, TEST_SCHEMA_FIELD]),
        # not a list
        TEST_SCHEMA_FIELD,
        # list contains raw types
        [StringType(), IntegerType(), TEST_SCHEMA_FIELD],
        # completely wrong!
        {"this": "str"},
    ],
)
def test_read_wrong_schema_format(schema_fields: Any, caplog: pytest.LogCaptureFixture) -> None:  # noqa: ANN401
    """Ensure that the schema is in the right form."""
    err_msg = "schema_fields must be specified as a list of StructFields"
    with pytest.raises(TypeError, match=err_msg):
        read(MagicMock(), "/some/path", schema_fields)

    assert len(caplog.records) == 1
    assert caplog.records[0].levelno == logging.ERROR
    assert caplog.records[0].message == err_msg


@pytest.mark.requires_spark
@pytest.mark.parametrize(("delimiter", "fmt"), [(",", "CSV"), ("\t", "TSV"), ("|", "DSV"), (None, "CSV")])
def test_read_errors(spark: SparkSession, delimiter: str | None, fmt: str, caplog: pytest.LogCaptureFixture) -> None:
    """Check error handling."""
    options = {"delimiter": delimiter} if delimiter else {}

    with pytest.raises(AnalysisException, match="Path does not exist:"):
        read(spark, "/path/to/nowhere", [TEST_SCHEMA_FIELD], options)

    assert len(caplog.records) == 1
    assert caplog.records[0].levelno == logging.ERROR
    assert caplog.records[0].message == f"Failed to load {fmt} from /path/to/nowhere"


@pytest.mark.parametrize("mode", INGEST_MODES)
@pytest.mark.parametrize("csv_lines", [VALID, MISSING_REQUIRED, TYPE_MISMATCH, TOO_FEW_COLS, TOO_MANY_COLS, ALL_LINES])
def test_csv_read_modes(  # noqa: PLR0913
    spark: SparkSession,
    mode: str,
    csv_lines: str,
    csv_schema: list[StructField],
    request: pytest.FixtureRequest,
    caplog: pytest.LogCaptureFixture,
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

    # spark will happily read in the data and report that it loaded the max number of lines
    test_df = read(spark, str(csv_lines_path), csv_schema, options=read_options)
    assert isinstance(test_df, DataFrame)
    # check logging
    assert len(caplog.records) == 1
    assert caplog.records[0].levelno == logging.INFO
    assert (
        caplog.records[0].message
        == f"Loaded {n_rows * 5 if csv_lines == ALL_LINES else n_rows} CSV records from {csv_lines_path!s}"
    )

    if mode == FAILFAST and csv_lines in (TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH, ALL_LINES):
        with pytest.raises(Py4JJavaError, match="An error occurred while calling "):
            test_df.collect()
        return

    read(spark, str(csv_lines_path), csv_schema, options=read_options)

    data_rows = [r.asDict() for r in test_df.collect()]

    # all modes should correctly parse the valid data
    # none of the modes GAF about missing required values, so all will be read
    if csv_lines in (VALID, MISSING_REQUIRED):
        assert len(data_rows) == n_rows
        return

    if csv_lines in (TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH):
        if mode == DROP:
            # dropmalformed will not parse any content from these files as all lines are invalid
            assert len(data_rows) == 0
        else:
            # permissive will pull in all the data
            assert len(data_rows) == n_rows
        return

    # ALL_LINES: permissive will pull in all, DROP will just pull in the VALID + MISSING_REQUIRED lines
    assert len(data_rows) == n_rows * 5 if mode == PERMISSIVE else n_rows * 2
