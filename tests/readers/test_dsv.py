"""Tests for parser error handling, schema compliance, and so on."""

import logging
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pytest
from pyspark.errors import AnalysisException
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import IntegerType, StringType, StructField, StructType
from pyspark.testing import assertSchemaEqual

from cdm_data_loader_utils.readers.dsv import read, read_csv, read_tsv
from tests.conftest import ALL_LINES, MISSING_REQUIRED, TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH, VALID
from tests.helpers import assertDataFrameEqual

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


@pytest.mark.requires_spark
def test_read_tsv_csv(spark: SparkSession, csv_schema: list[StructField], all_lines: Path, all_lines_tsv: Path) -> None:
    """Ensure that the TSV and CSV reader shortcuts work as expected."""
    read_options = {
        "header": False,
        "comment": "#",
        "dateFormat": "yyyyMMdd",
        "ignoreLeadingWhiteSpace": False,
        "ignoreTrailingWhiteSpace": False,
    }
    csv_options = {"delimiter": ","}
    tsv_options = {"delimiter": "\t"}

    test_df_csv = read(spark, str(all_lines), csv_schema, options={**read_options, **csv_options})
    test_df_tsv = read(spark, str(all_lines_tsv), csv_schema, options={**read_options, **tsv_options})

    csv_df = read_csv(spark, str(all_lines), csv_schema, read_options)
    tsv_df = read_tsv(spark, str(all_lines_tsv), csv_schema, read_options)

    for df in [test_df_tsv, csv_df, tsv_df]:
        assertSchemaEqual(test_df_csv.schema, df.schema)

    # TODO: compare the TSV and CSV versions?
    assertDataFrameEqual(test_df_tsv.collect(), tsv_df.collect())
    assertDataFrameEqual(test_df_csv.collect(), csv_df.collect())


@pytest.mark.requires_spark
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

    if mode in (DROP, FAILFAST):
        with pytest.raises(ValueError, match="The only permitted read mode is PERMISSIVE"):
            read(spark, str(csv_lines_path), csv_schema, options=read_options)
        return

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

    read(spark, str(csv_lines_path), csv_schema, options=read_options)

    data_rows = [r.asDict() for r in test_df.collect()]

    # all modes should correctly parse the valid data
    # none of the modes GAF about missing required values, so all will be read
    if csv_lines in (VALID, MISSING_REQUIRED):
        assert len(data_rows) == n_rows
        return

    if csv_lines in (TOO_FEW_COLS, TOO_MANY_COLS, TYPE_MISMATCH):
        # permissive will pull in all the data
        assert len(data_rows) == n_rows
        return

    # ALL_LINES: permissive will pull in all, DROP will just pull in the VALID + MISSING_REQUIRED lines
    assert len(data_rows) == n_rows * 5
