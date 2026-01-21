"""Generic DSV file reader with validation of incoming data."""

from typing import Any

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import StringType, StructField, StructType

from cdm_data_loader_utils.core.constants import INVALID_DATA_FIELD_NAME
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger, log_and_die

# mapping of delimiters to format names (for logging)
# spark defaults to separating on commas if nothing is specified
FORMAT_NAME = {None: "CSV", ",": "CSV", "\t": "TSV"}

# schema field to catch input errors
INVALID_DATA_FIELD = StructField(INVALID_DATA_FIELD_NAME, StringType(), nullable=True)

# read modes
PERMISSIVE = "PERMISSIVE"
FAILFAST = "FAILFAST"
DROP = "DROPMALFORMED"

# default options: enforce schema, catch parse errors in INVALID_DATA_FIELD_NAME column
DEFAULT_DSV_OPTIONS = {
    "inferSchema": False,
    "enforceSchema": True,
    "mode": PERMISSIVE,
    "columnNameOfCorruptRecord": INVALID_DATA_FIELD_NAME,
}


logger = get_cdm_logger()


def get_format_name(delimiter: str | None) -> str:
    """Get the nice name of the format being parsed, given the string used as delimiter."""
    return FORMAT_NAME.get(delimiter, "DSV")


def read(
    spark: SparkSession,
    path: str,
    schema_fields: list[StructField],
    options: dict[str, Any] | None = None,
) -> DataFrame:
    """Read in a delimiter-separated file with spark.

    :param spark: spark sesh
    :type spark: SparkSession
    :param path: location of the file to parse
    :type path: str
    :param schema_fields: list of StructFields describing the expected input
    :type schema_fields: list[StructField]
    :param options: dictionary of options
    :type options: dict[str, Any]
    :return: dataframe of parsed rows
    :rtype: DataFrame
    """
    if not isinstance(schema_fields, list) or not all(isinstance(field, StructField) for field in schema_fields):
        log_and_die("schema_fields must be specified as a list of StructFields", TypeError)

    if not options:
        options = {}

    dsv_schema = StructType([*schema_fields, INVALID_DATA_FIELD])
    dsv_options = {
        **DEFAULT_DSV_OPTIONS,
        **options,
    }

    # if dsv_options["mode"] not in (PERMISSIVE, FAILFAST):
    #     msg = "The only permitted read modes are PERMISSIVE and FAILFAST."
    #     log_and_die(msg, ValueError)

    format_name = get_format_name(options.get("delimiter", options.get("sep")))

    try:
        df = spark.read.options(**dsv_options).csv(path, schema=dsv_schema)
    except Exception:
        # Log the full stack trace and re-raise to be handled by the caller
        logger.exception("Failed to load %s from %s", format_name, path)
        raise

    # count will not trigger an error even if in FAILFAST mode and all records are corrupt.
    logger.info("Loaded %d %s records from %s", df.count(), format_name, path)
    return df


def read_tsv(
    spark: SparkSession, path: str, schema_fields: list[StructField], options: dict[str, Any] | None = None
) -> DataFrame:
    """Shortcut for reading in a tab-separated file.

    :param spark: spark sesh
    :type spark: SparkSession
    :param path: location of the file to parse
    :type path: str
    :param schema_fields: list of StructFields describing the expected input
    :type schema_fields: list[StructField]
    :param options: dictionary of options
    :type options: dict[str, Any]
    :return: dataframe of parsed rows
    :rtype: DataFrame
    """
    if not options:
        options = {}
    options["separator"] = "\t"
    return read(spark, path, schema_fields, options)


def read_csv(
    spark: SparkSession, path: str, schema_fields: list[StructField], options: dict[str, Any] | None = None
) -> DataFrame:
    """Shortcut for reading in a comma-separated file.

    :param spark: spark sesh
    :type spark: SparkSession
    :param path: location of the file to parse
    :type path: str
    :param schema_fields: list of StructFields describing the expected input
    :type schema_fields: list[StructField]
    :param options: dictionary of options
    :type options: dict[str, Any]
    :return: dataframe of parsed rows
    :rtype: DataFrame
    """
    if not options:
        options = {}
    options["separator"] = ","
    return read(spark, path, schema_fields, options)
