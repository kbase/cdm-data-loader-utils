"""Simple validator for checking for null columns in a dataframe."""

from pyspark.sql import DataFrame
from pyspark.sql import functions as sf
from pyspark.sql.types import StructField

from cdm_data_loader_utils.audit.schema import ROW_ERRORS

COLLECTED_ERRORS = "collected_errors"


def validate(
    df: DataFrame,
    schema_fields: list[StructField],
    invalid_col: str,
) -> DataFrame:
    """Validation function that ensures that nullability constraints in the schema are enforced.

    As of Jan 2026, Spark automagically converts the `nullable=False` attribute on StructFields to
    True when a schema is applied to a dataframe. This function should be supplied with the original
    schema, which it will use to check that null constraints are actually enforced!

    :param df: dataframe to test
    :type df: DataFrame
    :param schema_fields: schema, in the form of a list of SchemaFields
    :type schema_fields: list[StructField]
    :param invalid_col: name of the column where invalid data is stored
    :type invalid_col: str
    :return: df containing rejected data
    :rtype: DataFrame
    """
    # if a column is null but the schema says that it should not be, add the notation "missing_required:{col.name}"
    missing_fields = [
        sf.when(
            sf.col(col.name).isNull(),
            sf.lit(f"missing_required: {col.name}"),
        )
        for col in schema_fields
        # required columns only
        if not col.nullable
    ]

    return (
        df.withColumn(COLLECTED_ERRORS, sf.array(*missing_fields))
        .withColumn(
            ROW_ERRORS,
            sf.when(
                # from spark data ingestion: if the incoming row does not match the supplied schema,
                # incorrect data will go in invalid_col
                sf.col(invalid_col).isNotNull(),
                # mark this as a parse error
                sf.array(sf.lit("parse_error")),
            ).otherwise(
                sf.filter(
                    sf.col(COLLECTED_ERRORS),
                    lambda x: x.isNotNull(),
                )
            ),
        )
        .drop(COLLECTED_ERRORS)
    )
