"""Class providing an interface to validation with error and metrics auditing."""

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as sf
from pyspark.sql.types import StructField

from cdm_data_loader_utils.audit.metrics import write_metrics
from cdm_data_loader_utils.audit.rejects import write_rejects
from cdm_data_loader_utils.audit.schema import ROW_ERRORS
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger
from cdm_data_loader_utils.validation.validation_result import ValidationResult

logger = get_cdm_logger()


@dataclass
class Validator:
    """Base validator dataclass."""

    validation_fn: Callable
    args: dict[str, Any]


class DataFrameValidator:
    """Class for validating data."""

    def __init__(self, spark: SparkSession) -> None:
        """Instantiate an IngestValidator.

        :param spark: spark sesh
        :type spark: SparkSession
        """
        self.spark = spark

    def validate_dataframe(
        self,
        data_to_validate: DataFrame,
        schema: list[StructField],
        run: PipelineRun,
        validator: Validator,
        invalid_col: str,
    ) -> ValidationResult:
        """Validate a dataframe, outputting a ValidationResult.

        :param data_to_validate: dataframe to be validated
        :type data_to_validate: DataFrame
        :param schema: schema of the fields in the dataframe, expressed as a list of StructFields
        :type schema: list[StructField]
        :param run: the current pipeline run
        :type run: PipelineRun
        :param validation_fn: function for validating the dataframe
        :type validation_fn: Callable
        :return: data validation metrics and the valid data from the input dataframe
        :rtype: ValidationResult
        """
        if data_to_validate.count() == 0:
            logger.warning("%s %s: no data found to validate. Aborting.")
            return ValidationResult(
                valid_df=data_to_validate,
                records_read=0,
                records_valid=0,
                records_invalid=0,
                validation_errors=[],
            )

        # running the validator should produce a dataframe with a column called ROW_ERRORS
        annotated_df: DataFrame = validator.validation_fn(data_to_validate, schema, **validator.args)
        valid_df: DataFrame = annotated_df.filter(sf.size(ROW_ERRORS) == 0)
        write_rejects(
            run=run,
            annotated_df=annotated_df,
            schema_fields=schema,
            invalid_col=invalid_col,
        )
        metrics = write_metrics(self.spark, annotated_df, run)

        return ValidationResult(
            valid_df=valid_df.drop(ROW_ERRORS).drop(invalid_col),
            records_read=metrics.records_read,
            records_valid=metrics.records_valid,
            records_invalid=metrics.records_invalid,
            validation_errors=metrics.validation_errors,
        )
