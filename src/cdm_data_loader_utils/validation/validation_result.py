"""Dataclass for capturing validation results."""

from dataclasses import dataclass

from pyspark.sql import DataFrame


@dataclass(frozen=True)
class ValidationResult:
    """Dataclass for capturing validation results."""

    valid_df: DataFrame
    records_read: int
    records_valid: int
    records_invalid: int
    validation_errors: list[str]
