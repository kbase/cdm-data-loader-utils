"""Parametrized tests for the reusable Annotated field types in fields.py."""

from typing import Any
from unittest.mock import MagicMock

import dlt
import pytest
from pydantic import BaseModel, ConfigDict, ValidationError

from cdm_data_loaders.core.fields import (
    BATCH_SIZE,
    BUFFER_SIZE,
    DEFAULTS,
    DEV_MODE,
    INPUT_DIR,
    LOG_CONFIG_FILE,
    LOG_INTERVAL,
    MIN_START_AT,
    OUTPUT_DIR,
    START_AT,
    USE_DESTINATION,
    USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
    VALID_DESTINATIONS,
    BatchSize,
    BufferSize,
    DevMode,
    DltConfig,
    InputDir,
    LogConfigFile,
    LogInterval,
    OutputDir,
    StartAt,
    UseDestination,
    UseOutputDirForPipelineMetadata,
)


class _FieldsModel(BaseModel):
    """Throwaway model exercising each scalar field type in isolation."""

    batch_size: BatchSize
    buffer_size: BufferSize
    dev_mode: DevMode
    input_dir: InputDir
    log_config_file: LogConfigFile
    log_interval: LogInterval
    output_dir: OutputDir
    start_at: StartAt
    use_destination: UseDestination
    use_output_dir_for_pipeline_metadata: UseOutputDirForPipelineMetadata


VALID_CASES: list[tuple[str, Any]] = [
    ("batch_size", 1),
    ("batch_size", 1000),
    ("buffer_size", 1),
    ("buffer_size", 100),
    ("dev_mode", True),
    ("dev_mode", False),
    ("input_dir", "/input_dir"),
    ("input_dir", ""),
    ("log_config_file", None),
    ("log_config_file", "/etc/logging.conf"),
    ("log_interval", 1),
    ("log_interval", 1000),
    ("output_dir", "/output_dir"),
    ("output_dir", ""),
    ("start_at", 1),
    ("start_at", 42),
    ("use_destination", "local_fs"),
    ("use_destination", "s3"),
    ("use_output_dir_for_pipeline_metadata", True),
    ("use_output_dir_for_pipeline_metadata", False),
]

INVALID_CASES: list[Any] = [
    ("batch_size", 0),
    ("batch_size", -5),
    ("batch_size", "not-a-number"),
    ("buffer_size", 0),
    ("buffer_size", -1),
    ("dev_mode", "not-a-bool"),
    ("input_dir", 123),
    ("log_config_file", 123),
    ("log_interval", 0),
    ("log_interval", -1),
    ("output_dir", 123),
    ("start_at", "abc"),
    pytest.param(
        "start_at",
        0,
    ),
    pytest.param(
        "start_at",
        -1,
    ),
    ("use_destination", 123),
    pytest.param(
        "use_destination",
        "not_a_real_destination",
        marks=pytest.mark.xfail(reason="UseDestination is not restricted to VALID_DESTINATIONS", strict=True),
    ),
    ("use_output_dir_for_pipeline_metadata", "not-a-bool"),
]


@pytest.mark.parametrize(("field_name", "value"), VALID_CASES)
def test_valid_field_values(field_name: str, value: Any) -> None:  # noqa: ANN401
    """A known-good value for each field type is accepted and round-trips unchanged."""
    model = _FieldsModel(**{field_name: value})
    assert getattr(model, field_name) == value


@pytest.mark.parametrize(("field_name", "value"), INVALID_CASES)
def test_invalid_field_values(field_name: str, value: Any) -> None:  # noqa: ANN401
    """A known-bad value for each field type is rejected with a ValidationError."""
    with pytest.raises(ValidationError):
        _FieldsModel(**{field_name: value})


DEFAULT_DRIFT_CASES: list[tuple[str, str]] = [
    ("batch_size", BATCH_SIZE),
    ("buffer_size", BUFFER_SIZE),
    ("dev_mode", DEV_MODE),
    ("input_dir", INPUT_DIR),
    ("log_config_file", LOG_CONFIG_FILE),
    ("log_interval", LOG_INTERVAL),
    ("output_dir", OUTPUT_DIR),
    ("start_at", START_AT),
    ("use_destination", USE_DESTINATION),
    ("use_output_dir_for_pipeline_metadata", USE_OUTPUT_DIR_FOR_PIPELINE_METADATA),
]


@pytest.mark.parametrize(("attr_name", "defaults_key"), DEFAULT_DRIFT_CASES)
def test_field_default_matches_defaults_mapping(attr_name: str, defaults_key: str) -> None:
    """Each Annotated field's default value matches the corresponding entry in DEFAULTS."""
    model = _FieldsModel()  # pyright: ignore[reportCallIssue]
    assert getattr(model, attr_name) == DEFAULTS[defaults_key]


def test_output_dir_default_is_placeholder_empty_string() -> None:
    """OUTPUT_DIR's default is the empty-string placeholder meant to be overridden downstream."""
    assert DEFAULTS[OUTPUT_DIR] == ""


# Numeric constraint boundaries (PositiveInt fields)
POSITIVE_INT_FIELDS = ["batch_size", "buffer_size", "log_interval"]


@pytest.mark.parametrize("field_name", POSITIVE_INT_FIELDS)
def test_positive_int_field_rejects_zero(field_name: str) -> None:
    """Fields typed as PositiveInt reject zero."""
    with pytest.raises(ValidationError):
        _FieldsModel(**{field_name: 0})  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("field_name", POSITIVE_INT_FIELDS)
def test_positive_int_field_rejects_negative(field_name: str) -> None:
    """Fields typed as PositiveInt reject negative values."""
    with pytest.raises(ValidationError):
        _FieldsModel(**{field_name: -1})  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize("field_name", POSITIVE_INT_FIELDS)
def test_positive_int_field_accepts_one(field_name: str) -> None:
    """Fields typed as PositiveInt accept the smallest valid value, 1."""
    model = _FieldsModel(**{field_name: 1})  # pyright: ignore[reportArgumentType]
    assert getattr(model, field_name) == 1


# DltConfig: test to ensure that dlt's config doesn't change type
class _DltConfigModel(BaseModel):
    """DltConfig is arbitrary-typed and container-based, so test it separately."""

    model_config = ConfigDict(arbitrary_types_allowed=True)

    dlt_config: DltConfig


def test_dlt_config_accepts_dict() -> None:
    """DltConfig accepts a plain dict, the ease-of-testing fallback type."""
    model = _DltConfigModel(dlt_config={"key": "value"})
    assert model.dlt_config == {"key": "value"}


def test_dlt_config_accepts_real_dlt_accessor() -> None:
    """DltConfig accepts the actual dlt.config accessor (dlt's private _Accessor type)."""
    model = _DltConfigModel(dlt_config=dlt.config)
    assert model.dlt_config is dlt.config


def test_dlt_config_default_factory_produces_dlt_config(monkeypatch: pytest.MonkeyPatch) -> None:
    """DltConfig's default_factory produces dlt.config when no value is supplied."""
    dlt_conf = MagicMock()
    monkeypatch.setattr(dlt, "config", dlt_conf)

    class _WithDefault(BaseModel):
        model_config = ConfigDict(arbitrary_types_allowed=True)
        dlt_config: DltConfig

    assert _WithDefault().dlt_config is dlt_conf  # pyright: ignore[reportCallIssue]


def test_dlt_config_model_dump_does_not_raise() -> None:
    """A model containing DltConfig can be dumped without raising, guarding against a dlt upgrade breaking _Accessor."""
    model = _DltConfigModel(dlt_config=dlt.config)
    model.model_dump()


def test_dlt_config_rejects_non_dict_non_accessor() -> None:
    """DltConfig rejects a value that is neither a dict nor a dlt _Accessor."""
    with pytest.raises(ValidationError):
        _DltConfigModel(dlt_config=123)


# Module-level constants: immutability and consistency and other stuff that
# really doesn't need to be tested, but whatever
def test_defaults_mapping_is_immutable() -> None:
    """DEFAULTS is a frozendict and rejects item assignment."""
    with pytest.raises(TypeError):
        DEFAULTS[BATCH_SIZE] = 5  # type: ignore[index]


def test_valid_destinations_contains_expected_values() -> None:
    """VALID_DESTINATIONS currently contains exactly 'local_fs' and 's3'."""
    assert VALID_DESTINATIONS == ["local_fs", "s3"]


def test_valid_destinations_is_a_plain_mutable_list() -> None:
    """VALID_DESTINATIONS is (still) a plain, unprotected list rather than a tuple or Final."""
    assert isinstance(VALID_DESTINATIONS, list)


def test_min_start_at_matches_defaults_start_at() -> None:
    """MIN_START_AT and the START_AT default in DEFAULTS agree with each other."""
    assert DEFAULTS[START_AT] == MIN_START_AT
