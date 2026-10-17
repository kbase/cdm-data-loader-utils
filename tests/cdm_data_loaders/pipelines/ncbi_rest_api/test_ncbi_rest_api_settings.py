"""Tests for the NCBI datasets API pipeline functions."""

import re
from typing import Any

import pytest
from frozendict import frozendict
from pydantic import ValidationError
from pydantic_settings import CliApp

from cdm_data_loaders.pipelines import ncbi_rest_api as ncbi_module
from cdm_data_loaders.pipelines.ncbi_rest_api import (
    ANNOTATION,
    DATASET,
    MAX_IDS_PER_QUERY,
    NcbiRestApiSettings,
    get_settings,
    set_settings,
)
from tests.cdm_data_loaders.core.conftest import (
    TEST_CTS_SETTINGS,
    TEST_CTS_SETTINGS_RECONCILED,
    check_settings,
)
from tests.conftest import DEFAULT_VCR_CONFIG
from tests.helpers import assert_cli_field_roundtrips, assert_no_cli_clashes


@pytest.fixture(autouse=True)
def patch_rest_client_hooks(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that the REST_CLIENT_HOOKS dict is empty for tests."""
    monkeypatch.setattr("cdm_data_loaders.pipelines.ncbi_rest_api.REST_CLIENT_HOOKS", {})


ID_WITH_2K_ANNOTS = "GCF_000003135.1"
ID_WITH_500_ANNOTS = "GCF_000007725.1"
ID_TRIGGERS_500_ERR = "GCF_500_ERROR"
VALID_IDS = [ID_WITH_500_ANNOTS, ID_WITH_2K_ANNOTS]
INVALID_ID = "invalid_id"
ALL_IDS = [*VALID_IDS, INVALID_ID]
BATCH_SIZE = 200
BATCH_SIZE_STRING = str(BATCH_SIZE)

TEST_NCBI_SETTINGS = frozendict(**TEST_CTS_SETTINGS, batch_size=BATCH_SIZE_STRING, query_type="")

TEST_NCBI_SETTINGS_RECONCILED = frozendict(
    **TEST_CTS_SETTINGS_RECONCILED, batch_size=BATCH_SIZE, query_type=None, save_raw_responses=False
)

TEST_NCBI_SETTINGS_V1 = frozendict(**TEST_CTS_SETTINGS, batch_size=BATCH_SIZE_STRING, query_type=" ANNOTATION ")

TEST_NCBI_SETTINGS_RECONCILED_V1 = frozendict(
    **TEST_CTS_SETTINGS_RECONCILED, batch_size=BATCH_SIZE, query_type=ANNOTATION, save_raw_responses=False
)


@pytest.fixture(scope="module")
def vcr_config() -> dict[str, Any]:
    """VCR config for tests that make HTTP requests."""
    return {**DEFAULT_VCR_CONFIG}


@pytest.fixture(scope="module")
def valid_assembly_ids() -> list[str]:
    """A list of assembly IDs."""
    return VALID_IDS


@pytest.fixture(scope="module")
def assembly_id(valid_assembly_ids: list[str]) -> str:
    """Single valid assembly ID."""
    return valid_assembly_ids[0]


@pytest.fixture(scope="module")
def invalid_assembly_id() -> str:
    """Invalid assembly ID."""
    return INVALID_ID


@pytest.fixture(scope="module")
def assembly_ids(valid_assembly_ids: list[str], invalid_assembly_id: str) -> list[str]:
    """List of assembly IDs including both valid and invalid IDs."""
    return [*valid_assembly_ids, invalid_assembly_id]


@pytest.fixture(scope="module")
def test_settings() -> NcbiRestApiSettings:
    """Generate a test settings class."""
    return NcbiRestApiSettings()  # type: ignore[reportReturnType]


@pytest.fixture(autouse=True)
def reset_settings(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reset the pipeline settings prior to running tests.

    :param monkeypatch: monkeypatch
    :type monkeypatch: pytest.MonkeyPatch
    """
    monkeypatch.setattr(ncbi_module, "PIPELINE_SETTINGS", None)


def test_get_set_settings(test_settings: NcbiRestApiSettings) -> None:
    """Test setting and retrieval of settings."""
    assert ncbi_module.PIPELINE_SETTINGS is None
    with pytest.raises(RuntimeError, match="Pipeline settings have not been initialised"):
        get_settings()

    set_settings(test_settings)
    assert test_settings == ncbi_module.PIPELINE_SETTINGS
    assert get_settings() == test_settings

    set_settings(None)  # type: ignore[reportArgumentType]
    assert ncbi_module.PIPELINE_SETTINGS is None
    with pytest.raises(RuntimeError, match="Pipeline settings have not been initialised"):
        get_settings()


@pytest.mark.parametrize(
    ("batch_size", "parsed_batch_size"),
    [
        ("10", 10),
        (10, 10),
        (None, MAX_IDS_PER_QUERY),
    ],
)
def test_settings_valid_batch_size(
    batch_size: int | str | None, parsed_batch_size: int, test_settings: NcbiRestApiSettings
) -> None:
    """Ensure that a valid batch size is correctly parsed."""
    if batch_size is None:
        settings = test_settings
    else:
        settings: NcbiRestApiSettings = NcbiRestApiSettings(batch_size=batch_size)  # type: ignore[reportReturnType]
    assert settings.batch_size == parsed_batch_size


@pytest.mark.parametrize(
    ("bad_batch_size", "message"),
    [
        ("0", "Input should be greater than or equal to 1"),
        ("-1", "Input should be greater than or equal to 1"),
        ("1001", f"Input should be less than or equal to {MAX_IDS_PER_QUERY}"),
        ("notanint", "Input should be a valid integer"),
        ("", "Input should be a valid integer"),
        ("1.2345", "Input should be a valid integer"),
    ],
)
@pytest.mark.parametrize("use_cliapp", [True, False])
def test_cli_invalid_batch_size_via_cli_raises(bad_batch_size: str, message: str, use_cliapp: bool) -> None:
    """Ensure that an invalid batch size passed via CLI raises an error."""
    if use_cliapp:
        with pytest.raises(ValidationError, match=message):
            CliApp.run(NcbiRestApiSettings, cli_args=["--batch-size", bad_batch_size])
    else:
        with pytest.raises(ValidationError, match=message):
            NcbiRestApiSettings(batch_size=bad_batch_size)  # pyright: ignore[reportCallIssue]


@pytest.mark.parametrize(
    ("query_type", "parsed_query_type"),
    [
        (ANNOTATION, ANNOTATION),
        (DATASET, DATASET),
        ("  ANNOTATION  ", ANNOTATION),
        ("\n\nDataset\t\n", DATASET),
        ("   ", None),
        ("", None),
        (None, None),
    ],
)
def test_cli_valid_query_type(
    query_type: str | None, parsed_query_type: str | None, test_settings: NcbiRestApiSettings
) -> None:
    """Ensure that a valid batch size is correctly parsed."""
    settings: NcbiRestApiSettings = NcbiRestApiSettings(query_type=query_type)  # type: ignore[reportReturnType]
    assert settings.query_type == parsed_query_type

    if query_type is None:
        alt_settings: NcbiRestApiSettings = test_settings
        assert alt_settings.query_type == parsed_query_type


STRING_MATCH_MESSAGE = re.escape("String should match pattern '^(dataset|annotation)$'")


@pytest.mark.parametrize(
    "query_type",
    [
        "annot",
        "anotation",
        "data set",
    ],
)
@pytest.mark.parametrize("use_cliapp", [True, False])
def test_cli_invalid_query_type(query_type: str, use_cliapp: bool) -> None:
    """Ensure that an invalid batch size passed via CLI raises an error."""
    if use_cliapp:
        with pytest.raises(ValidationError, match=STRING_MATCH_MESSAGE):
            CliApp.run(NcbiRestApiSettings, cli_args=["--query-type", query_type])
    else:
        with pytest.raises(ValidationError, match=STRING_MATCH_MESSAGE):
            NcbiRestApiSettings(query_type=query_type)  # pyright: ignore[reportCallIssue]


@pytest.mark.parametrize(
    ("settings", "reconciled"),
    [
        (TEST_NCBI_SETTINGS, TEST_NCBI_SETTINGS_RECONCILED),
        (TEST_NCBI_SETTINGS_V1, TEST_NCBI_SETTINGS_RECONCILED_V1),
    ],
)
def test_settings_all_params_set(settings: frozendict, reconciled: frozendict) -> None:
    """Ensure that settings are set correctly when all args are specified."""
    s = NcbiRestApiSettings(**settings)
    check_settings(s, reconciled)


# model_fields
def test_no_cli_param_clashes() -> None:
    """Ensure that the settings object does not have any internal alias collisions."""
    assert_no_cli_clashes(NcbiRestApiSettings)


@pytest.mark.parametrize(
    ("settings", "reconciled"),
    [
        (TEST_NCBI_SETTINGS, TEST_NCBI_SETTINGS_RECONCILED),
        (TEST_NCBI_SETTINGS_V1, TEST_NCBI_SETTINGS_RECONCILED_V1),
    ],
)
@pytest.mark.parametrize("field_name", list(TEST_NCBI_SETTINGS))
def test_cli_fields_parse_correctly(field_name: str, settings: frozendict, reconciled: frozendict) -> None:
    """Ensure that all fields and variants are parsed correctly."""
    assert_cli_field_roundtrips(
        NcbiRestApiSettings,
        field_name,
        settings[field_name],
        reconciled[field_name],
    )
