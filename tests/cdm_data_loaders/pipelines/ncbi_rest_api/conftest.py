"""Shared fixtures and constants for NCBI REST API pipeline tests."""

from typing import Any

import pytest

from cdm_data_loaders.pipelines import ncbi_rest_api as ncbi_module
from cdm_data_loaders.pipelines.ncbi_rest_api import NcbiRestApiSettings
from tests.conftest import DEFAULT_VCR_CONFIG

ID_WITH_2K_ANNOTS = "GCF_000003135.1"
ID_WITH_500_ANNOTS = "GCF_000007725.1"
ID_TRIGGERS_500_ERR = "GCF_500_ERROR"
VALID_IDS = [ID_WITH_500_ANNOTS, ID_WITH_2K_ANNOTS]
INVALID_ID = "invalid_id"
ALL_IDS = [*VALID_IDS, INVALID_ID]


@pytest.fixture(autouse=True)
def patch_rest_client_hooks(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that the REST_CLIENT_HOOKS dict is empty for tests."""
    monkeypatch.setattr("cdm_data_loaders.pipelines.ncbi_rest_api.REST_CLIENT_HOOKS", {})
    monkeypatch.setattr("cdm_data_loaders.pipelines.ncbi_rest_api.ncbi_genome_client._config.http_show_error_body", True)


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


@pytest.fixture
def test_settings(test_defaults: dict[str, Any]) -> NcbiRestApiSettings:
    """Generate a test settings class."""
    return NcbiRestApiSettings(**test_defaults)


@pytest.fixture(autouse=True)
def reset_settings(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reset the pipeline settings prior to running tests.

    :param monkeypatch: monkeypatch
    :type monkeypatch: pytest.MonkeyPatch
    """
    monkeypatch.setattr(ncbi_module, "PIPELINE_SETTINGS", None)

