from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.pipelines.uniprot_kb import UNIPROT_LOG_INTERVAL, UniProtSettings


@pytest.fixture
def test_settings(test_defaults: dict[str, Any]) -> UniProtSettings:
    """Provide a minimal valid UniProtSettings object, with output_dir resolved against the test config."""
    return UniProtSettings(**test_defaults)


@pytest.fixture
def test_uniprot_settings(test_cts_settings: frozendict[str, Any]) -> frozendict[str, Any]:
    """Return uniprot settings dict extending test_cts_settings."""
    return frozendict({**test_cts_settings, "log_interval": UNIPROT_LOG_INTERVAL})


@pytest.fixture
def test_uniprot_settings_reconciled(test_cts_settings_reconciled: frozendict[str, Any]) -> frozendict[str, Any]:
    """Return reconciled uniprot settings dict extending test_cts_settings_reconciled."""
    return frozendict({**test_cts_settings_reconciled, "log_interval": UNIPROT_LOG_INTERVAL})
