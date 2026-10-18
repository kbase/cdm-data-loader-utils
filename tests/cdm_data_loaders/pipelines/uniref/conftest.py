from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.pipelines.uniref import (
    UNIREF_VARIANTS,
    VARIANT,
    UnirefSettings,
)

TEST_DEFAULT_UNIREF_VARIANT = "50"

UNIREF_VARIANT_ALIASES = ["v", "variant"]


@pytest.fixture
def test_uniref_settings(test_cts_settings: frozendict[str, Any]) -> frozendict[str, Any]:
    """Return uniref settings dict extending test_cts_settings."""
    return frozendict({**test_cts_settings, VARIANT: TEST_DEFAULT_UNIREF_VARIANT})


@pytest.fixture
def test_uniref_settings_reconciled(test_cts_settings_reconciled: frozendict[str, Any]) -> frozendict[str, Any]:
    """Return reconciled uniref settings dict extending test_cts_settings_reconciled."""
    return frozendict({**test_cts_settings_reconciled, VARIANT: TEST_DEFAULT_UNIREF_VARIANT})


@pytest.fixture(params=UNIREF_VARIANTS)
def uniref_variant_value(request: pytest.FixtureRequest) -> str:
    """Parametrized fixture over all valid uniref variants."""
    return request.param


@pytest.fixture
def test_settings(uniref_variant_value: str, test_defaults: dict[str, Any]) -> UnirefSettings:
    """A valid UnirefSettings object for each uniref variant, with output_dir resolved against the test config."""
    return UnirefSettings(variant=uniref_variant_value, **test_defaults)  # type: ignore[reportReturnType]
