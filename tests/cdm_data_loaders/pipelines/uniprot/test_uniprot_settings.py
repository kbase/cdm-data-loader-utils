"""Tests for the UniProt DLT pipeline."""

from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.pipelines.uniprot_kb import UniProtSettings
from tests.cdm_data_loaders.conftest import check_settings
from tests.helpers import assert_cli_field_roundtrips, assert_no_cli_clashes


def test_uniprot_settings_all_params_set(
    test_uniprot_settings: frozendict[str, Any],
    test_uniprot_settings_reconciled: frozendict[str, Any],
) -> None:
    """Ensure that settings are set correctly when all args are specified.

    Note that test_uniprot_settings includes a value for pipeline_dir.
    """
    s = UniProtSettings(**test_uniprot_settings)  # type: ignore[reportReturnType]
    check_settings(s, test_uniprot_settings_reconciled)


# model_fields
def test_no_cli_param_clashes() -> None:
    """Ensure that the settings object does not have any internal alias collisions."""
    assert_no_cli_clashes(UniProtSettings)


@pytest.mark.parametrize("field_name", list(UniProtSettings.model_fields))
def test_cli_fields_parse_correctly(
    field_name: str,
    test_uniprot_settings: frozendict[str, Any],
    test_uniprot_settings_reconciled: frozendict[str, Any],
    test_defaults: dict[str, Any],
) -> None:
    """Ensure that all fields and variants are parsed correctly."""
    assert_cli_field_roundtrips(
        UniProtSettings,
        field_name,
        test_uniprot_settings[field_name],
        test_uniprot_settings_reconciled[field_name],
        required=test_defaults,
    )
