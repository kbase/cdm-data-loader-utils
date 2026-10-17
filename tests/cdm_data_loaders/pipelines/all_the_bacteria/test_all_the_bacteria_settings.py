"""Tests for the all_the_bacteria pipeline."""

from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.pipelines.all_the_bacteria import (
    AtbSettings,
)
from tests.helpers import assert_cli_field_roundtrips, assert_no_cli_clashes


# model_fields
def test_no_cli_param_clashes() -> None:
    """Ensure that the settings object does not have any internal alias collisions."""
    assert_no_cli_clashes(AtbSettings)


@pytest.mark.parametrize("field_name", [f for f in AtbSettings.model_fields if f != "pattern_file"])
def test_cli_fields_parse_correctly(
    field_name: str,
    atb_settings_pair: tuple[frozendict[str, Any], frozendict[str, Any]],
    test_defaults: dict[str, Any],
) -> None:
    """Ensure that all fields and variants are parsed correctly."""
    settings, reconciled = atb_settings_pair
    assert_cli_field_roundtrips(
        AtbSettings,
        field_name,
        raw_value=settings[field_name],
        expected_value=reconciled[field_name],
        required=test_defaults,
    )


def test_cli_pattern_file_parses_correctly(
    atb_settings_pair: tuple[frozendict[str, Any], frozendict[str, Any]],
    test_defaults: dict[str, Any],
) -> None:
    """Ensure that pattern_file is parsed correctly when set."""
    settings, reconciled = atb_settings_pair
    if settings["pattern_file"] is None:
        return
    assert_cli_field_roundtrips(
        AtbSettings,
        "pattern_file",
        raw_value=settings["pattern_file"],
        expected_value=reconciled["pattern_file"],
        required=test_defaults,
    )
