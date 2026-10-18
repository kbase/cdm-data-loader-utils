"""Tests for the UniRef DLT pipeline."""

from typing import Any

import pytest
from frozendict import frozendict
from pydantic import ValidationError
from pydantic_settings import CliApp

from cdm_data_loaders.pipelines.uniref import (
    UNIREF_VARIANTS,
    VARIANT,
    UnirefSettings,
)
from tests.cdm_data_loaders.conftest import check_settings
from tests.cdm_data_loaders.pipelines.uniref.conftest import (
    UNIREF_VARIANT_ALIASES,
)
from tests.helpers import assert_cli_field_roundtrips, assert_no_cli_clashes, make_cli_arg


def test_uniref_settings_all_params_set(
    test_uniref_settings: frozendict[str, Any],
    test_uniref_settings_reconciled: frozendict[str, Any],
) -> None:
    """Ensure that settings are set correctly when all args are specified.

    Note that test_uniref_settings includes a value for pipeline_dir.
    """
    s = UnirefSettings(**test_uniref_settings)  # type: ignore[reportReturnType]
    check_settings(s, test_uniref_settings_reconciled)


# model_fields
def test_no_cli_param_clashes() -> None:
    """Ensure that the settings object does not have any internal alias collisions."""
    assert_no_cli_clashes(UnirefSettings)


@pytest.mark.parametrize("field_name", list(UnirefSettings.model_fields))
def test_cli_fields_parse_correctly(
    field_name: str,
    test_uniref_settings: frozendict[str, Any],
    test_uniref_settings_reconciled: frozendict[str, Any],
    test_defaults: dict[str, Any],
) -> None:
    """Ensure that all fields and variants are parsed correctly."""
    assert_cli_field_roundtrips(
        UnirefSettings,
        field_name,
        test_uniref_settings[field_name],
        test_uniref_settings_reconciled[field_name],
        {**test_defaults, VARIANT: test_uniref_settings[VARIANT]},
    )


@pytest.mark.parametrize("uniref_variant_value", UNIREF_VARIANTS)
def test_settings_valid_variants_accepted(uniref_variant_value: str, test_defaults: dict[str, Any]) -> None:
    """Ensure that each valid variant value is accepted without error."""
    s = UnirefSettings(VARIANT=uniref_variant_value, **test_defaults)  # pyright: ignore[reportCallIssue]
    assert isinstance(s, UnirefSettings)
    assert s.variant == uniref_variant_value


@pytest.mark.parametrize("uniref_variant_value", UNIREF_VARIANTS)
@pytest.mark.parametrize(VARIANT, UNIREF_VARIANT_ALIASES)
def test_cli_valid_variants_accepted(uniref_variant_value: str, variant: str, test_defaults: dict[str, Any]) -> None:
    """Ensure that each valid uniref variant value is accepted without error when passed via CLI."""
    default_args = [arg for k, v in test_defaults.items() for arg in (make_cli_arg(k), str(v))]
    s = CliApp.run(
        UnirefSettings,
        cli_args=[*default_args, make_cli_arg(variant), uniref_variant_value],
    )
    assert isinstance(s, UnirefSettings)
    assert s.variant == uniref_variant_value


@pytest.mark.parametrize("value", ["25", "75", "uniref50", "", "ALL"])
def test_invalid_variant_raises(value: str, test_defaults: dict[str, Any]) -> None:
    """Ensure that an unrecognised uniref variant raises a ValidationError."""
    with pytest.raises(ValidationError, match="1 validation error for UnirefSettings") as exc_info:
        UnirefSettings(VARIANT=value, **test_defaults)  # pyright: ignore[reportCallIssue]

    exc_message = str(exc_info.value)
    assert "Input should be '50', '90' or '100'" in exc_message


@pytest.mark.parametrize("value", ["25", "75", "uniref50", "", "ALL"])
@pytest.mark.parametrize(VARIANT, UNIREF_VARIANT_ALIASES)
def test_cli_invalid_variant_via_cli_raises(
    value: str,
    variant: str,
    test_defaults: dict[str, Any],
) -> None:
    """Ensure that an invalid uniref variant passed via CLI raises an error."""
    default_args = [arg for k, v in test_defaults.items() for arg in (make_cli_arg(k), str(v))]
    with pytest.raises(ValidationError, match="1 validation error for UnirefSettings") as exc_info:
        CliApp.run(UnirefSettings, cli_args=[*default_args, make_cli_arg(variant), value])

    exc_message = str(exc_info.value)
    assert "Input should be '50', '90' or '100'" in exc_message


def test_missing_required_uniref_variant_raises(test_defaults: dict[str, Any]) -> None:
    """Ensure that omitting the required uniref variant argument raises a ValidationError."""
    with pytest.raises(ValidationError, match="1 validation error for UnirefSettings") as exc_info:
        UnirefSettings(**test_defaults)  # pyright: ignore[reportCallIssue]

    exc_message = str(exc_info.value)
    assert "Field required" in exc_message


def test_cli_missing_required_uniref_variant_raises(test_defaults: dict[str, Any]) -> None:
    """Ensure that omitting the required uniref variant argument raises an error."""
    default_args = [arg for k, v in test_defaults.items() for arg in (make_cli_arg(k), str(v))]
    with pytest.raises(ValidationError, match="1 validation error for UnirefSettings") as exc_info:
        CliApp.run(UnirefSettings, cli_args=default_args)

    exc_message = str(exc_info.value)
    assert "Field required" in exc_message


@pytest.mark.parametrize("value", ["25", "75", "uniref50", "", "ALL"])
@pytest.mark.parametrize(VARIANT, UNIREF_VARIANT_ALIASES)
def test_cli_invalid_variant_and_destination_via_cli_raises(
    value: str, variant: str, test_defaults: dict[str, Any]
) -> None:
    """Ensure that invalid uniref variant and use_destination passed via CLI raises an error with both errors.

    N.b. Pydantic only reports the first error!!!
    """
    default_args = [
        arg for k, v in test_defaults.items() if k != "use_destination" for arg in (make_cli_arg(k), str(v))
    ]
    with pytest.raises(ValidationError, match="1 validation error for UnirefSettings") as exc_info:
        CliApp.run(
            UnirefSettings,
            cli_args=[
                *default_args,
                make_cli_arg(variant),
                value,
                "--use-destination",
                "some invalid destination",
            ],
        )

    # Check that both errors are present in the exception message
    exc_message = str(exc_info.value)
    assert "Input should be '50', '90' or '100'" in exc_message
