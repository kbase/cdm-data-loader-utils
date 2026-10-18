"""Tests for cdm_data_loaders.core.settings."""

from collections.abc import Mapping
from typing import cast

import pytest
from pydantic import AliasChoices, AliasPath, Field
from pydantic_settings import SettingsConfigDict, SettingsError

from cdm_data_loaders.core.settings import CdmDataLoadersBase, _alias_names


def define_settings(
    fields: Mapping[str, object],
    cli_shortcuts: Mapping[str, str | list[str]] | None = None,
    *,
    cli_kebab_case: bool = True,
) -> type[CdmDataLoadersBase]:
    """Define a str-field CdmDataLoadersBase subclass named 'MySettingsClass'; class creation runs check_aliases."""
    config = SettingsConfigDict(cli_kebab_case=cli_kebab_case)
    if cli_shortcuts is not None:
        config["cli_shortcuts"] = dict(cli_shortcuts)
    namespace = {"__annotations__": dict.fromkeys(fields, str), **fields, "model_config": config}
    return cast("type[CdmDataLoadersBase]", type("MySettingsClass", (CdmDataLoadersBase,), namespace))


# _alias_names
@pytest.mark.parametrize(
    ("alias", "expected"),
    [
        pytest.param(None, set(), id="none"),
        pytest.param("alpha", {"alpha"}, id="str"),
        pytest.param(AliasPath("outer", "inner", 0), {"outer"}, id="alias_path_str_first"),
        pytest.param(AliasPath(0), set(), id="alias_path_int_first"),  # pyright: ignore[reportArgumentType]
        pytest.param(AliasChoices("a", "b"), {"a", "b"}, id="alias_choices_str"),
        pytest.param(
            AliasChoices("a", AliasPath("b", 0), AliasPath(1, "c")),  # pyright: ignore[reportArgumentType]
            {"a", "b"},
            id="alias_choices_mixed",
        ),
    ],
)
def test_alias_names_pass(alias: str | AliasChoices | AliasPath | None, expected: set[str]) -> None:
    """Each alias form yields its top-level string names only."""
    assert _alias_names(alias) == expected


# CdmDataLoadersBase.check_aliases
@pytest.mark.parametrize(
    ("fields", "cli_shortcuts", "cli_kebab_case"),
    [
        pytest.param({}, None, True, id="no_fields"),
        pytest.param({"alpha": "a", "beta": "b"}, {"alpha": "a"}, True, id="str_shortcut"),
        pytest.param({"alpha": "a"}, {"alpha": ["a", "x"]}, True, id="list_shortcuts"),
        pytest.param(
            {"alpha": Field("a", validation_alias=AliasChoices("first", AliasPath("second", 0))), "beta": "b"},
            {"second": "s"},
            True,
            id="shortcut_targets_alias_path_name",
        ),
        pytest.param(
            {"alpha_beta": Field("x", validation_alias="alpha-beta")},
            {"alpha_beta": "a"},
            False,
            id="no_kebab-own_alias_differs_only_by_dash",
        ),
    ],
)
def test_check_aliases_pass_valid_names(
    fields: dict[str, object], cli_shortcuts: dict[str, str | list[str]] | None, cli_kebab_case: bool
) -> None:
    """Non-colliding names, aliases and shortcuts are accepted when the class is defined."""
    settings_cls = define_settings(fields, cli_shortcuts, cli_kebab_case=cli_kebab_case)
    assert settings_cls.check_aliases() is None


@pytest.mark.parametrize(
    ("fields", "cli_shortcuts", "expected"),
    [
        pytest.param(
            {"help": "x"},
            None,
            "MySettingsClass: CLI name 'help' for 'help' is reserved by argparse",
            id="reserved_field_name",
        ),
        pytest.param(
            {"alpha": "a"},
            {"alpha": "h"},
            "MySettingsClass: CLI name 'h' for 'alpha' is reserved by argparse",
            id="reserved_shortcut",
        ),
        pytest.param(
            {"alpha": "a", "beta": "b"},
            {"alpha": "beta"},
            "MySettingsClass: CLI name 'beta' is claimed by both 'alpha' and 'beta'",
            id="shortcut_is_other_field_name",
        ),
        pytest.param(
            {"alpha": "a", "beta": "b"},
            {"alpha": "x", "beta": "x"},
            "MySettingsClass: CLI name 'x' is claimed by both 'alpha' and 'beta'",
            id="same_shortcut_two_fields",
        ),
        pytest.param(
            {"log_config_file": "a", "other": "b"},
            {"other": "log_config_file"},
            "MySettingsClass: CLI name 'log_config_file' is claimed by both 'log_config_file' and 'other'",
            id="underscore_dash_collision",
        ),
        pytest.param(
            {"alpha": Field("a", alias="beta"), "beta": "b"},
            None,
            "MySettingsClass: CLI name 'beta' is claimed by both 'alpha' and 'beta'",
            id="alias_is_other_field_name",
        ),
        pytest.param(
            {"output_dir": "x"},
            {"output-dir": "output-dir"},
            "MySettingsClass: CLI name 'output-dir' is claimed twice by 'output_dir'",
            id="shortcut_is_own_name",
        ),
        pytest.param(
            {"alpha": "a"},
            {"alpha": ["a", "a"]},
            "MySettingsClass: CLI name 'a' is claimed twice by 'alpha'",
            id="repeated_shortcut",
        ),
        pytest.param(
            {"input_dir": "x"},
            {"input_dir": "i"},
            "MySettingsClass: cli_shortcuts target 'input_dir' does not match any CLI argument (did you mean 'input-dir'?)",
            id="target_not_kebab_case",
        ),
        pytest.param(
            {"alpha": "a"},
            {"nope": "n"},
            "MySettingsClass: cli_shortcuts target 'nope' does not match any CLI argument",
            id="unknown_target",
        ),
        pytest.param(
            {"alpha": "a"},
            {"alpha": "-a"},
            "MySettingsClass: invalid shortcut '-a' for 'alpha'; omit the leading dashes",
            id="leading_dash",
        ),
        pytest.param(
            {"alpha": "a"},
            {"alpha": ""},
            "MySettingsClass: invalid shortcut '' for 'alpha'; omit the leading dashes",
            id="empty",
        ),
        pytest.param(
            {"alpha": "a"},
            {"alpha": ["-a", "h"], "nope": "n"},
            "MySettingsClass: invalid shortcut '-a' for 'alpha'; omit the leading dashes\nMySettingsClass: CLI name 'h' for 'alpha' is reserved by argparse\nMySettingsClass: cli_shortcuts target 'nope' does not match any CLI argument",
            id="all_errors_reported",
        ),
    ],
)
def test_check_aliases_fail_invalid_names(
    fields: dict[str, object], cli_shortcuts: dict[str, str | list[str]] | None, expected: str
) -> None:
    """Reserved, colliding, mistargeted or malformed names raise a SettingsError when the class is defined."""
    with pytest.raises(SettingsError) as exc_info:
        define_settings(fields, cli_shortcuts)
    assert str(exc_info.value) == expected
