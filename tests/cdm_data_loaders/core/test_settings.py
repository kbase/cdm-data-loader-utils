"""Tests for the Settings objects used by DLT pipelines."""

from pathlib import Path
from typing import Any, Self

import dlt
import pytest
from frozendict import frozendict
from pydantic import ValidationError
from pydantic_settings import CliApp, SettingsConfigDict, SettingsError

from cdm_data_loaders.core.fields import (
    DEV_MODE,
    INPUT_DIR,
    OUTPUT_DIR,
    USE_DESTINATION,
    USE_OUTPUT_DIR_FOR_PIPELINE_METADATA,
    VALID_DESTINATIONS,
)
from cdm_data_loaders.core.settings import (
    CLI_SHORTCUTS,
    DEFAULT_SETTINGS_CONFIG_DICT,
    BatchedFileInputSettings,
    CdmDataLoadersBase,
    CtsSettings,
    InputOutputSettings,
    LoggerSettings,
)
from tests.cdm_data_loaders.core.conftest import (
    DEFAULT_BATCH_FILE_SETTINGS_RECONCILED,
    DEFAULT_CTS_SETTINGS_RECONCILED,
    DESTINATION_TO_OUTPUT,
    TEST_BATCH_FILE_SETTINGS,
    TEST_BATCH_FILE_SETTINGS_RECONCILED,
    TEST_CTS_SETTINGS,
    TEST_CTS_SETTINGS_RECONCILED,
    check_settings,
    make_settings,
    make_settings_autofill_config,
)
from tests.helpers import build_cli_arg_specs, make_cli_arg


def _raw_values(settings_cls: type[CtsSettings]) -> frozendict[str, Any] | dict[str, Any]:
    return TEST_BATCH_FILE_SETTINGS if settings_cls is BatchedFileInputSettings else TEST_CTS_SETTINGS


def _reconciled_values(settings_cls: type[CtsSettings]) -> frozendict[str, Any] | dict[str, Any]:
    return (
        TEST_BATCH_FILE_SETTINGS_RECONCILED
        if settings_cls is BatchedFileInputSettings
        else TEST_CTS_SETTINGS_RECONCILED
    )


SETTINGS_CLASSES = [CtsSettings, BatchedFileInputSettings]

INVALID_DESTINATIONS = ["gcs", "filesystem", "", "LocalFs", "S3"]
INVALID_BOOLEAN_VALUES = ["what", "yep", "nope", "2", -1, "", " ", "wtf", None]


def pytest_generate_tests(metafunc: pytest.Metafunc) -> None:
    """Dynamically parametrize CLI-alias tests from the real argparse parser.

    `@pytest.mark.settings_cls.with_args(*classes)` selects which settings classes to
    run against. `@pytest.mark.cli_fields.with_args(*field_names)` optionally restricts
    to a subset of fields (e.g. just the boolean or destination fields); without it,
    every CLI-exposed field on the class is used.
    """
    settings_marker = metafunc.definition.get_closest_marker("settings_cls")
    if settings_marker is None:
        return
    settings_classes = settings_marker.args

    field_marker = metafunc.definition.get_closest_marker("cli_fields")

    wants_option = {"settings_cls", "field_name", "cli_option"}.issubset(metafunc.fixturenames)
    wants_field_only = not wants_option and {"settings_cls", "field_name"}.issubset(metafunc.fixturenames)
    if not (wants_option or wants_field_only):
        return

    argvalues: list[tuple[Any, ...]] = []
    ids: list[str] = []
    for settings_cls in settings_classes:
        specs = build_cli_arg_specs(settings_cls)
        field_names = field_marker.args if field_marker else tuple(specs)
        for field_name in field_names:
            if field_name not in specs:
                continue  # this settings class doesn't define that field
            if wants_option:
                for option in specs[field_name].option_strings:
                    argvalues.append((settings_cls, field_name, option))
                    ids.append(f"{settings_cls.__name__}-{field_name}-{option}")
            else:
                argvalues.append((settings_cls, field_name))
                ids.append(f"{settings_cls.__name__}-{field_name}")

    if wants_option:
        metafunc.parametrize(("settings_cls", "field_name", "cli_option"), argvalues, ids=ids)
    else:
        metafunc.parametrize(("settings_cls", "field_name"), argvalues, ids=ids)


S3 = "is_s3"
OUT = OUTPUT_DIR
RAW = "raw_data_dir"
PIPE = "pipeline_dir"


# manually specify to avoid recapitulating logic
OUTPUT_PATHS: dict[str, dict[str, Any]] = {
    "": {S3: False, OUT: "", RAW: "raw_data", PIPE: ".dlt_conf"},
    "/": {S3: False, OUT: "/", RAW: "/raw_data", PIPE: "/.dlt_conf"},
    # from destination.local_fs
    "/output_dir": {
        S3: False,
        OUT: "/output_dir",
        RAW: "/output_dir/raw_data",
        PIPE: "/output_dir/.dlt_conf",
    },
    "/output/dir": {
        S3: False,
        OUT: "/output/dir",
        RAW: "/output/dir/raw_data",
        PIPE: "/output/dir/.dlt_conf",
    },
    "s3/some/path/": {
        S3: False,
        OUT: "s3/some/path",
        RAW: "s3/some/path/raw_data",
        PIPE: "s3/some/path/.dlt_conf",
    },
    # normalised form of the above
    "s3/some/path": {
        S3: False,
        OUT: "s3/some/path",
        RAW: "s3/some/path/raw_data",
        PIPE: "s3/some/path/.dlt_conf",
    },
    "s3a://bucket/key": {
        S3: True,
        OUT: "s3a://bucket/key",
        RAW: "s3a://bucket/key/raw_data",
        PIPE: None,
    },
    "s3://test/bucket/": {
        S3: True,
        OUT: "s3://test/bucket",
        RAW: "s3://test/bucket/raw_data",
        PIPE: None,
    },
    # normalised from above
    "s3://test/bucket": {
        S3: True,
        OUT: "s3://test/bucket",
        RAW: "s3://test/bucket/raw_data",
        PIPE: None,
    },
    # from destination.s3
    "s3://some/s3/bucket": {
        S3: True,
        OUT: "s3://some/s3/bucket",
        RAW: "s3://some/s3/bucket/raw_data",
        PIPE: None,
    },
}


# a whole load of values that Pydantic will coerce to a boolean
TRUE_FALSE_VALUES = [
    ("0", False),
    ("1", True),
    ("f", False),
    ("false", False),
    ("False", False),
    ("FALSE", False),
    ("n", False),
    ("no", False),
    ("off", False),
    ("on", True),
    ("t", True),
    ("true", True),
    ("True", True),
    ("TRUE", True),
    ("y", True),
    ("yes", True),
    (0, False),
    (1, True),
    (False, False),
    (True, True),
]


# Baseline: no aliases at all
def test_check_aliases_no_additional_fields() -> None:
    """A subclass that adds no fields beyond the base does not raise."""

    class EmptySettings(CdmDataLoadersBase):
        pass

    EmptySettings()


def test_check_aliases_fields_with_no_alias() -> None:
    """Fields with no alias are distinguished by field name only."""

    class PlainSettings(CdmDataLoadersBase):
        field_a: str = "a"
        field_b: str = "b"

        def cli_cmd(self) -> Self:
            return self

    ps = PlainSettings(field_a="whatever", field_b="something")
    assert ps.field_a == "whatever"
    assert ps.field_b == "something"

    ps_from_cmd_line = CliApp.run(PlainSettings, cli_args=["--field-a", "whatever", "--field-b", "something"])
    assert ps == ps_from_cmd_line


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_cli_alias_parses_expected_value(
    settings_cls: type[CtsSettings], field_name: str, cli_option: str
) -> None:
    """Every registered CLI flag for a field parses to the expected value.

    Includes cli_shortcuts, kebab-case versions, validation aliases, etc.
    """
    value = str(_raw_values(settings_cls)[field_name])
    settings = CliApp.run(settings_cls, cli_args=[cli_option, value])
    assert getattr(settings, field_name) == _reconciled_values(settings_cls)[field_name]


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_all_cli_aliases_equivalent(settings_cls: type[CtsSettings], field_name: str) -> None:
    """All CLI flags registered for a single field must produce identical settings objects.

    E.g. `--input-dir foo` and `-i foo` (if `-i` is a configured shortcut) must be
    indistinguishable to the resulting settings instance.
    """
    specs = build_cli_arg_specs(settings_cls)
    option_strings = specs[field_name].option_strings
    if len(option_strings) < 2:
        pytest.skip(f"{settings_cls.__name__}.{field_name} has only one registered CLI flag")

    value = str(_raw_values(settings_cls)[field_name])
    results = [CliApp.run(settings_cls, cli_args=[option, value]) for option in option_strings]
    assert all(r == results[0] for r in results), (
        f"{settings_cls.__name__}.{field_name}: flags {option_strings} produced different results: {results}"
    )


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_configured_cli_shortcuts_are_registered(settings_cls: type[CtsSettings]) -> None:
    """Every shortcut declared in CLI_SHORTCUTS must actually show up on the CLI parser.

    Guards against a `cli_shortcuts` target that doesn't match the (possibly
    kebab-cased) argument name pydantic-settings registers -- which would otherwise
    silently produce a documented shortcut nobody can use.
    """
    specs = build_cli_arg_specs(settings_cls)
    for field_name, shortcuts in CLI_SHORTCUTS.items():
        if field_name not in specs:
            continue
        for shortcut in [shortcuts] if isinstance(shortcuts, str) else shortcuts:
            expected_flag = make_cli_arg(shortcut)
            assert expected_flag in specs[field_name].option_strings, (
                f"shortcut {shortcut!r} for {settings_cls.__name__}.{field_name} was not registered "
                f"(got {specs[field_name].option_strings}); check the cli_shortcuts target matches "
                "the post-kebab-case argument name"
            )


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_env_var_parses_expected_value(
    monkeypatch: pytest.MonkeyPatch, settings_cls: type[CtsSettings], field_name: str
) -> None:
    """Every field can be populated via its env-prefixed environment variable."""
    env_prefix = settings_cls.model_config.get("env_prefix", "")
    monkeypatch.setenv(f"{env_prefix}{field_name}".upper(), str(_raw_values(settings_cls)[field_name]))
    settings = make_settings_autofill_config(settings_cls, {})
    assert getattr(settings, field_name) == _reconciled_values(settings_cls)[field_name]


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_env_var_overrides_dlt_config(
    monkeypatch: pytest.MonkeyPatch, settings_cls: type[CtsSettings], field_name: str
) -> None:
    """Environment variables should take precedence over dlt.config."""
    if field_name != OUTPUT_DIR:
        return

    env_prefix = settings_cls.model_config.get("env_prefix", "")
    override_value = "/env/override/path"
    monkeypatch.setenv(f"{env_prefix}{field_name}".upper(), override_value)

    settings = make_settings_autofill_config(settings_cls, {})
    assert settings.output_dir == override_value


# holy crap, that ain't right
def test_check_aliases_shortcut_matches_field_name_raises() -> None:
    """A shortcut identical to another field's real name is rejected."""
    with pytest.raises(SettingsError, match="' is claimed by both '"):

        class ShortcutMatchesFieldName(CdmDataLoadersBase):
            alpha: str = "a"
            beta: str = "b"

            model_config = SettingsConfigDict(
                **DEFAULT_SETTINGS_CONFIG_DICT,
                cli_shortcuts={"alpha": "beta"},
            )


def test_check_aliases_duplicate_shortcut_raises() -> None:
    """Two fields cannot be assigned the same shortcut."""
    with pytest.raises(SettingsError, match="' is claimed by both '"):

        class DuplicateShortcuts(CdmDataLoadersBase):
            alpha: str = "a"
            beta: str = "b"

            model_config = SettingsConfigDict(
                **DEFAULT_SETTINGS_CONFIG_DICT,
                cli_shortcuts={"alpha": "x", "beta": "x"},
            )


def test_check_aliases_self_referential_shortcut_is_allowed() -> None:
    """A shortcut that maps a field to its own kebab-case name is a harmless no-op."""

    class SelfShortcut(CdmDataLoadersBase):
        output_dir: str = "/tmp"

        model_config = SettingsConfigDict(
            **DEFAULT_SETTINGS_CONFIG_DICT,
            cli_shortcuts={"output_dir": "output-dir"},
        )

    SelfShortcut()


def test_check_aliases_kebab_case_collision_raises() -> None:
    """A field's kebab-cased CLI name cannot collide with another field's shortcut."""
    with pytest.raises(SettingsError, match="'log-config-file' is claimed by both 'log_config_file' and 'other'"):

        class KebabCollision(CdmDataLoadersBase):
            log_config_file: str = "log.json"
            other: str = "x"

            model_config = SettingsConfigDict(
                **DEFAULT_SETTINGS_CONFIG_DICT,
                cli_shortcuts={"other": "log-config-file"},
            )


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_settings_classes_have_no_cli_collisions(settings_cls: type[CtsSettings]) -> None:
    """Building the real parser for each production settings class must not raise an error.

    This is a regression guard against accidentally reusing a shortcut or alias.
    """
    # already ran at class-definition time; call again to be explicit
    settings_cls.check_aliases()
    # raises ArgumentError on any collision
    build_cli_arg_specs(settings_cls)


# CLI App: ignore extra properties
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.usefixtures("patch_dlt_config")
def test_cli_app_run_invalid_params_ignored(settings_cls: type[CtsSettings]) -> None:
    """Test that invalid parameter values are ignored."""
    s = CliApp.run(
        settings_cls,
        cli_args=[
            "--some_random_arg",
            "some value",
            "-q",
            "answer",
        ],
    )
    output = s.model_dump()

    assert "some value" not in output.values()
    assert "answer" not in output.values()


# LoggerSettings


@pytest.mark.parametrize(
    "env_var_name",
    ["CDL_log_config_file", "CDL_LOG_CONFIG_FILE"],
)
def test_logger_settings_alias_accepted_via_environment(
    env_var_name: str, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Both hyphenated and underscored environment variable aliases are resolved to log_config_file."""
    config_file = tmp_path / "log_conf.json"
    config_file.write_text("{}")
    monkeypatch.setenv(env_var_name, str(config_file))
    settings = LoggerSettings()  # pyright: ignore[reportCallIssue]
    assert settings.log_config_file == str(config_file)


# InputOutputSettings
def test_input_output_settings_validate_dir_path() -> None:
    """Test that InputOutputSettings correctly strips trailing slashes."""

    class TestIO(InputOutputSettings):
        log_config_file: str = "log.json"  # pyright: ignore[reportIncompatibleVariableOverride]
        input_dir: str = "/input/"
        output_dir: str = "/output//"

    s = TestIO()
    assert s.input_dir == "/input"
    assert s.output_dir == "/output"
    assert s.log_config_file == "log.json"

    t = InputOutputSettings(input_dir="/input/", output_dir="/output//", log_config_file="log.json")
    for attr in ["input_dir", "output_dir", "log_config_file"]:
        assert getattr(s, attr) == getattr(t, attr)


def test_input_output_settings_preserve_root() -> None:
    """Test that InputOutputSettings preserves the root directory slash."""

    class TestIO(InputOutputSettings):
        log_config_file: str = "log.json"  # pyright: ignore[reportIncompatibleVariableOverride]
        input_dir: str = "/"
        output_dir: str = "/"

    s = TestIO()
    assert s.input_dir == "/"
    assert s.output_dir == "/"
    assert s.log_config_file == "log.json"


# Generic settings tests
@pytest.mark.parametrize(
    ("settings_cls", "args", "expected"),
    [
        # default values
        (CtsSettings, {}, DEFAULT_CTS_SETTINGS_RECONCILED),
        (BatchedFileInputSettings, {}, DEFAULT_BATCH_FILE_SETTINGS_RECONCILED),
        # all args specified
        (CtsSettings, TEST_CTS_SETTINGS, TEST_CTS_SETTINGS_RECONCILED),
        (BatchedFileInputSettings, TEST_BATCH_FILE_SETTINGS, TEST_BATCH_FILE_SETTINGS_RECONCILED),
    ],
)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_all_settings_specified(
    settings_cls: type[CtsSettings], args: dict[str, Any], expected: dict[str, Any]
) -> None:
    """Ensure the CTS settings are set up correctly."""
    s = make_settings_autofill_config(settings_cls, args)
    check_settings(s, expected)


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.usefixtures("patch_dlt_config")
def test_cli_app_run_default_settings(settings_cls: type[CtsSettings]) -> None:
    """Ensure the CTS settings are set up correctly, CLI version."""
    s = CliApp.run(settings_cls)
    expected = (
        DEFAULT_CTS_SETTINGS_RECONCILED if settings_cls == CtsSettings else DEFAULT_BATCH_FILE_SETTINGS_RECONCILED
    )
    check_settings(s, expected)


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_settings_no_dlt_config_error(settings_cls: type[CtsSettings]) -> None:
    """Ensure an error is raised if there is no dlt_config."""
    with pytest.raises(ValidationError, match=r"dlt_config must be defined"):
        make_settings(settings_cls, dlt_config=None)


@pytest.mark.parametrize("invalid_destination_config", [{}, {"destination": {}}])
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_settings_no_destinations_set(
    settings_cls: type[CtsSettings], invalid_destination_config: dict[str, dict[Any, Any]]
) -> None:
    """Ensure that destinations are specified in the dlt config."""
    with pytest.raises(ValueError, match="No valid destinations found in dlt configuration"):
        make_settings(settings_cls, dlt_config=invalid_destination_config)


# same thing but via CliApp.run
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(
    ("dlt_config", "error", "err_msg"),
    [
        (None, ValidationError, "dlt_config must be defined"),
        ({}, ValueError, "No valid destinations found in dlt configuration"),
        ({"destination": {}}, ValueError, "No valid destinations found in dlt configuration"),
    ],
)
@pytest.mark.usefixtures("patch_dlt_config")
def test_cli_app_run_dlt_config_errors(
    settings_cls: type[CtsSettings],
    error: type[Exception],
    err_msg: str,
) -> None:
    """Test all the variants of the Settings fields."""
    with pytest.raises(error, match=err_msg):
        CliApp.run(settings_cls)


# destination tests
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(
    ("use_destination", "output_dir", "should_raise"),
    [
        ("local_fs", "s3://bucket/path", True),
        ("s3", "/local/path", True),
        ("local_fs", "/local/path", False),
        ("s3", "s3://bucket/path", False),
    ],
)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_destination_output_mismatch(
    settings_cls: type[CtsSettings], use_destination: str, output_dir: str, should_raise: bool
) -> None:
    """Mismatch between use_destination and output_dir should raise ValueError."""
    if should_raise:
        with pytest.raises(ValueError, match="Mismatch between output location and use_destination"):
            make_settings_autofill_config(settings_cls, {USE_DESTINATION: use_destination, OUTPUT_DIR: output_dir})
    else:
        s = make_settings_autofill_config(settings_cls, {USE_DESTINATION: use_destination, OUTPUT_DIR: output_dir})
        assert s.use_destination == use_destination
        assert s.output_dir == output_dir


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(USE_DESTINATION, VALID_DESTINATIONS)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_valid_destinations_accepted(use_destination: str, settings_cls: type[CtsSettings]) -> None:
    """Test valid destinations against the settings class."""
    s = make_settings_autofill_config(settings_cls, {USE_DESTINATION: use_destination})
    assert s.use_destination == use_destination


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(USE_DESTINATION, INVALID_DESTINATIONS)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_invalid_destination_raises(use_destination: str, settings_cls: type[CtsSettings]) -> None:
    """Ensure that an unrecognised use_destination raises a ValidationError."""
    with pytest.raises(ValidationError, match=r"use_destination must be one of \['local_fs', 's3'\]"):
        make_settings_autofill_config(settings_cls, {USE_DESTINATION: use_destination})


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_settings_destination_has_no_bucket_url(settings_cls: type[CtsSettings]) -> None:
    """Ensure that destinations have a bucket_url."""
    with pytest.raises(ValueError, match="No bucket_url specified for destination local_fs"):
        make_settings(
            settings_cls,
            dlt_config={"destination": {"local_fs": None}},
            kwargs={USE_DESTINATION: "local_fs"},
        )


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.cli_fields.with_args(DEV_MODE, USE_OUTPUT_DIR_FOR_PIPELINE_METADATA)
@pytest.mark.parametrize(("raw_value", "expected"), TRUE_FALSE_VALUES)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_boolean_cli_variants_accepted(
    settings_cls: type[CtsSettings], field_name: str, cli_option: str, raw_value: Any, expected: bool
) -> None:
    """Booleans can be correctly parsed from CLI args."""
    settings = CliApp.run(settings_cls, cli_args=[cli_option, str(raw_value)])
    assert getattr(settings, field_name) == expected


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.cli_fields.with_args(DEV_MODE, USE_OUTPUT_DIR_FOR_PIPELINE_METADATA)
@pytest.mark.parametrize("bad_value", INVALID_BOOLEAN_VALUES)
@pytest.mark.usefixtures("patch_dlt_config", "field_name")
def test_settings_boolean_cli_variants_rejected(
    settings_cls: type[CtsSettings], cli_option: str, bad_value: Any
) -> None:
    """Invalid booleans are rejected appropriately on the CLI."""
    with pytest.raises(ValidationError, match="Input should be a valid boolean"):
        CliApp.run(settings_cls, cli_args=[cli_option, str(bad_value)])


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.cli_fields.with_args(USE_DESTINATION)
@pytest.mark.parametrize(USE_DESTINATION, VALID_DESTINATIONS)
@pytest.mark.usefixtures("patch_dlt_config", "field_name")
def test_settings_cli_valid_destinations_accepted(
    settings_cls: type[CtsSettings], cli_option: str, use_destination: str
) -> None:
    """Valid destinations are accepted through the CLI."""
    settings = CliApp.run(settings_cls, cli_args=[cli_option, use_destination])
    assert settings.use_destination == use_destination


# boolean fields
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(("input_arg", "value"), TRUE_FALSE_VALUES)
@pytest.mark.parametrize("input_arg_name", [USE_OUTPUT_DIR_FOR_PIPELINE_METADATA, DEV_MODE])
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_boolean_variants_accepted(
    input_arg: str, value: bool, input_arg_name: str, settings_cls: type[CtsSettings]
) -> None:
    """Ensure that each valid boolean value is accepted without error."""
    s = make_settings_autofill_config(settings_cls, {input_arg_name: input_arg})  # type: ignore[reportArgumentType]
    assert getattr(s, input_arg_name) == value


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize("value", INVALID_BOOLEAN_VALUES)
@pytest.mark.parametrize("input_arg_name", [USE_OUTPUT_DIR_FOR_PIPELINE_METADATA, DEV_MODE])
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_invalid_boolean_variants_raises(
    value: bool, input_arg_name: str, settings_cls: type[CtsSettings]
) -> None:
    """Ensure that each invalid boolean value is throws an error."""
    with pytest.raises(ValidationError, match="Input should be a valid boolean"):
        make_settings_autofill_config(settings_cls, {input_arg_name: value})  # type: ignore[reportArgumentType]


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.cli_fields.with_args(USE_DESTINATION)
@pytest.mark.parametrize(USE_DESTINATION, INVALID_DESTINATIONS)
@pytest.mark.usefixtures("patch_dlt_config", "field_name")
def test_settings_cli_invalid_destinations_raises(
    settings_cls: type[CtsSettings], cli_option: str, use_destination: str
) -> None:
    """Invalid destinations are rejected on the CLI."""
    with pytest.raises(ValidationError, match="use_destination must be one of"):
        CliApp.run(settings_cls, cli_args=[cli_option, use_destination])


@pytest.mark.settings_cls.with_args(*SETTINGS_CLASSES)
@pytest.mark.cli_fields.with_args(USE_DESTINATION)
@pytest.mark.usefixtures("field_name")
def test_settings_cli_destination_has_no_bucket_url(
    settings_cls: type[CtsSettings], cli_option: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Incomplete configs are rejected."""
    monkeypatch.setattr(dlt, "config", {"destination": {"local_fs": None}})
    with pytest.raises(ValueError, match="No bucket_url specified for destination local_fs"):
        CliApp.run(settings_cls, cli_args=[cli_option, "local_fs"])


# input and output path coercion
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("/some/path/", "/some/path"),
        ("/some/path//", "/some/path"),
        ("/some/path", "/some/path"),
        ("///", "/"),
        ("/", "/"),
        ("", ""),
    ],
)
@pytest.mark.parametrize("field_name", [INPUT_DIR, OUTPUT_DIR])
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_trailing_slash_stripped(
    settings_cls: type[CtsSettings],
    raw: str,
    expected: str,
    field_name: str,
) -> None:
    """Ensure that validate_dir_path removes trailing slashes but leaves directory slashes intact."""
    s = make_settings_autofill_config(settings_cls, {field_name: raw})
    # output_dir gets filled in with the default if it is falsy
    if field_name == OUTPUT_DIR and raw == "":
        expected = "/output_dir"
    assert getattr(s, field_name) == expected


# values set during reconcile_with_dlt_config
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(USE_DESTINATION, VALID_DESTINATIONS)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_reconcile_with_dlt_config_output_resolved_from_dlt_config_bucket_url(
    settings_cls: type[CtsSettings],
    use_destination: str,
    dlt_config: dict[str, Any],
) -> None:
    """When output_dir is empty, it is populated from dlt config's bucket_url."""
    s = make_settings_autofill_config(settings_cls, {OUTPUT_DIR: "", USE_DESTINATION: use_destination})
    assert s.output_dir == dlt_config[f"destination.{use_destination}.bucket_url"]


# properties derived from self.output_dir: pipeline_dir and raw_data_dir
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(
    OUTPUT_DIR,
    list(OUTPUT_PATHS.keys()),
)
@pytest.mark.parametrize(USE_OUTPUT_DIR_FOR_PIPELINE_METADATA, [True, False])
@pytest.mark.parametrize(USE_DESTINATION, VALID_DESTINATIONS)
@pytest.mark.usefixtures("patch_dlt_config")
def test_settings_generate_pipeline_raw_data_dirs(
    settings_cls: type[CtsSettings],
    output_dir: str,
    use_output_dir_for_pipeline_metadata: bool,
    use_destination: str,
) -> None:
    """Ensure that the correct paths are generated for pipeline and raw data directories.

    Ensure that the destination set in `use_destination` concurs with any output_dir path set.

    Ensure that pipeline directories cannot be set if the output_dir is set to s3.
    """
    make_settings_args = {
        OUTPUT_DIR: output_dir,
        USE_DESTINATION: use_destination,
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: use_output_dir_for_pipeline_metadata,
    }

    expected = {
        **DEFAULT_CTS_SETTINGS_RECONCILED,
        USE_DESTINATION: use_destination,
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: use_output_dir_for_pipeline_metadata,
        OUTPUT_DIR: DESTINATION_TO_OUTPUT[use_destination] if output_dir == "" else OUTPUT_PATHS[output_dir][OUT],
    }
    if settings_cls == BatchedFileInputSettings:
        expected = {**DEFAULT_BATCH_FILE_SETTINGS_RECONCILED, **expected}

    if (OUTPUT_PATHS[expected[OUTPUT_DIR]][S3] and use_destination == "local_fs") or (
        OUTPUT_PATHS[expected[OUTPUT_DIR]][S3] is False and use_destination == "s3"
    ):
        with pytest.raises(ValueError, match="Mismatch between output location and use_destination"):
            make_settings_autofill_config(settings_cls, make_settings_args)
        return

    if use_output_dir_for_pipeline_metadata and OUTPUT_PATHS[expected[OUTPUT_DIR]][S3] is True:
        # can't have pipeline dir on s3
        with pytest.raises(ValueError, match="It is not currently possible to have the pipeline directory on s3"):
            make_settings_autofill_config(settings_cls, make_settings_args)
        return

    s = make_settings_autofill_config(settings_cls, make_settings_args)

    # get the pipeline and raw data dirs from OUTPUT_PATHS
    expected["raw_data_dir"] = OUTPUT_PATHS[expected[OUTPUT_DIR]][RAW]
    # No pipeline_dir if use_output_dir_for_pipeline_metadata is not set
    expected["pipeline_dir"] = (
        OUTPUT_PATHS[expected[OUTPUT_DIR]][PIPE] if use_output_dir_for_pipeline_metadata else None
    )
    check_settings(s, expected)
