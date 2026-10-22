"""Tests for the settings classes and helpers in cdm_data_loaders.core.settings."""

from pathlib import Path
from typing import Any

import pytest
from frozendict import frozendict
from pydantic import ValidationError
from pydantic_settings import CliApp, SettingsError

from cdm_data_loaders.core.settings import (
    DEFAULT_SETTINGS_CONFIG_DICT,
    FILE_FORMAT_CLI_SHORTCUTS,
    XML_CLI_SHORTCUTS,
    CdmDataLoadersBase,
    CtsSettings,
    FileFormatIngestSettings,
    LocalInputOutputSettings,
    LoggerSettings,
    XmlIngestSettings,
    default_settings_with_shortcuts,
)
from tests.cdm_data_loaders.conftest import (
    SETTINGS_SOURCES,
    SettingsFactory,
    check_settings,
    cts_expected,
    local_io_expected,
)
from tests.helpers import (
    assert_cli_field_roundtrips,
    assert_no_cli_clashes,
    settings_to_argv,
)


@pytest.mark.parametrize(
    "settings_cls",
    [LoggerSettings, LocalInputOutputSettings, CtsSettings, FileFormatIngestSettings, XmlIngestSettings],
)
def test_settings_pass_no_cli_clashes(settings_cls: type[CdmDataLoadersBase]) -> None:
    """Every core settings class has no CLI argument or alias collisions."""
    assert_no_cli_clashes(settings_cls)


def test_cts_settings_pass_all_params_set(
    test_cts_settings: frozendict[str, Any],
    test_cts_settings_reconciled: frozendict[str, Any],
) -> None:
    """CtsSettings populated with test_cts_settings matches test_cts_settings_reconciled."""
    check_settings(CtsSettings(**test_cts_settings), test_cts_settings_reconciled)


@pytest.mark.parametrize("field_name", list(CtsSettings.model_fields))
def test_cts_settings_pass_cli_field_roundtrips(
    field_name: str,
    test_cts_settings: frozendict[str, Any],
    test_cts_settings_reconciled: frozendict[str, Any],
    test_defaults: dict[str, Any],
) -> None:
    """Every registered CLI option and shortcut for each CtsSettings field roundtrips to the expected value."""
    assert_cli_field_roundtrips(
        CtsSettings,
        field_name,
        test_cts_settings[field_name],
        test_cts_settings_reconciled[field_name],
        test_defaults,
    )


@pytest.mark.parametrize(
    ("kwargs", "extra_expected"),
    [
        pytest.param({}, {}, id="no_args"),
        pytest.param({"cli_prog_name": "load-cts"}, {"cli_prog_name": "load-cts"}, id="prog_name"),
        pytest.param({"cli_prog_name": ""}, {"cli_prog_name": ""}, id="empty_prog_name"),
        pytest.param({"cli_shortcuts": {}}, {"cli_shortcuts": {}}, id="empty_shortcuts"),
        pytest.param(
            {
                "cli_shortcuts": {
                    "input_dir": "i",
                    "dlt_dev_mode": "dev_mode",
                    "use_output_dir_for_pipeline_metadata": ["p", "m"],
                }
            },
            {
                "cli_shortcuts": {
                    "input-dir": "i",
                    "dlt-dev-mode": "dev_mode",
                    "use-output-dir-for-pipeline-metadata": ["p", "m"],
                }
            },
            id="targets_kebab_cased_shortcuts_unchanged",
        ),
        pytest.param(
            {"cli_shortcuts": frozendict({"output_dir": "o"})},
            {"cli_shortcuts": {"output-dir": "o"}},
            id="frozendict_shortcuts",
        ),
        pytest.param(
            {"cli_prog_name": "load-cts", "cli_shortcuts": {"input_dir": "i"}},
            {"cli_prog_name": "load-cts", "cli_shortcuts": {"input-dir": "i"}},
            id="prog_name_and_shortcuts",
        ),
    ],
)
def test_default_settings_with_shortcuts_pass(kwargs: dict[str, Any], extra_expected: dict[str, Any]) -> None:
    """Build a config dict from the defaults plus the given prog name and kebab-cased shortcut targets."""
    assert default_settings_with_shortcuts(**kwargs) == {**DEFAULT_SETTINGS_CONFIG_DICT, **extra_expected}


@pytest.mark.parametrize("source", SETTINGS_SOURCES)
@pytest.mark.parametrize(
    ("values", "expected"),
    [
        pytest.param({}, {"log_config_file": None}, id="default"),
        pytest.param({"log_config_file": "logging.yaml"}, {"log_config_file": "logging.yaml"}, id="set"),
        pytest.param(
            {"log_config_file": "  logging.yaml\t"}, {"log_config_file": "logging.yaml"}, id="whitespace_stripped"
        ),
    ],
)
def test_logger_settings_pass_sources(
    make_settings: SettingsFactory, source: str, values: dict[str, str], expected: dict[str, Any]
) -> None:
    """Build LoggerSettings from each source, with surrounding whitespace stripped."""
    check_settings(make_settings(LoggerSettings, source, values), expected)


@pytest.mark.parametrize("source", SETTINGS_SOURCES)
def test_logger_settings_fail_whitespace_only_log_config_file(make_settings: SettingsFactory, source: str) -> None:
    """Reject a log_config_file that is empty after whitespace stripping."""
    with pytest.raises(ValidationError) as exc_info:
        make_settings(LoggerSettings, source, {"log_config_file": "   "})
    assert {err["loc"][0] for err in exc_info.value.errors()} == {"log_config_file"}


@pytest.mark.parametrize("source", SETTINGS_SOURCES)
def test_local_input_output_settings_pass_sources(
    make_settings: SettingsFactory, source: str, input_dir: Path, output_dir: Path
) -> None:
    """Build LocalInputOutputSettings from each source."""
    settings = make_settings(
        LocalInputOutputSettings, source, {"input_dir": str(input_dir), "output_dir": str(output_dir)}
    )
    check_settings(settings, local_io_expected(input_dir, output_dir))


def test_local_input_output_settings_pass_cli_shortcuts(input_dir: Path, output_dir: Path) -> None:
    """Set input_dir and output_dir with the -i and -o shortcuts."""
    settings = CliApp.run(LocalInputOutputSettings, cli_args=["-i", str(input_dir), "-o", str(output_dir)])
    check_settings(settings, local_io_expected(input_dir, output_dir))


@pytest.mark.parametrize("source", SETTINGS_SOURCES)
@pytest.mark.parametrize(
    ("values", "overrides"),
    [
        pytest.param({"use_destination": "cts"}, {"use_destination": "cts"}, id="defaults"),
        pytest.param(
            {
                "log_config_file": "logging.yaml",
                "buffer_size": "250",
                "dlt_dev_mode": "true",
                "log_interval": "50",
                "use_destination": "s3",
                "use_output_dir_for_pipeline_metadata": "true",
            },
            {
                "log_config_file": "logging.yaml",
                "buffer_size": 250,
                "dlt_dev_mode": True,
                "log_interval": 50,
                "use_destination": "s3",
                "use_output_dir_for_pipeline_metadata": True,
            },
            id="all_fields_set",
        ),
    ],
)
def test_cts_settings_pass_sources(
    make_settings: SettingsFactory,
    source: str,
    values: dict[str, str],
    overrides: dict[str, Any],
    test_defaults: dict[str, Any],
) -> None:
    """Build CtsSettings from each source and check every field, including the computed pipeline_dir."""
    settings = make_settings(
        CtsSettings,
        source,
        {**{k: str(v) for k, v in test_defaults.items()}, **values},
    )
    pipeline_dir = (
        test_defaults["output_dir"] / ".dlt_conf" if overrides.get("use_output_dir_for_pipeline_metadata") else None
    )
    check_settings(settings, cts_expected(**{**test_defaults, **overrides}, pipeline_dir=pipeline_dir))


def test_cts_settings_pass_cli_shortcuts(input_dir: Path, output_dir: Path) -> None:
    """Set fields with the -i, -o, -d and -p shortcuts, as CliApp.run does in production."""
    settings = CliApp.run(
        CtsSettings,
        cli_args=["-i", str(input_dir), "-o", str(output_dir), "-d", "s3", "-p", "true"],
    )
    check_settings(
        settings,
        cts_expected(
            input_dir,
            output_dir,
            use_destination="s3",
            use_output_dir_for_pipeline_metadata=True,
            pipeline_dir=output_dir.resolve() / ".dlt_conf",
        ),
    )


@pytest.mark.parametrize(
    "extra_args",
    [
        pytest.param(["--not-a-setting", "x"], id="unknown_arg"),
        pytest.param(["--buffer_size", "5"], id="snake_case_arg"),
    ],
)
def test_cts_settings_fail_unrecognised_cli_args(test_defaults: dict[str, Any], extra_args: list[str]) -> None:
    """Raise SettingsError, not SystemExit, for CLI args that are not registered."""
    base = CtsSettings(**test_defaults)
    with pytest.raises(SettingsError, match="unrecognized arguments"):
        CliApp.run(CtsSettings, cli_args=[*settings_to_argv(base), *extra_args])


@pytest.mark.parametrize(
    "settings_cls",
    [LocalInputOutputSettings, CtsSettings],
    ids=["local_input_output_settings", "cts_settings"],
)
@pytest.mark.parametrize("set_unprefixed_env", [False, True], ids=["nothing_set", "env_vars_without_prefix"])
def test_settings_fail_missing_dirs(
    monkeypatch: pytest.MonkeyPatch,
    input_dir: Path,
    output_dir: Path,
    settings_cls: type[CdmDataLoadersBase],
    set_unprefixed_env: bool,
) -> None:
    """Fail when input_dir and output_dir are missing; env vars without the CDL_ prefix are ignored."""
    if set_unprefixed_env:
        monkeypatch.setenv("INPUT_DIR", str(input_dir))
        monkeypatch.setenv("OUTPUT_DIR", str(output_dir))
    expected = {
        (("input_dir",), "missing"),
        (("output_dir",), "missing"),
    }
    if settings_cls is CtsSettings:
        expected.add((("use_destination",), "missing"))
    with pytest.raises(ValidationError) as exc_info:
        CliApp.run(settings_cls, cli_args=[])
    assert {(err["loc"], err["type"]) for err in exc_info.value.errors()} == expected


def test_file_format_ingest_settings_pass_cli_shortcuts(test_defaults: dict[str, Any]) -> None:
    """FileFormatIngestSettings supports FILE_FORMAT_CLI_SHORTCUTS (-f) in addition to CtsSettings shortcuts."""
    assert FILE_FORMAT_CLI_SHORTCUTS["loader_file_format"] == "f"
    base = FileFormatIngestSettings(**test_defaults, dataset_name="my_dataset")
    settings = CliApp.run(
        FileFormatIngestSettings,
        cli_args=[
            *settings_to_argv(base, exclude=["loader_file_format"]),
            "-f",
            "jsonl",
        ],
    )
    check_settings(
        settings,
        {
            **cts_expected(**test_defaults, dataset_name="my_dataset"),
            "loader_file_format": "jsonl",
        },
    )


@pytest.mark.parametrize("source", SETTINGS_SOURCES)
def test_xml_ingest_settings_pass_defaults_and_cli_shortcuts(
    make_settings: SettingsFactory, source: str, test_defaults: dict[str, Any]
) -> None:
    """XmlIngestSettings defaults file_glob to '*.xml*' across sources and supports XML_CLI_SHORTCUTS (-f, -g)."""
    assert XML_CLI_SHORTCUTS["loader_file_format"] == "f"
    assert XML_CLI_SHORTCUTS["file_glob"] == "g"
    default_glob_settings = make_settings(
        XmlIngestSettings,
        source,
        {
            **{k: str(v) for k, v in test_defaults.items()},
            "dataset_name": "xml_dataset",
        },
    )
    check_settings(
        default_glob_settings,
        {
            **cts_expected(**test_defaults, dataset_name="xml_dataset"),
            "loader_file_format": "parquet",
            "file_glob": "*.xml*",
        },
    )
    cli_settings = CliApp.run(
        XmlIngestSettings,
        cli_args=[
            *settings_to_argv(default_glob_settings, exclude=["loader_file_format", "file_glob"]),
            "-f",
            "jsonl",
            "-g",
            "*.xml.gz",
        ],
    )
    check_settings(
        cli_settings,
        {
            **cts_expected(**test_defaults, dataset_name="xml_dataset"),
            "loader_file_format": "jsonl",
            "file_glob": "*.xml.gz",
        },
    )
