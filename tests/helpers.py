"""Utilities for running tests."""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml
from pydantic_settings import BaseSettings, CliApp, CliSettingsSource
from pyspark.sql import Row, SparkSession
from pyspark.sql.types import StructType

from cdm_data_loaders.core.settings import CdmDataLoadersBase


def make_cli_arg(arg: str) -> str:
    """Generate the appropriate CLI form for an argument."""
    return f"{'--' if len(arg) > 1 else '-'}{arg.replace('_', '-')}"


@dataclass(frozen=True)
class CliArgSpec:
    """Everything a test needs to drive a single field from the CLI."""

    field_name: str
    option_strings: tuple[str, ...]
    is_flag: bool
    nargs: Any


def build_cli_arg_specs(settings_cls: type[BaseSettings]) -> dict[str, CliArgSpec]:
    """Introspect the argparse parser pydantic-settings builds for `settings_cls`.

    This only constructs the parser; it never parses any arguments, so it's safe to
    call at collection time. Reading directly off the live parser -- rather than from
    a hand-maintained alias table -- means the specs always reflect whatever mix of
    field names, `cli_shortcuts`, `cli_kebab_case`, and validation aliases is currently
    configured on the settings class.
    """
    parser = CliSettingsSource(settings_cls).root_parser
    specs: dict[str, CliArgSpec] = {}
    for action in parser._actions:  # noqa: SLF001 -- argparse has no public introspection API
        if not action.option_strings or action.dest == "help" or action.dest.endswith(":subcommand"):
            continue
        field_info = settings_cls.model_fields.get(action.dest)
        if field_info is not None and field_info.exclude:
            # pydantic-settings always registers a real argparse action per field -- there is no
            # way to omit one outright -- so Field(exclude=True) fields (e.g. CtsSettings.dlt_config)
            # still show up here unless filtered explicitly. Such fields aren't real user-facing CLI
            # options (see DltConfig in cdm_data_loaders.core.fields), so exclude them from the specs.
            continue
        specs[action.dest] = CliArgSpec(
            field_name=action.dest,
            option_strings=tuple(action.option_strings),
            is_flag=action.nargs == 0,
            nargs=action.nargs,
        )
    return specs


def assert_no_cli_clashes(settings_cls: type[CdmDataLoadersBase]) -> None:
    """Fail if any two CLI names on `settings_cls` collide.

    `check_aliases()` already ran at class-definition time and covers field
    names / kebab-case / cli_shortcuts collisions. Building the real argparse
    spec additionally catches collisions arising from validation_alias /
    AliasChoices / AliasPath, which check_aliases doesn't model.
    """
    settings_cls.check_aliases()
    build_cli_arg_specs(settings_cls)  # raises argparse.ArgumentError on any residual clash


def assert_cli_field_roundtrips(
    settings_cls: type[BaseSettings],
    field_name: str,
    raw_value: Any,  # noqa: ANN401
    expected_value: Any,  # noqa: ANN401
    required: dict[str, Any] | None = None,
) -> None:
    """Every registered CLI flag for `field_name` must parse to the same, correct value.

    For models that have required field(s) with no default value, the `required` dictionary should be used
    to populate key-value pairs. The field_name from the model should be used as the key so that it can be
    checked against the field_name under test.
    """
    specs = build_cli_arg_specs(settings_cls)
    option_strings = specs[field_name].option_strings

    required_args = {k: str(v) for k, v in (required or {}).items() if k != field_name}
    required_list = [
        element for item in {make_cli_arg(k): v for k, v in required_args.items()}.items() for element in item
    ]

    results = [
        CliApp.run(model_cls=settings_cls, cli_args=[*required_list, option, str(raw_value)])
        for option in option_strings
    ]

    assert all(getattr(r, field_name) == getattr(results[0], field_name) for r in results), (
        f"{settings_cls.__name__}.{field_name}: flags {option_strings} disagree: {results}"
    )
    assert getattr(results[0], field_name) == expected_value


def create_empty_delta_table(
    spark: SparkSession,
    db: str,
    table: str,
    schema: StructType,
) -> None:
    """Create an empty delta table, initialising the db namespace first."""
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {db}")
    df = spark.createDataFrame([], schema)
    df.write.format("delta").mode("error").saveAsTable(f"{db}.{table}")
    assert df.count() == 0


def assertDataFrameEqual(result_rows: list[Row], expected_rows: list[Row]) -> None:  # noqa: N802
    """Workaround for assertDataFrameEqual from pyspark.testing being broken by pandas 3.0.

    :param result_df: list of dataframe rows, as returned by df.collect()
    :type result_df: list[Row]
    :param expected_df: expected dataframe rows, as returned by df.collect()
    :type expected_df: list[Row]
    """
    results_dict = [r.asDict() for r in result_rows]
    expected_dict = [r.asDict() for r in expected_rows]
    assert len(results_dict) == len(expected_dict)
    for row in results_dict:
        assert row in expected_dict
    for row in expected_dict:
        assert row in results_dict


def write_json_or_yaml(path: Path, data_structure: dict[str, Any]) -> Path:
    """Write ``data_structure`` to ``path`` as JSON or YAML, depending on the file suffix."""
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(data_structure) if path.suffix.lower() == ".json" else yaml.safe_dump(data_structure)
    path.write_text(text, encoding="utf-8")
    return path


def jsonlines_text_to_objects(text: str) -> list[dict[str, Any]]:
    """Parse all non-empty lines of a text string containing JSONlines records.

    :param text: text string with line break-delimited JSONl records
    :type  text: str
    :return: list of JSONL objects
    :rtype : list[dict[str, Any]]
    """
    return [json.loads(line) for line in text.splitlines() if line.strip()]
