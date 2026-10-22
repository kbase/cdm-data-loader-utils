"""Tests for XmlToDictSettings field validation and defaults."""

from collections.abc import Callable
from pathlib import Path

import pytest
from pydantic import ValidationError
from pydantic_settings import CliApp

from cdm_data_loaders.core.fields import DEFAULTS
from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings
from tests.conftest import S3


def test_xmltodict_settings_pass_defaults(settings_factory: Callable[..., XmlToDictSettings]) -> None:
    """file_glob defaults to '*.xml*'. buffer_size and log_interval default to the common CTS defaults."""
    settings = settings_factory()
    assert settings.file_glob == "*.xml*"
    assert settings.buffer_size == DEFAULTS["buffer_size"]
    assert settings.log_interval == DEFAULTS["log_interval"]


def test_xmltodict_settings_pass_custom_file_glob(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """file_glob accepts a custom glob pattern."""
    assert settings_factory(file_glob="*.rdf").file_glob == "*.rdf"


def test_xmltodict_settings_pass_custom_xml_tag(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """xml_tag accepts any tag name."""
    settings = settings_factory(xml_tag="entry")
    assert settings.xml_tag == "entry"


def test_xmltodict_settings_fail_missing_xml_tag(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """Omitting xml_tag raises ValidationError."""
    with pytest.raises(ValidationError):
        settings_factory(xml_tag=None)


def test_xmltodict_settings_fail_missing_table_name(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """Omitting table_name raises ValidationError."""
    with pytest.raises(ValidationError):
        settings_factory(table_name=None)


def test_xmltodict_settings_fail_missing_dataset_name(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """Omitting dataset_name raises ValidationError."""
    with pytest.raises(ValidationError):
        settings_factory(dataset_name=None)


def test_xmltodict_settings_fail_non_positive_buffer_size(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """buffer_size must be a positive integer."""
    with pytest.raises(ValidationError):
        settings_factory(buffer_size=0)


def test_xmltodict_settings_fail_non_positive_log_interval(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """log_interval must be a positive integer."""
    with pytest.raises(ValidationError):
        settings_factory(log_interval=0)


def test_xmltodict_settings_fail_local_destination_with_s3_output_dir(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """An s3:// output_dir raises ValidationError at construction because output_dir must be a local path."""
    with pytest.raises(ValidationError):
        settings_factory(output_dir="s3://some-bucket/path", use_destination="local_fs")


def test_xmltodict_settings_fail_local_destination_with_s3_flag(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """An s3 destination with a local output_dir is accepted at construction."""
    settings = settings_factory(use_destination=S3)
    assert settings.use_destination == S3


def test_xmltodict_settings_pass_unknown_use_destination_rejected(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """use_destination not present in the dlt config is accepted at construction."""
    settings = settings_factory(use_destination="not_a_destination")
    assert settings.use_destination == "not_a_destination"


@pytest.mark.parametrize("file_glob_arg", ["--file-glob", "-g"])
def test_xmltodict_settings_pass_xml_tag_shortcut(tmp_path: Path, file_glob_arg: str) -> None:
    """The file_glob field is exposed as both --file-glob and -g on the command line."""
    log_config_file = tmp_path / "logging.json"
    log_config_file.write_text('{"version": 1}')
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    argv = [
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--use-destination",
        "local_fs",
        "--log-config-file",
        str(log_config_file),
        file_glob_arg,
        "*.xml.gz",
        "--dataset-name",
        "cli_dataset",
        "--table-name",
        "entry",
        "--xml-tag",
        "entry",
    ]

    settings = CliApp.run(XmlToDictSettings, cli_args=argv)
    assert settings.file_glob == "*.xml.gz"
    assert settings.xml_tag == "entry"
