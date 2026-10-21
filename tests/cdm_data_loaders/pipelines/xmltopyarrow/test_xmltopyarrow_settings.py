"""Tests for XmlToPyArrowSettings field validation and defaults."""

from collections.abc import Callable
from pathlib import Path

import pyarrow as pa
import pytest
import xmlschema
from pydantic import ValidationError
from pydantic_settings import CliApp

from cdm_data_loaders.core.fields import DEFAULTS
from cdm_data_loaders.pipelines.xmltopyarrow.settings import XmlToPyArrowSettings
from cdm_data_loaders.readers.xsd import DEFAULT_MAX_RECURSION_DEPTH, generate_pyarrow_schema, load_schema
from tests.conftest import S3

UNIREF_LIKE_XSD: str = "uniref_like.xsd"


def test_xmltopyarrow_settings_pass_defaults(settings_factory: Callable[..., XmlToPyArrowSettings]) -> None:
    """file_glob defaults to '*.xml*'. buffer_size, log_interval, and max_recursion_depth default to expected values."""
    settings = settings_factory()
    assert settings.file_glob == "*.xml*"
    assert settings.buffer_size == DEFAULTS["buffer_size"]
    assert settings.log_interval == DEFAULTS["log_interval"]
    assert settings.max_recursion_depth == DEFAULT_MAX_RECURSION_DEPTH


def test_xmltopyarrow_settings_pass_custom_file_glob(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """file_glob accepts a custom glob pattern."""
    assert settings_factory(file_glob="*.rdf").file_glob == "*.rdf"


def test_xmltopyarrow_settings_pass_custom_xml_tag(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    copy_xsd_fixture: Callable[[str], str],
) -> None:
    """xml_tag accepts any tag name present in the XSD."""
    xsd_file = copy_xsd_fixture(UNIREF_LIKE_XSD)
    settings = settings_factory(xsd_file=xsd_file, xml_tag="entry")
    assert settings.xml_tag == "entry"


def test_xmltopyarrow_settings_fail_missing_xml_tag(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """Omitting xml_tag raises ValidationError."""
    with pytest.raises(ValidationError):
        settings_factory(xml_tag=None)


def test_xmltopyarrow_settings_fail_missing_table_name(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """Omitting table_name raises ValidationError."""
    with pytest.raises(ValidationError):
        settings_factory(table_name=None)


def test_xmltopyarrow_settings_fail_missing_dataset_name(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """Omitting dataset_name raises ValidationError."""
    with pytest.raises(ValidationError):
        settings_factory(dataset_name=None)


def test_xmltopyarrow_settings_fail_missing_xsd_file(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """Omitting xsd_file raises ValidationError."""
    with pytest.raises(ValidationError):
        settings_factory(xsd_file=None)


def test_xmltopyarrow_settings_fail_non_positive_buffer_size(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """buffer_size must be a positive integer."""
    with pytest.raises(ValidationError):
        settings_factory(buffer_size=0)


def test_xmltopyarrow_settings_fail_non_positive_log_interval(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """log_interval must be a positive integer."""
    with pytest.raises(ValidationError):
        settings_factory(log_interval=0)


def test_xmltopyarrow_settings_fail_local_destination_with_s3_output_dir(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """An s3:// output_dir raises ValidationError at construction because output_dir must be a local path."""
    with pytest.raises(ValidationError):
        settings_factory(output_dir="s3://some-bucket/path", use_destination="local_fs")


def test_xmltopyarrow_settings_fail_local_destination_with_s3_flag(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """An s3 destination with a local output_dir is accepted at construction."""
    settings = settings_factory(use_destination=S3)
    assert settings.use_destination == S3


def test_xmltopyarrow_settings_pass_unknown_use_destination_rejected(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """use_destination not present in the dlt config is accepted at construction."""
    settings = settings_factory(use_destination="not_a_destination")
    assert settings.use_destination == "not_a_destination"


def test_xmltopyarrow_settings_pass_schema_properties_populated_when_xsd_file_set(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    copy_xsd_fixture: Callable[[str], str],
) -> None:
    """With xsd_file set, xsd_schema and pyarrow_schema are populated for xml_tag."""
    xsd_file = copy_xsd_fixture(UNIREF_LIKE_XSD)
    settings = settings_factory(xsd_file=xsd_file, xml_tag="entry", table_name="entry")

    assert isinstance(settings.xsd_schema, xmlschema.XMLSchema)
    assert isinstance(settings.pyarrow_schema, pa.Schema)
    expected_schema = generate_pyarrow_schema(
        load_schema(xsd_file),
        "entry",
        max_recursion_depth=DEFAULT_MAX_RECURSION_DEPTH,
    )
    assert settings.pyarrow_schema.equals(expected_schema, check_metadata=True)


def test_xmltopyarrow_settings_pass_custom_max_recursion_depth_controls_pyarrow_schema(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    copy_xsd_fixture: Callable[[str], str],
) -> None:
    """Custom max_recursion_depth is forwarded to generate_pyarrow_schema."""
    xsd_file = copy_xsd_fixture("self_referential.xsd")
    settings_d1 = settings_factory(xsd_file=xsd_file, xml_tag="node", table_name="node", max_recursion_depth=1)
    settings_d2 = settings_factory(xsd_file=xsd_file, xml_tag="node", table_name="node", max_recursion_depth=2)

    assert settings_d1.pyarrow_schema is not None
    assert settings_d2.pyarrow_schema is not None
    assert not settings_d1.pyarrow_schema.equals(settings_d2.pyarrow_schema)


def test_xmltopyarrow_settings_fail_non_positive_max_recursion_depth(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """max_recursion_depth must be a positive integer."""
    with pytest.raises(ValidationError):
        settings_factory(max_recursion_depth=0)


def test_xmltopyarrow_settings_fail_xsd_file_not_found_raises_validation_error(
    settings_factory: Callable[..., XmlToPyArrowSettings],
) -> None:
    """An xsd_file that does not exist raises ValidationError at construction."""
    with pytest.raises(ValidationError, match="Path does not point to a file"):
        settings_factory(xsd_file="does_not_exist.xsd")


def test_xmltopyarrow_settings_fail_xml_tag_not_in_xsd_raises_validation_error(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    copy_xsd_fixture: Callable[[str], str],
) -> None:
    """An xml_tag not present in xsd_file raises ValidationError at construction."""
    xsd_file = copy_xsd_fixture(UNIREF_LIKE_XSD)
    with pytest.raises(ValidationError, match="Could not generate PyArrow schema from XSD"):
        settings_factory(xsd_file=xsd_file, xml_tag="nonExistentTag")


@pytest.mark.parametrize("file_glob_arg", ["--file-glob", "-g"])
def test_xmltopyarrow_settings_pass_xml_tag_shortcut(
    tmp_path: Path,
    library_xsd: Path,
    file_glob_arg: str,
) -> None:
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
        "book",
        "--xml-tag",
        "book",
        "--xsd-file",
        str(library_xsd),
    ]

    settings = CliApp.run(XmlToPyArrowSettings, cli_args=argv)
    assert settings.file_glob == "*.xml.gz"
    assert settings.xml_tag == "book"
