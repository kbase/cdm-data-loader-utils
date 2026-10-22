"""Tests for XmlToPyArrowSettings field validation and defaults."""

from collections.abc import Callable
from pathlib import Path
from typing import Any, ClassVar
from unittest.mock import patch

import pyarrow as pa
import pytest
import xmlschema
from pydantic import ValidationError
from pydantic_settings import CliApp

import cdm_data_loaders.pipelines.xmltopyarrow.settings as xmltopyarrow_settings_module
from cdm_data_loaders.core.fields import DEFAULTS
from cdm_data_loaders.core.settings import CtsSettings
from cdm_data_loaders.pipelines.xmltopyarrow.settings import (
    ENTRY_TABLE_NAME,
    ENTRY_XML_TAG,
    UNIPROT_PIPELINE_NAME,
    UNIPROT_XSD_FILE,
    UNIREF_PIPELINE_NAME,
    UNIREF_XSD_FILE,
    PresetXmlToPyArrowSettings,
    UniProtXmlToPyArrowSettings,
    UniRefXmlToPyArrowSettings,
    XmlToPyArrowSettings,
)
from cdm_data_loaders.readers.xsd import DEFAULT_MAX_RECURSION_DEPTH, generate_pyarrow_schema, load_schema
from cdm_data_loaders.schemas.uniprot.schema import SCHEMA as UNIPROT_SCHEMA
from cdm_data_loaders.schemas.uniref.schema import SCHEMA as UNIREF_SCHEMA
from tests.conftest import S3
from tests.helpers import settings_to_argv

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
    settings_factory: Callable[..., Any],
    file_glob_arg: str,
) -> None:
    """The file_glob field is exposed as both --file-glob and -g on the command line."""
    base_settings = settings_factory()
    argv = [*settings_to_argv(base_settings, exclude=["file_glob"]), file_glob_arg, "*.xml.gz"]

    settings = CliApp.run(XmlToPyArrowSettings, cli_args=argv)
    assert settings.file_glob == "*.xml.gz"
    assert settings.xml_tag == "book"


def test_uniprot_xmltopyarrow_settings_pass_defaults(
    settings_factory: Callable[..., Any],
) -> None:
    """UniProtXmlToPyArrowSettings defaults dataset_name, table_name, xml_tag, xsd_file and uses precomputed SCHEMA."""
    with patch.object(xmltopyarrow_settings_module, "generate_pyarrow_schema") as mock_gen_schema:
        settings = settings_factory(UniProtXmlToPyArrowSettings)

    mock_gen_schema.assert_not_called()
    assert settings.dataset_name == UNIPROT_PIPELINE_NAME
    assert settings.table_name == "entry"
    assert settings.xml_tag == "entry"
    assert settings.xsd_file == UNIPROT_XSD_FILE
    assert isinstance(settings.xsd_schema, xmlschema.XMLSchema)
    assert settings.pyarrow_schema is UNIPROT_SCHEMA


def test_uniref_xmltopyarrow_settings_pass_defaults(
    settings_factory: Callable[..., Any],
) -> None:
    """UniRefXmlToPyArrowSettings defaults dataset_name, table_name, xml_tag, xsd_file and uses precomputed SCHEMA."""
    with patch.object(xmltopyarrow_settings_module, "generate_pyarrow_schema") as mock_gen_schema:
        settings = settings_factory(UniRefXmlToPyArrowSettings)

    mock_gen_schema.assert_not_called()
    assert settings.dataset_name == UNIREF_PIPELINE_NAME
    assert settings.table_name == "entry"
    assert settings.xml_tag == "entry"
    assert settings.xsd_file == UNIREF_XSD_FILE
    assert isinstance(settings.xsd_schema, xmlschema.XMLSchema)
    assert settings.pyarrow_schema is UNIREF_SCHEMA


@pytest.mark.parametrize(
    ("settings_cls", "expected_dataset_name", "expected_xsd_file", "expected_schema"),
    [
        (UniProtXmlToPyArrowSettings, UNIPROT_PIPELINE_NAME, UNIPROT_XSD_FILE, UNIPROT_SCHEMA),
        (UniRefXmlToPyArrowSettings, UNIREF_PIPELINE_NAME, UNIREF_XSD_FILE, UNIREF_SCHEMA),
    ],
)
def test_subclass_xmltopyarrow_settings_pass_cli_defaults(
    settings_factory: Callable[..., Any],
    *,
    settings_cls: type[XmlToPyArrowSettings],
    expected_dataset_name: str,
    expected_xsd_file: Path,
    expected_schema: pa.Schema,
) -> None:
    """UniProt and UniRef settings subclasses can be parsed via CLI without specifying defaulted arguments."""
    base_settings = settings_factory(CtsSettings)
    settings = CliApp.run(settings_cls, cli_args=settings_to_argv(base_settings))
    assert settings.dataset_name == expected_dataset_name
    assert settings.table_name == "entry"
    assert settings.xml_tag == "entry"
    assert settings.xsd_file == expected_xsd_file
    assert settings.pyarrow_schema is expected_schema


def test_preset_xmltopyarrow_settings_pass_base_defaults(
    settings_factory: Callable[..., Any],
    copy_xsd_fixture: Callable[[str], str],
) -> None:
    """PresetXmlToPyArrowSettings subclasses inherit default table_name and xml_tag as 'entry' and use _PRECOMPUTED_SCHEMA."""

    class CustomPresetSettings(PresetXmlToPyArrowSettings):
        _PRECOMPUTED_SCHEMA: ClassVar[pa.Schema] = UNIREF_SCHEMA

    xsd_path = copy_xsd_fixture(UNIREF_LIKE_XSD)
    settings = settings_factory(
        CustomPresetSettings,
        dataset_name="preset_ds",
        xsd_file=xsd_path,
    )
    assert ENTRY_TABLE_NAME == "entry"
    assert ENTRY_XML_TAG == "entry"
    assert settings.table_name == ENTRY_TABLE_NAME
    assert settings.xml_tag == ENTRY_XML_TAG
    assert settings.pyarrow_schema is UNIREF_SCHEMA
    assert isinstance(settings.xsd_schema, xmlschema.XMLSchema)

    with pytest.raises(ValidationError, match="_PRECOMPUTED_SCHEMA"):
        settings_factory(
            PresetXmlToPyArrowSettings,
            dataset_name="preset_ds",
            xsd_file=xsd_path,
        )
