"""Tests for XmlToDictSettings field validation and defaults."""

import shutil
from collections.abc import Callable
from pathlib import Path

import pytest
from pydantic import ValidationError
from pydantic_settings import CliApp

from cdm_data_loaders.core.fields import DEFAULTS, S3
from cdm_data_loaders.pipelines.core import resolve_cts_settings
from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings

# All schema-driven tests below use tests/data/xsd/uniref_like.xsd, which is already exercised
# (and its list_paths/single_paths independently verified) by test_xsd.py:
# test_find_list_and_single_child_paths_pass_uniref_like_schema.
UNIREF_LIKE_XSD: str = "uniref_like.xsd"


def _copy_xsd_fixture(test_data_dir: Path, input_dir: Path, filename: str) -> str:
    """Copy an XSD fixture from tests/data/xsd into input_dir; return its path."""
    input_dir.mkdir(parents=True, exist_ok=True)
    dest = input_dir / filename
    shutil.copy(test_data_dir / "xsd" / filename, dest)
    return str(dest)


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
    """A local_fs destination with an s3:// output_dir raises ValueError at resolution.

    The check lives in resolve_output_dir and runs when run_cli resolves the
    settings, so constructing them still succeeds.
    """
    settings = settings_factory(output_dir="s3://some-bucket/path", use_destination="local_fs")
    with pytest.raises(ValueError, match="uses protocol 's3', but destination 'local_fs' is configured for 'file'"):
        resolve_cts_settings(settings, {"destination": {"local_fs": {"bucket_url": "/out"}}})


def test_xmltodict_settings_fail_local_destination_with_s3_flag(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """An s3 destination with a local output_dir raises ValueError at resolution."""
    settings = settings_factory(output_dir="/local/out", use_destination=S3)
    with pytest.raises(ValueError, match="uses protocol 'file', but destination 's3' is configured for 's3'"):
        resolve_cts_settings(settings, {"destination": {"s3": {"bucket_url": "s3://bucket/prefix"}}})


def test_xmltodict_settings_pass_unknown_use_destination_rejected(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """use_destination not present in the dlt config raises ValueError at resolution."""
    settings = settings_factory(use_destination="not_a_destination")
    with pytest.raises(ValueError, match="use_destination must be one of"):
        resolve_cts_settings(settings, {"destination": {"local_fs": {"bucket_url": "/out"}}})


def test_xmltodict_settings_pass_xmltodict_args_empty_without_xsd_file(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """Without xsd_file (the default), xmltodict_args contributes no extra xmltodict.parse() kwargs."""
    settings = settings_factory()
    assert settings.xsd_file is None
    assert settings.xmltodict_args == {}


def test_xmltodict_settings_pass_xmltodict_args_has_force_list_callable_when_xsd_file_set(
    tmp_path: Path, test_data_dir: Path, settings_factory: Callable[..., XmlToDictSettings]
) -> None:
    """With xsd_file set, xmltodict_args exposes exactly one force_list callable kwarg."""
    xsd_file = _copy_xsd_fixture(test_data_dir, tmp_path / "input", UNIREF_LIKE_XSD)
    settings = settings_factory(xsd_file=xsd_file)

    assert settings.xmltodict_args.keys() == {"force_list"}
    assert callable(settings.xmltodict_args["force_list"])


def test_xmltodict_settings_pass_xmltodict_args_is_cached(
    tmp_path: Path, test_data_dir: Path, settings_factory: Callable[..., XmlToDictSettings]
) -> None:
    """xmltodict_args is a cached_property: repeated access returns the same dict, not a rebuild."""
    xsd_file = _copy_xsd_fixture(test_data_dir, tmp_path / "input", UNIREF_LIKE_XSD)
    settings = settings_factory(xsd_file=xsd_file)
    assert settings.xmltodict_args is settings.xmltodict_args


def test_xmltodict_settings_fail_xsd_file_not_found_raises_runtime_error(
    settings_factory: Callable[..., XmlToDictSettings],
) -> None:
    """An xsd_file that does not exist raises ValidationError at construction."""
    with pytest.raises(ValidationError, match="Path does not point to a file"):
        settings_factory(xsd_file="does_not_exist.xsd")


@pytest.fixture
def force_list_from_uniref_like_xsd(
    tmp_path: Path, test_data_dir: Path, settings_factory: Callable[..., XmlToDictSettings]
) -> Callable[..., bool]:
    """Build a real XmlToDictSettings from uniref_like.xsd and return its force_list callable."""
    xsd_file = _copy_xsd_fixture(test_data_dir, tmp_path / "input", UNIREF_LIKE_XSD)
    settings = settings_factory(xsd_file=xsd_file)
    return settings.xmltodict_args["force_list"]


@pytest.mark.parametrize(
    ("path", "key", "expected"),
    [
        pytest.param([], "UniRef50", False, id="document-root-is-never-forced-into-a-list"),
        pytest.param([("UniRef50", None)], "entry", True, id="direct-child-of-root-list-path"),
        pytest.param([("UniRef50", None), ("entry", None)], "property", True, id="two-levels-deep-list-path"),
        pytest.param(
            [("UniRef50", None), ("entry", None)],
            "representativeMember",
            False,
            id="two-levels-deep-single-path",
        ),
        pytest.param(
            [("UniRef50", None), ("entry", None), ("representativeMember", None)],
            "dbReference",
            False,
            id="three-levels-deep-single-path-uses-immediate-parent-not-full-ancestor-chain",
        ),
        pytest.param(
            [("UniRef50", None), ("entry", None), ("representativeMember", None), ("dbReference", None)],
            "property",
            True,
            id="four-levels-deep-list-path-uses-immediate-parent-not-full-ancestor-chain",
        ),
        pytest.param(
            [("UniRef50", None)],
            "totallyUnmodelledChild",
            True,
            id="child-not-covered-by-schema-defaults-to-list",
        ),
    ],
)
def test_xmltodict_settings_pass_force_list_classifies_by_immediate_parent(
    force_list_from_uniref_like_xsd: Callable[..., bool],
    path: list[tuple[str, dict[str, str] | None]],
    key: str,
    expected: bool,
) -> None:
    """force_list classifies every (parent, key) pair correctly, using only the immediate parent.

    xmltodict calls force_list(path, key, value) with `path` set to the full ancestor chain (a
    list of (name, attrs) tuples) down to and including the immediate parent -- not just the
    parent's own name -- and with an empty path for the top-level parsed element itself. A
    correct implementation must:

    - Extract the immediate parent's name from the end of that chain, not compare the whole
      chain (which would never match at any nesting depth beyond the root's direct children,
      since find_list_and_single_child_paths stores plain (parent_name, child_name) pairs).
    - Never force the top-level element into a list, since a parsed document only ever has one
      root.
    - Not raise: the ancestor chain is a list, so comparing it directly (unconverted) against a
      set of (parent, key) string pairs raises `TypeError: unhashable type: 'list'`.
    """
    assert force_list_from_uniref_like_xsd(path, key, "some-value") is expected


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
