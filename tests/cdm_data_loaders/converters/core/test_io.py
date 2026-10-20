"""Tests for schema text and file parsing policies."""

import json
from pathlib import Path

import pytest
import yaml

from cdm_data_loaders.converters.core.io import load_schema_file, load_schema_text


@pytest.mark.parametrize("allow_yaml", [False, True], ids=["json-only", "json-or-yaml"])
@pytest.mark.parametrize("text", ['{"type": "object"}', b'{"type": "object"}'], ids=["string", "bytes"])
def test_load_schema_text_pass_json(text: str | bytes, allow_yaml: bool) -> None:
    """Parse JSON strings and bytes under either text policy."""
    assert load_schema_text(text, allow_yaml=allow_yaml) == {"type": "object"}


@pytest.mark.parametrize("text", ["type: object\n", b"type: object\n"], ids=["string", "bytes"])
def test_load_schema_text_pass_yaml_enabled(text: str | bytes) -> None:
    """Parse YAML only when fallback is enabled."""
    assert load_schema_text(text, allow_yaml=True) == {"type": "object"}


@pytest.mark.parametrize("text", ["type: object\n", "", "{invalid"], ids=["yaml", "empty", "malformed"])
def test_load_schema_text_fail_json_only(text: str) -> None:
    """Propagate JSON errors when YAML fallback is disabled."""
    with pytest.raises(json.JSONDecodeError):
        load_schema_text(text)


@pytest.mark.parametrize(
    ("text", "allow_yaml"),
    [("[1, 2, 3]", False), ("- 1\n- 2\n", True), ("", True), ('"just a string"', False)],
    ids=["json-list", "yaml-list", "empty-yaml-is-none", "json-scalar"],
)
def test_load_schema_text_fail_non_dict_top_level(text: str, allow_yaml: bool) -> None:
    """Reject a parsed document whose top level is not an object."""
    with pytest.raises(TypeError, match="must be a JSON/YAML object"):
        load_schema_text(text, allow_yaml=allow_yaml)


def test_load_schema_text_fail_invalid_yaml() -> None:
    """Propagate YAML syntax errors after JSON parsing fails."""
    with pytest.raises(yaml.YAMLError):
        load_schema_text("type: [", allow_yaml=True)


@pytest.mark.parametrize("as_string", [False, True], ids=["path-object", "path-string"])
@pytest.mark.parametrize(
    ("filename", "text"),
    [
        ("schema.json", '{"type": "object"}'),
        ("schema.yaml", "type: object"),
        ("schema.yml", "type: object"),
        ("schema.txt", "type: object"),
        ("schema.JSON", "type: object"),
        ("schema", "type: object"),
    ],
    ids=["json", "yaml", "yml", "other-extension", "uppercase-extension", "no-extension"],
)
def test_load_schema_file_pass_format_policy(tmp_path: Path, filename: str, text: str, as_string: bool) -> None:
    """Use JSON for lowercase .json files and YAML for all other extensions."""
    path = tmp_path / filename
    path.write_text(text, encoding="utf-8")
    assert load_schema_file(str(path) if as_string else path) == {"type": "object"}


@pytest.mark.parametrize(
    ("filename", "text", "error"),
    [
        ("schema.json", "type: object", json.JSONDecodeError),
        (".json", "type: object", json.JSONDecodeError),
        ("schema.json", "", json.JSONDecodeError),
        ("schema.yaml", "type: [", yaml.YAMLError),
    ],
    ids=["yaml-in-json-file", "json-dotfile", "empty-json", "invalid-yaml"],
)
def test_load_schema_file_fail_invalid_format(
    tmp_path: Path,
    filename: str,
    text: str,
    error: type[Exception],
) -> None:
    """Preserve parser errors instead of changing file formats after failure."""
    path = tmp_path / filename
    path.write_text(text, encoding="utf-8")
    with pytest.raises(error):
        load_schema_file(path)


@pytest.mark.parametrize(
    ("filename", "text"),
    [("schema.json", "[1, 2, 3]"), ("schema.yaml", "- 1\n- 2\n")],
    ids=["json-list", "yaml-list"],
)
def test_load_schema_file_fail_non_dict_top_level(tmp_path: Path, filename: str, text: str) -> None:
    """Reject a file whose parsed top level is not an object."""
    path = tmp_path / filename
    path.write_text(text, encoding="utf-8")
    with pytest.raises(TypeError, match="must be a JSON/YAML object"):
        load_schema_file(path)


def test_load_schema_file_fail_missing_file(tmp_path: Path) -> None:
    """Propagate missing input files."""
    with pytest.raises(FileNotFoundError):
        load_schema_file(tmp_path / "missing.json")
