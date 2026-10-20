"""Schema input parsing with explicit JSON-only and JSON/YAML text policies."""

import json
from pathlib import Path
from typing import Any

import yaml


def _require_mapping(document: object, source: str) -> dict[str, Any]:
    """Reject a parsed document that is not a top-level object."""
    if not isinstance(document, dict):
        msg = f"Parsed schema {source} must be a JSON/YAML object, got {type(document).__name__}"
        raise TypeError(msg)
    return document


def load_schema_text(text: str | bytes, *, allow_yaml: bool = False) -> dict[str, Any]:
    """Parse JSON text, falling back to YAML only when explicitly enabled."""
    try:
        document = json.loads(text)
    except json.JSONDecodeError:
        if not allow_yaml:
            raise
        document = yaml.safe_load(text)
    return _require_mapping(document, "text")


def load_schema_file(path: str | Path) -> dict[str, Any]:
    """Load a JSON or YAML schema from a file, choosing by the `.json` extension.

    :param path: path to the schema file
    :type path: str | Path
    :return: the parsed schema document
    :rtype: dict[str, Any]
    """
    path = Path(path)
    loader = json.loads if path.name.endswith(".json") else yaml.safe_load
    return _require_mapping(loader(path.read_bytes()), f"from {path}")
