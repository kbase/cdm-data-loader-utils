"""Reader for XSV data already cleaned and validated by `xsv_validator.qsv`.

Reads rows from a file produced by `xsv_validator.qsv.clean_validate_file`, casting each field
to the python type declared for it in the schema's `properties` block. Also resolves the
qsv cleaning/validation parameters (delimiter, quote, escape, comment character, null handling)
from a schema's `x-xsv-config` extension block.
"""

import csv
from collections.abc import Callable, Generator
from pathlib import Path
from typing import Any, Final

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import (
    ValidatedSchema,
    get_schema_parsing_metadata,
)

JSON_SCHEMA_TYPE_CASTERS: Final[dict[str, Callable[[str], Any]]] = {
    "integer": int,
    "number": float,
    "boolean": lambda value: value.strip().lower() in {"true", "1"},
}

# maps x-xsv-config metadata keys (after get_schema_parsing_metadata strips the "x-" prefix) to
# the CleanerValidatorArgs kwarg they configure. x-quoting-policy is omitted: it configures XSV
# output, not the cleaning/validation pass driven by this config.
PARSING_CONFIG_KEYS: Final[dict[str, str]] = {
    "delimiter": "delimiter",
    "comment_char": "comment_char",
    "quote": "quote",
    "escape": "escape",
    "null_regex": "null_regex",
    "null_cols": "null_regex_cols",
}


def resolve_xsv_parsing_config(validated_schema: ValidatedSchema) -> dict[str, Any]:
    """Derive `CleanerValidatorArgs` kwargs from a schema's `x-xsv-config` extension block.

    Schemas without an `x-xsv-config` block resolve to an empty dict, letting
    `CleanerValidatorArgs` fall back to its own defaults (tab-delimited, `#`-commented, header
    present).

    :param validated_schema: parsed, validated JSON Schema for the XSV data
    :type validated_schema: ValidatedSchema
    :return: kwargs accepted by `CleanerValidatorArgs`
    :rtype: dict[str, Any]
    """
    if not validated_schema.has_xsv_parser_config:
        return {}

    metadata = get_schema_parsing_metadata(validated_schema)
    config = {arg_name: metadata[key] for key, arg_name in PARSING_CONFIG_KEYS.items() if key in metadata}
    if "has_header" in metadata:
        config["missing_header"] = not metadata["has_header"]
    return config


def cast_xsv_value(value: str, type_spec: str | list[str] | None) -> Any:  # noqa: ANN401
    """Cast one raw string field value to the python type declared for it in a JSON Schema.

    An empty string is treated as null whenever "null" is one of the declared types, or no type
    is declared at all. Types without a caster (`string`, or unrecognised) pass through unchanged.

    :param value: raw string field value from a validated XSV row
    :type value: str
    :param type_spec: the property's declared JSON Schema type(s)
    :type type_spec: str | list[str] | None
    :return: the cast value
    :rtype: Any
    """
    types = [type_spec] if isinstance(type_spec, str) else list(type_spec or [])
    if value == "" and ("null" in types or not types):
        return None
    for candidate in types:
        caster = JSON_SCHEMA_TYPE_CASTERS.get(candidate)
        if caster is not None:
            return caster(value)
    return value


def cast_xsv_row(row: dict[str, str], properties: dict[str, Any]) -> dict[str, Any]:
    """Cast every field in a raw XSV row to the type declared for it in `properties`.

    :param row: one row of raw string values, keyed by column name
    :type row: dict[str, str]
    :param properties: the schema's `properties` block, keyed by column name
    :type properties: dict[str, Any]
    :return: the row with each value cast to its declared type
    :rtype: dict[str, Any]
    """
    return {column: cast_xsv_value(value, properties.get(column, {}).get("type")) for column, value in row.items()}


def read_validated_rows(
    xsv_path: Path,
    properties: dict[str, Any],
    delimiter: str = "\t",
    quote: str | None = None,
    escape: str | None = None,
) -> Generator[dict[str, Any], Any, Any]:
    """Read a cleaned, validated XSV file, casting each row to its declared schema types.

    The file is expected to carry a header row matching the schema's column names, as produced
    by `xsv_validator.qsv.clean_validate_file`.

    :param xsv_path: path to the cleaned, validated XSV file
    :type xsv_path: Path
    :param properties: the schema's `properties` block, keyed by column name
    :type properties: dict[str, Any]
    :param delimiter: delimiter used in the file, defaults to tab
    :type delimiter: str, optional
    :param quote: quote character used in the file; defaults to `"` when unset
    :type quote: str | None, optional
    :param escape: escape character used in the file; when unset, quotes are escaped by doubling
    :type escape: str | None, optional
    :yield: rows cast to their declared schema types
    :rtype: Generator[dict[str, Any], Any, Any]
    """
    with xsv_path.open(newline="") as fh:
        reader = csv.DictReader(
            fh,
            delimiter=delimiter,
            quotechar=quote or '"',
            escapechar=escape,
            doublequote=escape is None,
        )
        for row in reader:
            yield cast_xsv_row(row, properties)
