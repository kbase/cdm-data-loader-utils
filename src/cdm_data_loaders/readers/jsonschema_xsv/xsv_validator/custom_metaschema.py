"""Metaschema for JSON Schemas describing XSV file format."""

from typing import Any, Final

import jsonschema
import jsonschema.protocols
from frozendict import frozendict

CUSTOM_META_SCHEMA_URI: Final[str] = "https://json-schema.kbase.us/draft/2020-12-xsv/schema"
X_XSV_CONFIG_KEY: Final[str] = "x-xsv-config"

# schema for the "x-xsv-config" extension block
X_XSV_CONFIG_SCHEMA: frozendict[str, Any] = frozendict(
    {
        "type": "object",
        "properties": {
            "x-has-header": {"type": "boolean", "description": "whether or not the file has a header line"},
            # qsv default: `,`
            "x-delimiter": {
                "type": "string",
                "minLength": 1,
                "maxLength": 1,
                "description": "the delimiter used to demarcate columns in the xsv file",
            },
            # qsv default: none
            "x-comment-char": {
                "type": "string",
                "minLength": 1,
                "maxLength": 1,
                "description": "the comment character",
            },
            # qsv default: `"`
            "x-quote": {"type": "string", "minLength": 1, "maxLength": 1, "description": "quote character"},
            # qsv default: \
            "x-escape": {"type": "string", "minLength": 1, "maxLength": 1, "description": "escape character"},
            "x-quoting-policy": {
                "type": "string",
                "enum": ["all", "necessary", "nonnumeric", "never"],
                "description": "Quoting style to use when writing XSV using qsv. Possible values: all, necessary, nonnumeric and never.\n\nAll: Quotes all fields. \nNecessary: Quotes fields only when necessary - when fields contain a quote, delimiter or record terminator. Quotes are also necessary when writing an empty record (which is indistinguishable from a record with one empty field).\nNonNumeric: Quotes all fields that are non-numeric.\nNever: Never write quotes, even if the resulting XSV file is invalid.",
            },
            "x-null-regex": {
                "type": "string",
                "description": "regular expression capturing the strings that are used to represent NULL values",
            },
            "x-null-cols": {
                "type": "array",
                "items": {"type": "string"},
                "minItems": 1,
                "description": "list of columns in which the null regex is applicable; if not supplied, it is assumed that all columns may contain the NULL value",
            },
        },
        "additionalProperties": False,
        "anyOf": [
            {"required": ["x-has-header"]},
            {"required": ["x-delimiter"]},
            {"required": ["x-comment-char"]},
            {"required": ["x-quote"]},
            {"required": ["x-escape"]},
            {"required": ["x-quoting-policy"]},
            {"required": ["x-null-regex"]},
            {"required": ["x-null-cols"]},
        ],
    }
)

CUSTOM_META_SCHEMA = frozendict(
    {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": CUSTOM_META_SCHEMA_URI,
        "title": "Draft 2020-12 meta-schema extended with x-xsv-config",
        "allOf": [
            # pull in everything the official 2020-12 meta-schema requires/allows
            {"$ref": "https://json-schema.org/draft/2020-12/schema"}
        ],
        "properties": {
            X_XSV_CONFIG_KEY: X_XSV_CONFIG_SCHEMA,
        },
    }
)


# ensure that the metaschema is valid
jsonschema.protocols.Validator.check_schema({**CUSTOM_META_SCHEMA})
