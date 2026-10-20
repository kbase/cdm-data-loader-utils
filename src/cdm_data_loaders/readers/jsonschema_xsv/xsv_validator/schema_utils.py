"""JSON Schema-related functions."""

import json
from logging import Logger, getLogger
from pathlib import Path
from typing import Annotated, Any

import jsonschema
import jsonschema.exceptions
import jsonschema.validators
from jsonschema import Draft202012Validator, validators
from pydantic import (
    BaseModel,
    DirectoryPath,
    StringConstraints,
    computed_field,
    field_validator,
    validate_call,
)
from referencing import Registry, Resource

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.custom_metaschema import (
    CUSTOM_META_SCHEMA,
    CUSTOM_META_SCHEMA_URI,
    X_XSV_CONFIG_KEY,
)

NonEmptyStr = Annotated[str, StringConstraints(min_length=1)]
CharStr = Annotated[str, StringConstraints(min_length=1, max_length=1)]

logger: Logger = getLogger(__name__)


class ValidatedSchema(BaseModel):
    """A JSON Schema to be used for parsing and validating XSV data."""

    jsonschema: dict[str, Any]

    @computed_field
    @property
    def has_xsv_parser_config(self) -> bool:
        """Whether or not this schema has custom fields for parsing XSV files."""
        return len(self.jsonschema.get(X_XSV_CONFIG_KEY, {}).keys()) > 0

    @computed_field
    @property
    def has_xsv_metaschema(self) -> bool:
        """Whether or not this schema has the custom metaschema URI indicating it holds XSV metadata."""
        return self.jsonschema.get("$schema") == CUSTOM_META_SCHEMA_URI

    @computed_field
    @property
    def required_cols(self) -> list[str]:
        """Retrieve the list of headers from the top-level `required` field."""
        return self.jsonschema.get("required", [])

    @field_validator("jsonschema", mode="after")
    @classmethod
    def validate_jsonschema(cls, schema: dict[str, Any]) -> dict[str, Any]:
        """Validate the schema in the jsonschema field.

        :param schema: a JSON schema
        :type schema: dict[str, Any]
        :raises ValueError: if there is no `$schema` field to identify the metaschema used
        :raises ValueError: if there are no `required` fields to indicate the CSV columns and order
        :return: validated schema
        :rtype: dict[str, Any]
        """
        register_xsv_validator()
        # not required by the spec, but we have higher standards
        if "$schema" not in schema:
            err_msg = "JSON Schema is missing the $schema keyword"
            raise ValueError(err_msg)

        # retrieve the appropriate validator for the schema and ensure it is valid
        # if the $schema value is invalid, jsonschema will use the most recent draft by default and emit a warning
        validator = jsonschema.validators.validator_for(schema)
        try:
            validator.check_schema(schema)
        except (jsonschema.exceptions.SchemaError, jsonschema.exceptions.ValidationError):
            logger.exception("Error validating JSON Schema")
            raise

        required_cols = schema.get("required")
        if not isinstance(required_cols, list) or not required_cols:
            err_msg = "Could not find any required cols in schema"
            raise ValueError(err_msg)
        return schema


def register_xsv_validator() -> None:
    """Register the custom validator for validating JSON Schema with embedded XSV parsing metadata."""
    # register the new meta-schema and build a validator class for it
    Registry().with_resource(
        CUSTOM_META_SCHEMA_URI,
        Resource.from_contents(CUSTOM_META_SCHEMA),
    )

    # extend the Draft202012Validator metaschema to include the extra XSV parsing info
    # the metaschema only affects validation of schemas, not of data conforming to the schema
    xsv_validator = validators.extend(Draft202012Validator)
    xsv_validator.META_SCHEMA = CUSTOM_META_SCHEMA

    # Register it so `jsonschema.validate(...)`/`validator_for(...)` automatically
    # pick xsv_validator whenever a schema declares this "$schema" URI.
    validators.validates(CUSTOM_META_SCHEMA_URI)(xsv_validator)


def validate_jsonschema(schema_path: Path) -> ValidatedSchema:
    """Ensure that a given data structure is a valid JSON schema.

    :param schema: JSON schema, loaded as a python data structure
    :type schema: dict[str, Any]
    :raises jsonschema.exceptions.SchemaError: if the schema is invalid
    :return: validated JSON Schema
    :rtype: dict[str, Any]
    """
    schema = json.loads(schema_path.read_bytes())
    register_xsv_validator()
    return ValidatedSchema(jsonschema=schema)


@validate_call
def generate_first_pass_schema(validated_schema: ValidatedSchema) -> ValidatedSchema:
    """Given a full schema file, generate a schema to perform a loose first-pass validation with.

    The first pass schema is used for verifying that the top-level columns are correct; no further checks
    are performed.

    :param schema: parsed, validated JSON schema
    :type schema: ValidatedSchema
    :return: the validated first pass schema
    :rtype: ValidatedSchema
    """
    first_pass_schema = {
        k: v for k, v in validated_schema.jsonschema.items() if k in ["$schema", "$id", "title", "required"]
    }

    first_pass_schema["type"] = "object"
    first_pass_schema["properties"] = {req: {"type": ["string", "null"]} for req in validated_schema.required_cols}

    # retrieve the appropriate validator for the schema and ensure it is valid
    validator = jsonschema.validators.validator_for(first_pass_schema)
    validator.check_schema(first_pass_schema)

    return ValidatedSchema(jsonschema=first_pass_schema)


@validate_call
def generate_header(
    validated_schema: ValidatedSchema,
    target_dir: DirectoryPath,
    header_file_name: NonEmptyStr = "header.txt",
    delimiter: CharStr = "\t",
) -> Path:
    r"""Generate a header file for an xSV file from a JSON Schema.

    If no file name is supplied, the header file is saved as `header.txt`.

    The top level required properties are assumed to be the columns of the xSV file.

    :param schema: parsed, validated JSON schema
    :type schema: ValidatedSchema
    :param target_dir: directory in which to save the file
    :type target_dir: DirectoryPath
    :param delimiter: delimiter to use for the headers, defaults to "\t"
    :type delimiter: CharStr, optional
    :param header_file_name: name for the header file, defaults to "header.txt"
    :type header_file_name: str, optional
    :raises TypeError: if the schema is not in the correct format (dictionary)
    :raises ValueError: if the schema file has no `required` cols
    :return: path to the newly-created header.txt file
    :rtype: Path
    """
    header_file_path = target_dir / header_file_name
    header_file_path.write_text(delimiter.join(validated_schema.required_cols) + "\n")
    return header_file_path


def get_schema_parsing_metadata(validated_schema: ValidatedSchema) -> dict[str, Any]:
    """Read in the schema parser config information from the supplied schema.

    :param schema: parsed, validated JSON schema
    :type schema: ValidatedSchema
    :return: parsed xsv schema metadata
    :rtype: dict[str, Any]
    """
    if not validated_schema.has_xsv_parser_config:
        err_msg = "No xsv config information found in schema"
        raise ValueError(err_msg)

    return {
        k.replace("x-", "").replace("-", "_"): v
        for k, v in validated_schema.jsonschema.get(X_XSV_CONFIG_KEY, {}).items()
    }
