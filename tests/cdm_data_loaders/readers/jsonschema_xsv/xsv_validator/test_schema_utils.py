"""Tests for JSON schema-related functions."""

import json
import re
from collections.abc import Callable
from pathlib import Path
from typing import Any, Final

import jsonschema
import jsonschema.exceptions
import jsonschema.validators
import pytest
from pydantic import ValidationError

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.custom_metaschema import (
    CUSTOM_META_SCHEMA,
    CUSTOM_META_SCHEMA_URI,
)
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import (
    ValidatedSchema,
    generate_first_pass_schema,
    generate_header,
    get_schema_parsing_metadata,
    register_xsv_validator,
    validate_jsonschema,
)
from tests.cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.conftest import (
    COLUMNS,
    DELIMITERS,
)

VALID_SCHEMA_URI: Final[str] = "https://json-schema.org/draft/2019-09/schema"
SCHEMA_KEY_VALUE: dict[str, str] = {"$schema": VALID_SCHEMA_URI}

INVALID_TOP_LEVEL_SCHEMA_LIST = [
    pytest.param(["a", "b"], id="schema-is-a-list"),
    pytest.param("just a string", id="schema-is-string"),
    pytest.param(42, id="schema-is-a-number"),
    pytest.param(None, id="schema-is-null"),
    pytest.param(True, id="schema-is-bool"),
]
MISSING_EMPTY_REQUIRED_LIST = [
    pytest.param({}, id="missing-required-key"),
    pytest.param({"required": []}, id="empty-required-list"),
]
INVALID_REQUIRED_LIST = [
    pytest.param({"required": "not-a-list"}, id="required-is-a-string"),
    pytest.param({"required": None}, id="required-is-null"),
    pytest.param({"required": {"a": 1}}, id="required-is-a-dict"),
]


@pytest.fixture
def valid_schema(tmp_path: Path) -> ValidatedSchema:
    """Minimal ValidatedSchema object."""
    return ValidatedSchema(
        jsonschema={"$schema": VALID_SCHEMA_URI, "required": COLUMNS},
        path=tmp_path / "some_file.json",
    )


# def schema_key_value(tmp_path: Path) -> ValidatedSchema:
#     """Schema consisting of just the schema declaration."""
#     return ValidatedSchema(
#         jsonschema=SCHEMA_KEY_VALUE,
#         path=tmp_path / "some_file.json",
#     )


"""validate_jsonschema"""


def test_validate_jsonschema_pass_valid_schema(make_schema_file: Callable[..., Path]) -> None:
    """A well-formed schema, complete with $schema keyword, is returned unchanged."""
    schema = {
        "$schema": VALID_SCHEMA_URI,
        "type": "object",
        "required": ["a"],
        "properties": {"a": {"type": "string"}},
    }
    schema_path = make_schema_file(schema)

    result = validate_jsonschema(schema_path)
    assert isinstance(result, ValidatedSchema)
    assert {**result.jsonschema} == schema


def test_validate_jsonschema_fail_missing_file(tmp_path: Path) -> None:
    """A schema_path that doesn't exist on disk raises FileNotFoundError when read."""
    with pytest.raises(FileNotFoundError, match="No such file or directory"):
        validate_jsonschema(tmp_path / "does-not-exist.json")


def test_validate_jsonschema_fail_invalid_json_content(tmp_path: Path) -> None:
    """A file that isn't valid JSON raises a JSONDecodeError."""
    bad_json_path = tmp_path / "not-json.json"
    bad_json_path.write_text("{this is not: valid json,,,")

    with pytest.raises(json.JSONDecodeError):
        validate_jsonschema(bad_json_path)


@pytest.mark.parametrize("schema", INVALID_TOP_LEVEL_SCHEMA_LIST)
def test_validate_jsonschema_fail_non_dict_schema_or_invalid_schema_raises_error(
    make_schema_file: Callable[..., Path], schema: list | str | int | None
) -> None:
    """A schema whose top-level JSON value isn't a dict raises TypeError."""
    schema_path = make_schema_file(schema)

    with pytest.raises(ValidationError, match="Input should be a valid dictionary"):
        validate_jsonschema(schema_path)


@pytest.mark.parametrize("schema", INVALID_REQUIRED_LIST + MISSING_EMPTY_REQUIRED_LIST)
def test_validate_jsonschema_fail_no_schema_keyword_missing_or_empty_required_raises_value_error(
    make_schema_file: Callable[..., Path], schema: dict[str, Any]
) -> None:
    """A schema with no `$schema` keyword raises a ValueError."""
    schema_path = make_schema_file(schema)

    with pytest.raises(ValueError, match=r"JSON Schema is missing the \$schema keyword"):
        validate_jsonschema(schema_path)


@pytest.mark.parametrize("schema", INVALID_REQUIRED_LIST)
def test_validate_jsonschema_fail_invalid_required_format_error(
    make_schema_file: Callable[..., Path], schema: dict[str, Any]
) -> None:
    """A schema whose required field isn't a list of strings raises SchemaError."""
    schema_path = make_schema_file({**SCHEMA_KEY_VALUE, **schema})

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not of type 'array'"):
        validate_jsonschema(schema_path)


@pytest.mark.parametrize("schema", MISSING_EMPTY_REQUIRED_LIST)
def test_validate_jsonschema_pass_missing_empty_required_list_is_fine_fine_fine(
    make_schema_file: Callable[..., Path], schema: dict[str, Any]
) -> None:
    """A schema whose required field is missing or empty throws a validation error."""
    schema_path = make_schema_file({**SCHEMA_KEY_VALUE, **schema})
    with pytest.raises(ValidationError, match="Could not find any required cols in schema"):
        validate_jsonschema(schema_path)


def test_validate_jsonschema_pass_validator_dgaf_about_unresolvable_schema_uri(
    make_schema_file: Callable[..., Path],
) -> None:
    """A `$schema` value that doesn't correspond to a known JSON Schema draft does not throw an error.

    The jsonschema module uses the most recent draft to validate against if the metaschema is invalid or absent.
    """
    schema = {"$schema": "https://not-a-real-schema.uri", "required": ["a"]}
    schema_path = make_schema_file(schema)
    validated_schema = validate_jsonschema(schema_path)
    assert isinstance(validated_schema, ValidatedSchema)
    assert {**validated_schema.jsonschema} == schema


def test_validate_jsonschema_pass_validator_does_gaf_about_invalid_schema_uri(
    make_schema_file: Callable[..., Path],
) -> None:
    """A `$schema` value that isn't a valid URI throws a SchemaError."""
    schema = {"$schema": "not-a-real-schema-uri", "required": ["a"]}
    schema_path = make_schema_file(schema)
    with pytest.raises(jsonschema.exceptions.SchemaError, match="'not-a-real-schema-uri' is not a 'uri'"):
        validate_jsonschema(schema_path)


def test_validate_jsonschema_fail_invalid_schema_content_raises_schema_error(
    make_schema_file: Callable[[Any, str], Path],
) -> None:
    """A structurally invalid schema fails with a SchemaError."""
    schema = {"$schema": VALID_SCHEMA_URI, "type": "not-a-real-type", "required": ["a"]}
    schema_path = make_schema_file(schema, "schema.json")

    with pytest.raises(
        jsonschema.exceptions.SchemaError, match="'not-a-real-type' is not valid under any of the given schemas"
    ):
        validate_jsonschema(schema_path)


def test_validate_jsonschema_fail_invalid_schema_content_logs_and_raises_schema_error(
    make_schema_file: Callable[..., Path], caplog: pytest.LogCaptureFixture
) -> None:
    """A schema with other invalid fields raises a SchemaError."""
    schema = {"$schema": VALID_SCHEMA_URI, "type": "not-a-real-type", "required": ["a"]}
    schema_path = make_schema_file(schema)

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not valid under any of the given schemas"):
        validate_jsonschema(schema_path)

    assert len(caplog.records) == 1
    assert caplog.records[0].message.startswith("Error validating JSON Schema")


@pytest.mark.parametrize(
    "fixture_name", ["first_pass_schema_path", "post_norm_schema_path", "derived_first_pass_schema_path"]
)
def test_validate_jsonschema_fixtures_are_valid_schemas(fixture_name: str, request: pytest.FixtureRequest) -> None:
    """Ensure that all the JSON Schemas used as fixtures are valid!"""
    schema = request.getfixturevalue(fixture_name)
    validate_jsonschema(schema)


"""get_schema_parsing_metadata"""


def test_get_schema_parsing_metadata_pass_transforms_keys(tmp_path: Path) -> None:
    """Keys prefixed with "x-" are stripped of the prefix and hyphens become underscores."""
    vs = ValidatedSchema(
        jsonschema={
            "$schema": VALID_SCHEMA_URI,
            "required": ["a"],
            "x-xsv-config": {
                "x-has-header": True,
                "x-delimiter": ",",
                "x-comment-char": "#",
                "x-quote": '"',
                "x-escape": "\\",
                "x-quoting-policy": "necessary",
                "x-null-regex": "^(NA|NULL)$",
                "x-null-cols": ["a", "b"],
            },
        },
        path=tmp_path / "schema.json",
    )

    result = get_schema_parsing_metadata(vs)

    assert result == {
        "has_header": True,
        "delimiter": ",",
        "comment_char": "#",
        "quote": '"',
        "escape": "\\",
        "quoting_policy": "necessary",
        "null_regex": "^(NA|NULL)$",
        "null_cols": ["a", "b"],
    }


def test_get_schema_parsing_metadata_fail_no_xsv_config_key(tmp_path: Path) -> None:
    """A schema with no x-xsv-config key raises ValueError."""
    vs = ValidatedSchema(
        jsonschema={"$schema": VALID_SCHEMA_URI, "required": ["a"]},
        path=tmp_path / "schema.json",
    )

    with pytest.raises(ValueError, match="No xsv config information found in schema"):
        get_schema_parsing_metadata(vs)


def test_get_schema_parsing_metadata_fail_empty_xsv_config(tmp_path: Path) -> None:
    """A schema with an empty x-xsv-config object raises ValueError."""
    vs = ValidatedSchema(
        jsonschema={"$schema": VALID_SCHEMA_URI, "required": ["a"], "x-xsv-config": {}},
        path=tmp_path / "schema.json",
    )

    with pytest.raises(ValueError, match="No xsv config information found in schema"):
        get_schema_parsing_metadata(vs)


"""register_xsv_validator"""


def test_register_xsv_validator_pass_registers_validator_for_custom_uri() -> None:
    """After registration, `validator_for` returns a validator class using the custom meta-schema."""
    register_xsv_validator()

    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    assert validator_cls.META_SCHEMA == CUSTOM_META_SCHEMA


def test_register_xsv_validator_pass_idempotent() -> None:
    """Calling register_xsv_validator multiple times does not raise and stays consistent."""
    register_xsv_validator()
    register_xsv_validator()

    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})
    assert validator_cls.META_SCHEMA == CUSTOM_META_SCHEMA


def test_register_xsv_validator_pass_valid_xsv_config_conforms_to_metaschema() -> None:
    """A schema with a well-formed x-xsv-config block validates cleanly against the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {
            "x-has-header": True,
            "x-delimiter": ",",
            "x-quote": '"',
            "x-escape": "\\",
            "x-quoting-policy": "necessary",
        },
    }

    validator_cls.check_schema(schema)  # should not raise


def test_register_xsv_validator_fail_additional_property_not_allowed() -> None:
    """A schema whose x-xsv-config block has an unrecognised key fails the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-not-a-real-option": True},
    }

    with pytest.raises(jsonschema.exceptions.SchemaError, match="Additional properties are not allowed"):
        validator_cls.check_schema(schema)


def test_register_xsv_validator_fail_field_exceeds_max_length() -> None:
    """A schema whose x-delimiter value is longer than one character fails the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-delimiter": ",;"},
    }

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is too long"):
        validator_cls.check_schema(schema)


"""register_xsv_validator - x-null-cols"""


def test_register_xsv_validator_pass_valid_null_cols_conforms_to_metaschema() -> None:
    """A schema with a well-formed x-null-cols array of strings validates against the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": ["col_a", "col_b"]},
    }

    validator_cls.check_schema(schema)  # should not raise


def test_register_xsv_validator_fail_null_cols_wrong_type() -> None:
    """A schema whose x-null-cols value isn't an array fails the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": "col_a"},
    }

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not of type 'array'"):
        validator_cls.check_schema(schema)


def test_register_xsv_validator_fail_null_cols_non_string_items() -> None:
    """A schema whose x-null-cols array contains non-string items fails the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": ["col_a", 42]},
    }

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not of type 'string'"):
        validator_cls.check_schema(schema)


"""validate_jsonschema - x-null-cols"""


def test_validate_jsonschema_pass_valid_null_cols(make_schema_file: Callable[..., Path]) -> None:
    """A schema with a well-formed x-null-cols array of strings validates."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a", "b"],
        "x-xsv-config": {"x-null-cols": ["a", "b"]},
    }
    schema_path = make_schema_file(schema)

    result = validate_jsonschema(schema_path)

    assert isinstance(result, ValidatedSchema)
    assert result.jsonschema == schema


def test_validate_jsonschema_fail_null_cols_wrong_type(make_schema_file: Callable[..., Path]) -> None:
    """A schema whose x-null-cols value isn't an array fails validation."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": "a"},
    }
    schema_path = make_schema_file(schema)

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not of type 'array'"):
        validate_jsonschema(schema_path)


def test_validate_jsonschema_fail_null_cols_non_string_items(make_schema_file: Callable[..., Path]) -> None:
    """A schema whose x-null-cols array contains non-string items fails validation."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": [1, 2, 3]},
    }
    schema_path = make_schema_file(schema)

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not of type 'string'"):
        validate_jsonschema(schema_path)


"""register_xsv_validator - x-null-cols minItems"""


def test_register_xsv_validator_fail_null_cols_empty_array() -> None:
    """A schema whose x-null-cols array is empty fails the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": []},
    }

    with pytest.raises(jsonschema.exceptions.SchemaError, match="should be non-empty"):
        validator_cls.check_schema(schema)


def test_register_xsv_validator_pass_null_cols_single_item() -> None:
    """A schema whose x-null-cols array has exactly one item conforms to the custom meta-schema."""
    register_xsv_validator()
    validator_cls = jsonschema.validators.validator_for({"$schema": CUSTOM_META_SCHEMA_URI})

    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": ["col_a"]},
    }

    validator_cls.check_schema(schema)  # should not raise


"""validate_jsonschema - x-null-cols minItems"""


def test_validate_jsonschema_fail_null_cols_empty_array(make_schema_file: Callable[..., Path]) -> None:
    """A schema whose x-null-cols array is empty fails validation."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": []},
    }
    schema_path = make_schema_file(schema)

    with pytest.raises(jsonschema.exceptions.SchemaError, match="should be non-empty"):
        validate_jsonschema(schema_path)


def test_validate_jsonschema_pass_null_cols_single_item(make_schema_file: Callable[..., Path]) -> None:
    """A schema whose x-null-cols array has exactly one item validates."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-null-cols": ["a"]},
    }
    schema_path = make_schema_file(schema)

    result = validate_jsonschema(schema_path)

    assert isinstance(result, ValidatedSchema)
    assert result.jsonschema == schema


"""validate_jsonschema - custom xsv metaschema"""


def test_validate_jsonschema_pass_conforms_to_custom_xsv_metaschema(
    make_schema_file: Callable[..., Path],
) -> None:
    """A schema declaring the custom xsv $schema URI with a well-formed x-xsv-config block validates."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a", "b"],
        "x-xsv-config": {
            "x-has-header": True,
            "x-delimiter": ",",
            "x-quoting-policy": "necessary",
        },
    }
    schema_path = make_schema_file(schema)

    result = validate_jsonschema(schema_path)

    assert isinstance(result, ValidatedSchema)
    assert result.has_xsv_metaschema
    assert result.has_xsv_parser_config
    assert result.jsonschema == schema


def test_validate_jsonschema_fail_violates_custom_xsv_metaschema_additional_property(
    make_schema_file: Callable[..., Path],
) -> None:
    """A schema whose x-xsv-config block has an unrecognised key fails validation."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-this-is-not-a-real-option": "nope"},
    }
    schema_path = make_schema_file(schema)

    with pytest.raises(jsonschema.exceptions.SchemaError, match="Additional properties are not allowed"):
        validate_jsonschema(schema_path)


def test_validate_jsonschema_fail_violates_custom_xsv_metaschema_bad_type(
    make_schema_file: Callable[..., Path],
) -> None:
    """A schema whose x-has-header value is not a boolean fails validation."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-has-header": "yes"},
    }
    schema_path = make_schema_file(schema)

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not of type 'boolean'"):
        validate_jsonschema(schema_path)


def test_validate_jsonschema_fail_violates_custom_xsv_metaschema_bad_enum_value(
    make_schema_file: Callable[..., Path],
) -> None:
    """A schema whose x-quoting-policy value isn't one of the allowed enum values fails validation."""
    schema = {
        "$schema": CUSTOM_META_SCHEMA_URI,
        "required": ["a"],
        "x-xsv-config": {"x-quoting-policy": "sometimes"},
    }
    schema_path = make_schema_file(schema)

    with pytest.raises(jsonschema.exceptions.SchemaError, match="is not one of"):
        validate_jsonschema(schema_path)


"""generate_header"""


def test_generate_header_pass_uses_defaults(tmp_path: Path, valid_schema: ValidatedSchema) -> None:
    """With no overrides, the header is written to `header.txt` using a tab delimiter."""
    result = generate_header(valid_schema, tmp_path)
    assert result == tmp_path / "header.txt"
    assert result.read_text() == "\t".join(COLUMNS) + "\n"


@pytest.mark.parametrize("delimiter", DELIMITERS)
def test_generate_header_pass_with_various_delimiters(
    tmp_path: Path, delimiter: str, valid_schema: ValidatedSchema
) -> None:
    """Every delimiter supported by SEP_TO_EXT produces a correctly-joined header row."""
    result = generate_header(valid_schema, tmp_path, delimiter=delimiter)
    assert result.read_text() == delimiter.join(COLUMNS) + "\n"


def test_generate_header_pass_with_custom_file_name(tmp_path: Path, valid_schema: ValidatedSchema) -> None:
    """A custom header_file_name is respected in both the returned path and the file written."""
    result = generate_header(valid_schema, tmp_path, header_file_name="custom-header.tsv")
    assert result == tmp_path / "custom-header.tsv"
    assert result.read_text() == "\t".join(COLUMNS) + "\n"


def test_generate_header_pass_overwrites_existing_file(tmp_path: Path, valid_schema: ValidatedSchema) -> None:
    """An existing header file at the target path is overwritten."""
    (tmp_path / "header.txt").write_text("stale content\n")
    result = generate_header(valid_schema, tmp_path)
    assert result.read_text() == "\t".join(COLUMNS) + "\n"


def test_generate_header_fail_missing_target_dir(tmp_path: Path, valid_schema: ValidatedSchema) -> None:
    """A target_dir that doesn't exist fails DirectoryPath validation."""
    missing_dir = tmp_path / "no-such-dir"

    with pytest.raises(ValidationError, match="Path does not point to a directory"):
        generate_header(valid_schema, missing_dir)


def test_generate_header_fail_empty_header_file_name(tmp_path: Path, valid_schema: ValidatedSchema) -> None:
    """An empty header_file_name fails the NonEmptyStr constraint."""
    with pytest.raises(ValidationError, match="String should have at least 1 character"):
        generate_header(valid_schema, tmp_path, header_file_name="")


@pytest.mark.parametrize(
    "delimiter",
    [
        pytest.param("", id="empty-delimiter"),
        pytest.param(",;", id="multi-character-delimiter"),
    ],
)
def test_generate_header_fail_invalid_delimiter_length(
    tmp_path: Path, make_schema_file: Callable[..., Path], delimiter: str
) -> None:
    """A delimiter that isn't exactly one character fails the CharStr constraint."""
    schema_path = make_schema_file({**SCHEMA_KEY_VALUE, "required": COLUMNS})
    err_msg = re.compile("String should have at (mo|lea)st 1 character")
    vs = validate_jsonschema(schema_path)
    with pytest.raises(ValidationError, match=err_msg):
        generate_header(vs, tmp_path, delimiter=delimiter)


"""generate_first_pass_schema"""


def test_generate_first_pass_schema_pass_happy_path(
    post_norm_schema_path: Path, derived_first_pass_schema: dict[str, Any]
) -> None:
    """A first pass schema can be generated from an existing schema file."""
    vs = validate_jsonschema(post_norm_schema_path)
    validated_schema = generate_first_pass_schema(vs)
    assert validated_schema.jsonschema == derived_first_pass_schema


def test_generate_first_pass_schema_pass_bare_minimum(make_schema_file: Callable[[Any, str], Path]) -> None:
    """Test that a schema can be generated from the very barest of bare minimum JSON schema."""
    bare_minimum = {"required": ["this", "that"], **SCHEMA_KEY_VALUE}
    type_dict = {"type": ["string", "null"]}
    schema_path = make_schema_file(bare_minimum, "schema.json")
    vs = validate_jsonschema(schema_path)

    validated_schema = generate_first_pass_schema(vs)
    assert validated_schema.jsonschema == {
        **bare_minimum,
        "type": "object",
        "properties": {"this": type_dict, "that": type_dict},
    }


def test_generate_first_pass_schema_pass_validator_dgaf_about_unresolvable_schema_uri(
    make_schema_file: Callable[[Any, str], Path],
) -> None:
    """An unrecognised `$schema` draft URI is ignored by the validator."""
    schema = {"$schema": "https://not-a-real.schema-uri", "required": ["a"]}
    schema_path = make_schema_file(schema, "schema.json")
    vs = validate_jsonschema(schema_path)

    validated_schema = generate_first_pass_schema(vs)
    assert validated_schema.jsonschema == {
        "$schema": "https://not-a-real.schema-uri",
        "required": ["a"],
        "type": "object",
        "properties": {
            "a": {"type": ["string", "null"]},
        },
    }


def test_generate_first_pass_schema_pass_retains_only_allowed_top_level_keys() -> None:
    """Only $schema, $id, title, and required are copied from the original schema."""
    schema = {
        "$schema": VALID_SCHEMA_URI,
        "$id": "https://example.com/schemas/my-schema.json",
        "title": "My Schema",
        "description": "This should not appear in the first-pass schema",
        "type": "object",
        "additionalProperties": False,
        "required": ["a", "b"],
        "properties": {
            "a": {"type": "string"},
            "b": {"type": "integer"},
        },
    }
    vs = ValidatedSchema(jsonschema=schema)

    validated_schema = generate_first_pass_schema(vs)
    assert validated_schema.jsonschema == {
        "$schema": VALID_SCHEMA_URI,
        "$id": "https://example.com/schemas/my-schema.json",
        "title": "My Schema",
        "type": "object",
        "required": ["a", "b"],
        "properties": {
            "a": {"type": ["string", "null"]},
            "b": {"type": ["string", "null"]},
        },
    }
