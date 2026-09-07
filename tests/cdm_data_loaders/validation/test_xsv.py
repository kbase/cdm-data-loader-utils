"""Tests of the xsv-validate adapter."""

import json
import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from cdm_data_loaders.validation.xsv import Status, ValidationResults, validate

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PATCH = "cdm_data_loaders.validation.xsv.subprocess.run"

_CHAR_ERROR_REGEX = r"String should have at (?:least|most) 1 character"


def _data_file(tmp_path: Path, content: str = "col1,col2\nval1,val2\n") -> Path:
    p = tmp_path / "data.csv"
    p.write_text(content)
    return p


def _schema_file(tmp_path: Path, schema: dict | None = None) -> Path:
    schema = schema or {"fields": [{"name": "col1"}, {"name": "col2"}]}
    p = tmp_path / "schema.json"
    p.write_text(json.dumps(schema))
    return p


def _mock_result(returncode: int = 0) -> MagicMock:
    m = MagicMock()
    m.returncode = returncode
    return m


def _get_args(mock_run: MagicMock) -> list[str]:
    return mock_run.call_args.args[0]


# ---------------------------------------------------------------------------
# Input validation — comment_char
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("comment_char", ["", "##", "abc"])
def test_validate_invalid_comment_char(tmp_path: Path, comment_char: str) -> None:
    """validate() raises ValueError when comment_char is not exactly one character."""
    with pytest.raises(ValueError, match=_CHAR_ERROR_REGEX):
        validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path, comment_char=comment_char)


# ---------------------------------------------------------------------------
# Input validation — delimiter
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("delimiter", ["", ",,", "ab"])
def test_validate_invalid_delimiter(tmp_path: Path, delimiter: str) -> None:
    """validate() raises ValueError when delimiter is not exactly one character."""
    with pytest.raises(ValueError, match=_CHAR_ERROR_REGEX):
        validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path, delimiter=delimiter)


# ---------------------------------------------------------------------------
# Input validation — null_strings
# ---------------------------------------------------------------------------


def test_validate_null_strings_contains_empty_string(tmp_path: Path) -> None:
    """validate() raises ValueError when null_strings contains an empty string."""
    with pytest.raises(ValueError, match=_CHAR_ERROR_REGEX):
        validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path, null_strings={""})


def test_validate_null_strings_mixed_valid_and_empty(tmp_path: Path) -> None:
    """validate() raises ValueError even when the empty string is mixed with valid strings."""
    with pytest.raises(ValueError, match=_CHAR_ERROR_REGEX):
        validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path, null_strings={"NA", ""})


# ---------------------------------------------------------------------------
# Input validation — file existence
# ---------------------------------------------------------------------------


def test_validate_missing_file(tmp_path: Path) -> None:
    """validate() raises ValueError when the input file does not exist."""
    with pytest.raises(ValueError, match="Path does not point to a file"):
        validate(tmp_path / "missing.csv", schema=_schema_file(tmp_path), output_path=tmp_path)


# ---------------------------------------------------------------------------
# Status mapping from subprocess return code
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("returncode", "expected_status"),
    [
        (0, Status.VALID),
        (1, Status.INVALID),
        (2, Status.ERROR),
        (127, Status.ERROR),
    ],
)
def test_validate_status_from_returncode(tmp_path: Path, returncode: int, expected_status: Status) -> None:
    """validate() maps subprocess return codes to the correct Status values."""
    with patch(_PATCH, return_value=_mock_result(returncode)):
        result = validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path)
    assert result.status == expected_status


# ---------------------------------------------------------------------------
# Default message and empty results when no summary file is present
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("returncode", "expected_message"),
    [
        (0, "All data are valid"),
        (1, "At least one invalid record"),
        (2, "An error occurred during validation"),
    ],
)
def test_validate_default_message_no_summary_file(tmp_path: Path, returncode: int, expected_message: str) -> None:
    """validate() falls back to the generic status message when no summary file exists."""
    with patch(_PATCH, return_value=_mock_result(returncode)):
        result = validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path)
    assert result.message == expected_message
    assert result.valid_rows is None
    assert result.invalid_rows is None
    assert result.valid_records_file is None
    assert result.invalid_records_file is None
    assert result.errors_file is None


# ---------------------------------------------------------------------------
# Summary file is read and populates results
# ---------------------------------------------------------------------------


def test_validate_full_summary_data(tmp_path: Path) -> None:
    """validate() populates all result fields from a complete summary JSON file."""
    num_valid = 10
    num_invalid = 2
    data_file = _data_file(tmp_path)
    summary = {
        "status_message": "Custom message",
        "valid_rows": num_valid,
        "invalid_rows": num_invalid,
        "valid_records_file": str(tmp_path / "valid.csv"),
        "invalid_records_file": str(tmp_path / "invalid.csv"),
        "errors_file": str(tmp_path / "errors.json"),
    }
    (tmp_path / f"{data_file.name}.summary.json").write_text(json.dumps(summary))

    with patch(_PATCH, return_value=_mock_result(0)):
        result = validate(data_file, schema=_schema_file(tmp_path), output_path=tmp_path, summary=True)

    assert result.message == "Custom message"
    assert result.valid_rows == num_valid
    assert result.invalid_rows == num_invalid
    assert result.valid_records_file == tmp_path / "valid.csv"
    assert result.invalid_records_file == tmp_path / "invalid.csv"
    assert result.errors_file == tmp_path / "errors.json"


@pytest.mark.parametrize(
    ("summary_data", "expected_valid_rows", "expected_invalid_rows"),
    [
        ({"valid_rows": 7}, 7, None),
        ({"invalid_rows": 3}, None, 3),
        ({}, None, None),
    ],
)
def test_validate_partial_summary_data(
    tmp_path: Path,
    summary_data: dict,
    expected_valid_rows: int | None,
    expected_invalid_rows: int | None,
) -> None:
    """validate() handles summary files with missing optional fields."""
    data_file = _data_file(tmp_path)
    (tmp_path / f"{data_file.name}.summary.json").write_text(json.dumps(summary_data))

    with patch(_PATCH, return_value=_mock_result(0)):
        result = validate(data_file, schema=_schema_file(tmp_path), output_path=tmp_path)

    assert result.valid_rows == expected_valid_rows
    assert result.invalid_rows == expected_invalid_rows


# ---------------------------------------------------------------------------
# Return type
# ---------------------------------------------------------------------------


def test_validate_returns_validation_results(tmp_path: Path) -> None:
    """validate() always returns a ValidationResults instance."""
    with patch(_PATCH, return_value=_mock_result()):
        result = validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path)
    assert isinstance(result, ValidationResults)


# ---------------------------------------------------------------------------
# Subprocess args construction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("kwargs", "expected_present", "expected_absent"),
    [
        ({}, [], ["--delimiter", "--missing-header", "--summary-file"]),
        ({"delimiter": ","}, ["--delimiter", ","], []),
        ({"missing_header": True}, ["--missing-header"], []),
        ({"missing_header": False}, [], ["--missing-header"]),
        ({"skip_lines": 5}, ["--skip-lines", "5"], []),
        ({"skip_lines": 0}, ["--skip-lines", "0"], []),
        ({"summary": True}, ["--summary-file"], []),
        ({"summary": False}, [], ["--summary-file"]),
        ({"comment_char": "!"}, ["--comment", "!"], []),
    ],
)
def test_validate_args_flags(
    tmp_path: Path,
    kwargs: dict,
    expected_present: list[str],
    expected_absent: list[str],
) -> None:
    """validate() includes and excludes the correct flags in the subprocess args."""
    with patch(_PATCH, return_value=_mock_result()) as mock_run:
        validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path, **kwargs)
    args = _get_args(mock_run)
    for arg in expected_present:
        assert arg in args
    for arg in expected_absent:
        assert arg not in args


@pytest.mark.parametrize(
    "null_strings",
    [
        {"NA"},
        {"NA", "NULL"},
        {"NA", "NULL", "N/A"},
    ],
)
def test_validate_null_strings_args(tmp_path: Path, null_strings: set[str]) -> None:
    """Each null string produces a --null <value> pair in the subprocess args."""
    with patch(_PATCH, return_value=_mock_result()) as mock_run:
        validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path, null_strings=null_strings)
    args = _get_args(mock_run)
    found = {args[i + 1] for i, a in enumerate(args) if a == "--null"}
    assert found == null_strings


@pytest.mark.parametrize("null_strings", [None, set()])
def test_validate_no_null_args_when_empty(tmp_path: Path, null_strings: set[str] | None) -> None:
    """No --null flags are added when null_strings is None or empty."""
    with patch(_PATCH, return_value=_mock_result()) as mock_run:
        validate(_data_file(tmp_path), schema=_schema_file(tmp_path), output_path=tmp_path, null_strings=null_strings)
    args = _get_args(mock_run)
    assert "--null" not in args


# ---------------------------------------------------------------------------
# Schema as Path vs dict
# ---------------------------------------------------------------------------


def test_validate_schema_as_path_passed_directly(tmp_path: Path) -> None:
    """When schema is a Path, it is passed directly via -s in the subprocess args."""
    schema_file = _schema_file(tmp_path)
    with patch(_PATCH, return_value=_mock_result()) as mock_run:
        validate(_data_file(tmp_path), schema=schema_file, output_path=tmp_path)
    args = _get_args(mock_run)
    assert Path(args[args.index("-s") + 1]) == schema_file


def test_validate_schema_as_dict_written_to_temp_json(tmp_path: Path) -> None:
    """When schema is a dict, it is written to a temporary schema.json and passed via -s."""
    schema_dict = {"fields": [{"name": "col1"}]}
    with patch(_PATCH, return_value=_mock_result()) as mock_run:
        validate(_data_file(tmp_path), schema=schema_dict, output_path=tmp_path)
        args = _get_args(mock_run)
    schema_arg = Path(args[args.index("-s") + 1])
    assert schema_arg.name == "schema.json"


def test_validate_schema_as_dict_content_is_correct(tmp_path: Path) -> None:
    """When schema is a dict, the written temp file contains the correct JSON."""
    schema_dict = {"fields": [{"name": "col1"}]}
    written_content: dict = {}

    def capture_run(args: list[str], **_: object) -> MagicMock:
        schema_path = Path(args[args.index("-s") + 1])
        nonlocal written_content
        written_content = json.loads(schema_path.read_text())
        return _mock_result()

    with patch(_PATCH, side_effect=capture_run):
        validate(_data_file(tmp_path), schema=schema_dict, output_path=tmp_path)

    assert written_content == schema_dict


# =============================================================================
# Integration tests — require xsv-validate.sh on PATH
# =============================================================================


@pytest.fixture(autouse=True)
def _xsv_required_check(request: pytest.FixtureRequest) -> None:
    """Fail for tests that require xsv-validate.sh when it is not installed."""
    if "requires_xsv" not in request.node.keywords:
        return
    if shutil.which("xsv-validate.sh"):
        return
    pytest.fail("xsv-validate.sh not installed; use `not requires_xsv` to exclude", pytrace=False)


# Schema matching the xsv-validator upstream fixture format.
# All columns are strings; age must be digits only; email must be RFC-valid.
_SCHEMA: dict = {
    "$schema": "http://json-schema.org/draft-07/schema#",
    "type": "object",
    "properties": {
        "id": {"type": "string", "minLength": 1},
        "name": {"type": "string", "minLength": 1},
        "email": {"type": "string", "format": "email"},
        "age": {"type": "string", "pattern": "^[0-9]+$"},
    },
    "required": ["id", "name", "email", "age"],
}

_HEADER = "id,name,email,age"


def _write(tmp_path: Path, name: str, content: str) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


# ---------------------------------------------------------------------------
# Status and row counts for different data inputs
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
@pytest.mark.parametrize(
    ("csv_rows", "expected_status", "expected_valid", "expected_invalid"),
    [
        # all rows conform to the schema
        (
            "1,Alice,alice@example.com,30\n2,Bob,bob@example.com,25\n",
            Status.VALID,
            2,
            0,
        ),
        # name is empty (violates minLength:1) on all rows
        (
            "1,,alice@example.com,30\n2,,bob@example.com,25\n",
            Status.INVALID,
            0,
            2,
        ),
        # first row valid, second has bad email and non-numeric age
        (
            "1,Alice,alice@example.com,30\n2,Bob,not-an-email,abc\n",
            Status.INVALID,
            1,
            1,
        ),
        # single valid row
        (
            "1,Alice,alice@example.com,30\n",
            Status.VALID,
            1,
            0,
        ),
        # age is not numeric
        (
            "1,Alice,alice@example.com,thirty\n",
            Status.INVALID,
            0,
            1,
        ),
    ],
    ids=["all_valid", "all_invalid_empty_name", "mixed", "single_valid", "bad_age"],
)
def test_e2e_status_and_counts(
    tmp_path: Path,
    csv_rows: str,
    expected_status: Status,
    expected_valid: int,
    expected_invalid: int,
) -> None:
    """validate() returns the correct status and row counts for various inputs."""
    data_file = _write(tmp_path, "data.csv", f"{_HEADER}\n{csv_rows}")
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, summary=True)

    assert result.status == expected_status
    assert result.valid_rows == expected_valid
    assert result.invalid_rows == expected_invalid


# ---------------------------------------------------------------------------
# Output files exist / don't exist based on status
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
@pytest.mark.parametrize(
    ("csv_rows", "valid_has_rows", "invalid_exists"),
    [
        ("1,Alice,alice@example.com,30\n", True, False),
        ("1,,alice@example.com,30\n", False, True),  # .valid has header only
        ("1,Alice,alice@example.com,30\n1,,x,abc\n", True, True),
    ],
    ids=["all_valid", "all_invalid", "mixed"],
)
def test_e2e_output_files(
    tmp_path: Path,
    csv_rows: str,
    valid_has_rows: bool,
    invalid_exists: bool,
) -> None:
    """validate() creates .valid and .invalid output files as appropriate."""
    data_file = _write(tmp_path, "data.csv", f"{_HEADER}\n{csv_rows}")
    output_dir = tmp_path / "out"

    validate(data_file, schema=_SCHEMA, output_path=output_dir)

    valid_path = output_dir / f"{data_file.name}.valid"
    assert valid_path.exists()
    with valid_path.open("r") as f:
        row_count = sum(1 for _ in f)
    assert (row_count > 1) == valid_has_rows
    assert (output_dir / f"{data_file.name}.invalid").exists() == invalid_exists


# ---------------------------------------------------------------------------
# comment_char — comment lines are stripped before validation
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
@pytest.mark.parametrize(
    ("comment_char", "csv_content"),
    [
        (None, f"# generated by system\n{_HEADER}\n1,Alice,alice@example.com,30\n"),
        ("!", f"! generated by system\n{_HEADER}\n1,Alice,alice@example.com,30\n"),
    ],
    ids=["default_hash", "bang"],
)
def test_e2e_comment_char(tmp_path: Path, comment_char: str | None, csv_content: str) -> None:
    """validate() strips comment lines and validates the remaining data."""
    data_file = _write(tmp_path, "data.csv", csv_content)
    output_dir = tmp_path / "out"

    if comment_char:
        result = validate(data_file, schema=_SCHEMA, output_path=output_dir, comment_char=comment_char)
    else:
        result = validate(data_file, schema=_SCHEMA, output_path=output_dir)
    assert result.status == Status.VALID

    if comment_char:
        result = validate(data_file, schema=_SCHEMA, output_path=output_dir)
    else:
        result = validate(data_file, schema=_SCHEMA, output_path=output_dir, comment_char="?")
    assert result.status == Status.INVALID


# ---------------------------------------------------------------------------
# skip_lines — preamble lines before the header are skipped
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
@pytest.mark.parametrize(
    ("skip_lines", "preamble_lines"),
    [
        (1, "Report generated on 2026-07-09\n"),
        (2, "Report generated on 2026-07-09\nVersion: 1.0\n"),
    ],
    ids=["skip_1", "skip_2"],
)
def test_e2e_skip_lines(tmp_path: Path, skip_lines: int, preamble_lines: str) -> None:
    """validate() skips the specified number of preamble lines before the header."""
    csv_content = f"{preamble_lines}{_HEADER}\n1,Alice,alice@example.com,30\n"
    data_file = _write(tmp_path, "data.csv", csv_content)
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, skip_lines=skip_lines)
    assert result.status == Status.VALID

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, summary=True)
    assert result.status == Status.INVALID


# ---------------------------------------------------------------------------
# delimiter — TSV and explicit comma
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
@pytest.mark.parametrize(
    ("delimiter", "auto_detect", "csv_content"),
    [
        ("\t", True, "id\tname\temail\tage\n1\tAlice\talice@example.com\t30\n"),
        (",", True, f"{_HEADER}\n1,Alice,alice@example.com,30\n"),
        ("&", False, "id&name&email&age\n1&Alice&alice@example.com&30\n"),
    ],
    ids=["tsv", "csv_explicit", "asv_unusual"],
)
def test_e2e_delimiter(tmp_path: Path, delimiter: str, auto_detect: bool, csv_content: str) -> None:
    """validate() correctly parses files with an explicit delimiter."""
    suffix = ".tsv" if delimiter == "\t" else ".csv"
    data_file = _write(tmp_path, f"data{suffix}", csv_content)
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, delimiter=delimiter)
    assert result.status == Status.VALID

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir)
    if auto_detect:
        assert result.status == Status.VALID
    else:
        assert result.status == Status.INVALID

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, delimiter="?")
    assert result.status == Status.INVALID


# ---------------------------------------------------------------------------
# null_strings — custom null values are replaced and trigger schema failures
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
@pytest.mark.parametrize(
    ("null_token", "auto_detect"),
    [
        ("SENTINEL", False),
        ("MISSING", True),
        ("N.D.", False),
    ],
    ids=["sentinel", "missing", "nd"],
)
def test_e2e_custom_null_string_causes_invalid(tmp_path: Path, null_token: str, auto_detect: bool) -> None:
    """A custom null string in a required field makes the row invalid."""
    csv_content = f"{_HEADER}\n1,{null_token},alice@example.com,30\n"
    data_file = _write(tmp_path, "data.csv", csv_content)
    output_dir = tmp_path / "out"

    result = validate(
        data_file,
        schema=_SCHEMA,
        output_path=output_dir,
        null_strings={null_token},
        summary=True,
    )

    assert result.status == Status.INVALID
    assert result.invalid_rows == 1

    result = validate(
        data_file,
        schema=_SCHEMA,
        output_path=output_dir,
        null_strings={"foo", null_token, "bar"},
        summary=True,
    )

    assert result.status == Status.INVALID
    assert result.invalid_rows == 1

    result = validate(
        data_file,
        schema=_SCHEMA,
        output_path=output_dir,
        summary=True,
    )

    if auto_detect:
        assert result.status == Status.INVALID
        assert result.invalid_rows == 1
    else:
        assert result.status == Status.VALID


# ---------------------------------------------------------------------------
# missing_header — header is injected from schema properties
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
def test_e2e_missing_header_valid_data(tmp_path: Path) -> None:
    """validate() injects a header from the schema when missing_header=True."""
    # No header row — just data values in schema property order
    csv_content = "1,Alice,alice@example.com,30\n2,Bob,bob@example.com,25\n"
    data_file = _write(tmp_path, "data.csv", csv_content)
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, missing_header=True)

    assert result.status == Status.VALID
    valid_file = output_dir / f"{data_file.name}.valid"
    assert valid_file.is_file()
    assert valid_file.read_text().splitlines()[0] == _HEADER


@pytest.mark.requires_xsv
def test_e2e_missing_header_invalid_data(tmp_path: Path) -> None:
    """validate() correctly identifies invalid rows when missing_header=True."""
    # age field is not numeric
    csv_content = "1,Alice,alice@example.com,not-a-number\n"
    data_file = _write(tmp_path, "data.csv", csv_content)
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, missing_header=True)

    assert result.status == Status.INVALID
    invalid_file = output_dir / f"{data_file.name}.invalid"
    assert invalid_file.is_file()
    assert invalid_file.read_text().splitlines()[0] == _HEADER


# ---------------------------------------------------------------------------
# schema as dict — same outcome as providing a schema Path
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
def test_e2e_schema_as_dict_produces_same_result(tmp_path: Path) -> None:
    """validate() produces the same status whether schema is a Path or a dict."""
    csv_content = f"{_HEADER}\n1,Alice,alice@example.com,30\n2,Bob,bob@example.com,25\n"
    data_file = _write(tmp_path, "data.csv", csv_content)

    schema_file = tmp_path / "schema.json"
    schema_file.write_text(json.dumps(_SCHEMA))

    result_path = validate(data_file, schema=schema_file, output_path=tmp_path / "out_path")
    result_dict = validate(data_file, schema=_SCHEMA, output_path=tmp_path / "out_dict")

    assert result_path.status == result_dict.status


# ---------------------------------------------------------------------------
# summary file — counts and paths populated correctly
# ---------------------------------------------------------------------------


@pytest.mark.requires_xsv
@pytest.mark.parametrize(
    ("csv_rows", "expected_valid", "expected_invalid"),
    [
        ("1,Alice,alice@example.com,30\n2,Bob,bob@example.com,25\n", 2, 0),
        ("1,,alice@example.com,30\n2,,bob@example.com,25\n", 0, 2),
        ("1,Alice,alice@example.com,30\n2,,not-an-email,abc\n", 1, 1),
    ],
    ids=["all_valid", "all_invalid", "mixed"],
)
def test_e2e_summary_counts(
    tmp_path: Path,
    csv_rows: str,
    expected_valid: int,
    expected_invalid: int,
) -> None:
    """With summary=True, valid_rows and invalid_rows reflect the actual data."""
    data_file = _write(tmp_path, "data.csv", f"{_HEADER}\n{csv_rows}")
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, summary=True)

    assert result.valid_rows == expected_valid
    assert result.invalid_rows == expected_invalid


@pytest.mark.requires_xsv
def test_e2e_summary_file_paths_point_to_existing_files(tmp_path: Path) -> None:
    """With summary=True and mixed data, all path fields in the result point to real files."""
    csv_content = f"{_HEADER}\n1,Alice,alice@example.com,30\n2,,not-an-email,abc\n"
    data_file = _write(tmp_path, "data.csv", csv_content)
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, summary=True)

    assert result.valid_records_file is not None
    assert result.valid_records_file.is_file()
    assert result.invalid_records_file is not None
    assert result.invalid_records_file.is_file()
    assert result.errors_file is not None
    assert result.errors_file.is_file()


@pytest.mark.requires_xsv
def test_e2e_no_summary_file_when_flag_not_set(tmp_path: Path) -> None:
    """With summary=False (default), no summary.json is written and counts are None."""
    data_file = _write(tmp_path, "data.csv", f"{_HEADER}\n1,Alice,alice@example.com,30\n")
    output_dir = tmp_path / "out"

    result = validate(data_file, schema=_SCHEMA, output_path=output_dir, summary=False)

    summary_file = output_dir / f"{data_file.name}.summary.json"
    assert not summary_file.is_file()
    assert result.valid_rows is None
    assert result.invalid_rows is None
