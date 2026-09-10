"""Tests for the xsv_validator.qsv module."""

import csv
import re
import shutil
import subprocess
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Final, NamedTuple
from unittest.mock import MagicMock

import pytest
from pydantic import ValidationError

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator import qsv
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.helpers import (
    CLEANED_SUFFIX,
    HEADER_SUFFIX,
    NORM_SUFFIX,
    SEP_TO_EXT,
    VALID_SUFFIX,
    CleanerValidatorArgs,
    generate_qsv_validate_file_names,
    non_header_lines_present,
)
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.qsv import (
    clean_validate_file,
    qsv_check,
    run_qsv_input,
    run_qsv_null_replacement,
    run_qsv_validate,
)
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import ValidatedSchema
from tests.cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.conftest import (
    COLUMNS,
    DELIMITERS,
    FAKE_QSV_CMD,
    VALID_ROWS,
    WriteFile,
    _write_cleaned_input_file,
    _write_header_input_file,
    _write_validator_input_file,
    build_xsv_content,
    interleave_comments,
    parse_xsv,
    snapshot_dir,
)

# One ragged row (missing the trailing `string` field) among otherwise-valid rows: recoverable.
PARTIAL_RAGGED_ROWS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "3.14", "true", "key:value1"],
    ["3", "2023-02-20", "1.11", "false"],
    ["4", "2023-03-25", "2.22", "true", "key:value3"],
]

# Every row here is ragged (too few or too many fields): unrecoverable.
ALL_RAGGED_ROWS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "3.14", "true"],
    ["3", "2023-02-20", "1.11", "false", "key:value2", "extra"],
]


ROWS_INVALID_NUMBER: Final[list[list[str]]] = [
    ["10", "2023-01-15", "3.14", "true", "key:value1"],
    ["-1", "2022-01-01", "", "false", "a:b"],
]
ROWS_INVALID_DATE: Final[list[list[str]]] = [
    ["2", "not-a-date", "3.14", "true", "key:value1"],
    ["4", "2501-51-61", "1", "false", "nug:get"],
]
ROWS_INVALID_STRING: Final[list[list[str]]] = [
    ["2", "2023-01-15", "3.14", "true", "badstring"],
    ["4", "2020-05-15", "", "false", "whatever"],
]
ROWS_PADDED_WHITESPACE: Final[list[list[str]]] = [
    [" 2 ", " 2023-01-15 ", "3.14", "true", " key:value1 "],
    [" 2", " 2023-01-15    ", "   ", "   true", " key:value2                            \r"],
]

INVALID_UTF8_WITH_NULLS: Final[bytes] = b"col1\tcol2\tNA\tNA\t\xff\n\xff\xfe\tcol2\t\xff\xfe\tNA\t\xee\n"
INVALID_UTF8_NO_NULLS: Final[bytes] = b"col1\tcol2\t\xfe\t \t\xff\n\xff\xfe\tcol2\t\xff\xfe\t1234\t\xee\n"
INVALID_UTF8_MIXED: Final[bytes] = INVALID_UTF8_NO_NULLS + INVALID_UTF8_WITH_NULLS


class NullValueCase(NamedTuple):
    """A null placeholder value paired with the regex that should match (and only match) it.

    :param placeholder: literal text appearing in the source data that represents a null value
    :param regex: Rust-compatible regex passed as `null_regex` to match `placeholder` exactly
    """

    placeholder: str
    regex: str


# Generic null value
NULL_PLACEHOLDER: Final[str] = "NA"
NULL_REGEX: Final[str] = r"^NA$"


# Anchored so that each regex matches only a field consisting entirely of the placeholder text,
# not fields that merely contain it as a substring (e.g. "national" should not match "na").
NULL_REGEX_CASES: Final[list[NullValueCase]] = [
    NullValueCase(placeholder="null", regex=r"^null$"),
    NullValueCase(placeholder="NULL", regex=r"^NULL$"),
    NullValueCase(placeholder="None", regex=r"^None$"),
    NullValueCase(placeholder="N/A", regex=r"^N/A$"),
    NullValueCase(placeholder="na", regex=r"^na$"),
    NullValueCase(placeholder="NA", regex=r"^NA$"),
    # Contains a regex metacharacter ('.') that must be escaped for the match to be exact.
    NullValueCase(placeholder="n.a.", regex=r"^n\.a\.$"),
]

# A single regex, case-insensitively matching several differently-cased spellings of "na" that
# might appear across different rows of the same file.
MIXED_CASE_NULL_REGEX: Final[str] = r"(?i)^na$"
ROWS_WITH_MIXED_CASE_NULLS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "NA", "true", "key:value1"],
    ["3", "2023-02-20", "na", "false", "key:value2"],
    ["4", "2023-03-25", "Na", "true", "key:value3"],
]


def rows_with_null_placeholder(placeholder: str) -> list[list[str]]:
    """Build a small set of otherwise-valid rows using `placeholder` as the (to-be-replaced) float value.

    :param placeholder: literal text to insert into the `float` column of each row
    :return: a list of rows suitable for passing to `build_xsv_content`
    """
    return [
        ["2", "2023-01-15", placeholder, "true", "key:value1"],
        ["3", "2023-02-20", placeholder, "false", "key:value2"],
    ]


ROWS_WITH_NULLS: Final[list[list[str]]] = rows_with_null_placeholder(NULL_PLACEHOLDER)
# A single regex with alternation, matching several *different* null-value spellings (not just
# case variants of the same word) that might appear together in one file.
MULTI_SPELLING_NULL_REGEX: Final[str] = r"^n(a|ull|one)$"
ROWS_WITH_MULTIPLE_NULL_SPELLINGS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "na", "true", "key:value1"],
    ["3", "2023-02-20", "null", "false", "key:value2"],
    ["4", "2023-02-20", "none", "false", "key:value3"],
]

# Malformed schema content used to trigger qsv's own "bad config" exit codes (>1), as opposed to
# the module's own pre-flight JSON validation in `generate_header` (which is tested separately).
MALFORMED_SCHEMA_CASES: Final[list[tuple[str, str]]] = [
    ("truncated-json", '{"type": "object", "required": ['),
    ("not-json-at-all", "this is definitely not json"),
]

# Regexes with invalid syntax, used to trigger qsv replace's own "bad pattern" exit code (>1).
INVALID_NULL_REGEXES: Final[list[str]] = [
    "[unclosed",
    "(unclosed",
    "*leading-quantifier",
]


"""qsv_check"""


@pytest.mark.usefixtures("qsv_cmd")
def test_qsv_check_pass_qsv_found() -> None:
    """Qsv can be used for validation if the binary can be found and it outputs the expected stuff."""
    output = qsv_check()
    assert output


@pytest.mark.usefixtures("qsv_cmd")
def test_qsv_check_fail_qsv_version_unexpected_output(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that unexpected output from qsv version raises a runtime error."""
    original_subprocess_fn = subprocess.run

    def _patched(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str] | None:  # noqa: ANN003
        corrupted_cmd = [*cmd[:2], "--some-invalid-flag", *cmd[2:]]
        return original_subprocess_fn(corrupted_cmd, **kwargs)

    monkeypatch.setattr(subprocess, "run", _patched)

    err_regex = re.compile("Cannot perform validation with qsv: Command .*? returned non-zero exit status")
    with pytest.raises(RuntimeError, match=err_regex):
        qsv_check()


@pytest.mark.usefixtures("qsv_cmd")
def test_qsv_check_fail_qsv_version_throws_error(mock_qsv_run: Callable[..., MagicMock]) -> None:
    """Ensure that unexpected output from qsv version raises a runtime error."""
    mock_qsv_run(returncode=1, stderr="This wasn't supposed to happen")
    with pytest.raises(
        RuntimeError, match="`qsv --version` exited with code 1; STDERR: This wasn't supposed to happen"
    ):
        qsv_check()


def test_qsv_check_fail_qsv_not_found(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that running qsv_check when qsv is nowhere to be found results in a runtime error."""
    monkeypatch.setattr(shutil, "which", lambda _: None)
    monkeypatch.delenv("QSV_BIN", raising=False)

    with pytest.raises(RuntimeError, match="Could not locate the qsv binary"):
        qsv_check()


# oh crap, no QSV!
def test_run_qsv_input_fail_qsv_binary_not_found(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
) -> None:
    """Pin down current behaviour when qsv_cmd points at a nonexistent executable.

    subprocess.run raises FileNotFoundError from within _run_qsv_step, which is converted into
    an ErrorRecord.
    """
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(build_xsv_content(VALID_ROWS, header=COLUMNS), input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert "No such file or directory: 'mock-qsv'" in error.message


# run_qsv_input tests with mocked qsv
def test_run_qsv_input_mocked_qsv_pass_returns_cleaned_file_name(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """On a simulated exit code 0, the expected -cleaned file name is returned and populated."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_header_input_file(args, write_working_file)
    mock_qsv_run(returncode=0, output_content="mocked cleaned output\n")

    result = run_qsv_input(args, input_file_name)

    expected_output = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    assert result == expected_output
    assert (args.tmp_dir_path / expected_output).read_text() == "mocked cleaned output\n"
    assert args.errors == []


def test_run_qsv_input_mocked_qsv_pass_invokes_expected_command(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """The constructed qsv command and subprocess kwargs match what run_qsv_input is documented to build."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_header_input_file(args, write_working_file)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    run_qsv_input(args, input_file_name)

    output_file_name = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    mock.assert_called_once()
    (cmd,), kwargs = mock.call_args
    assert cmd == [
        FAKE_QSV_CMD,
        "input",
        "--trim-headers",
        "--trim-fields",
        "--encoding-errors",
        "strict",
        "--delimiter",
        args.delimiter,
        "--output",
        str(args.tmp_dir_path / output_file_name),
        str(args.tmp_dir_path / input_file_name),
    ]
    assert kwargs["text"] is True
    assert kwargs["env"] == args.qsv_env


@pytest.mark.parametrize(
    ("quote", "escape", "expected_extra_flags"),
    [
        pytest.param("'", None, ["--quote", "'"], id="quote-only"),
        pytest.param(None, "\\", ["--escape", "\\"], id="escape-only"),
        pytest.param("'", "\\", ["--quote", "'", "--escape", "\\"], id="quote-and-escape"),
    ],
)
def test_run_qsv_input_mocked_qsv_pass_invokes_expected_command_with_quote_and_escape(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    quote: str | None,
    escape: str | None,
    expected_extra_flags: list[str],
) -> None:
    """When quote and/or escape are set, --quote/--escape are appended (in that order) right after --delimiter."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    args = args.model_copy(update={"quote": quote, "escape": escape})
    input_file_name = _write_header_input_file(args, write_working_file)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    run_qsv_input(args, input_file_name)

    output_file_name = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    mock.assert_called_once()
    (cmd,), kwargs = mock.call_args
    assert cmd == [
        FAKE_QSV_CMD,
        "input",
        "--trim-headers",
        "--trim-fields",
        "--encoding-errors",
        "strict",
        "--delimiter",
        args.delimiter,
        *expected_extra_flags,
        "--output",
        str(args.tmp_dir_path / output_file_name),
        str(args.tmp_dir_path / input_file_name),
    ]
    assert kwargs["text"] is True
    assert kwargs["env"] == args.qsv_env


@pytest.mark.parametrize("delimiter", DELIMITERS)
def test_run_qsv_input_mocked_qsv_pass_passes_correct_delimiter(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    delimiter: str,
) -> None:
    """Every delimiter in SEP_TO_EXT is passed through to qsv verbatim via --delimiter."""
    source = write_source_file(
        build_xsv_content(VALID_ROWS, header=COLUMNS, delimiter=delimiter), f"data{SEP_TO_EXT[delimiter]}"
    )
    args = make_mock_args(source, delimiter=delimiter)
    input_file_name = _write_header_input_file(args, write_working_file)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    run_qsv_input(args, input_file_name)

    (cmd,), _ = mock.call_args
    assert cmd[cmd.index("--delimiter") + 1] == delimiter


def test_run_qsv_input_mocked_qsv_fail_records_error_without_invoking_qsv_when_input_file_missing(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """A missing input file records an error and the (mocked) subprocess.run is never called."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    result = run_qsv_input(args, "does-not-exist.tsv")

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.message.startswith(f"Input file not found at {args.tmp_dir_path / 'does-not-exist.tsv'!s}")

    mock.assert_not_called()


@pytest.mark.parametrize("returncode", [1, 2, 127, 255])
@pytest.mark.parametrize(
    "output_content",
    [None, "unexpected output written despite failure\n"],
    ids=["no-output-file-written", "output-file-written-anyway"],
)
def test_run_qsv_input_mocked_qsv_fail_records_error_and_returns_none_on_nonzero_exit(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    returncode: int,
    output_content: str | None,
) -> None:
    """Any non-zero exit code is records a single ErrorRecord and returns None."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_header_input_file(args, write_working_file)
    mock_qsv_run(returncode=returncode, stderr="simulated failure", output_content=output_content)

    result = run_qsv_input(args, input_file_name)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == returncode
    assert error.message == "simulated failure"


# run_qsv_input, real qsv!
def test_run_qsv_input_pass_returns_cleaned_file_name(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """On success, the -cleaned file is written to tmp_dir_path and its name is returned."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(build_xsv_content(VALID_ROWS, header=COLUMNS), input_file_name)

    result = run_qsv_input(args, input_file_name)

    expected_output = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    assert result == expected_output
    assert (args.tmp_dir_path / expected_output).is_file()
    assert args.errors == []


def test_run_qsv_input_pass_trims_header_and_field_whitespace(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """--trim-headers/--trim-fields strip surrounding whitespace from every field, including the header."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    padded_header = [f" {c} " for c in COLUMNS]
    write_working_file(build_xsv_content(ROWS_PADDED_WHITESPACE, header=padded_header), input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is not None
    rows = parse_xsv((args.tmp_dir_path / result).read_text())
    assert rows[0] == COLUMNS
    for row in rows[1:]:
        for field in row:
            assert field == field.strip()


@pytest.mark.parametrize("delimiter", DELIMITERS)
def test_run_qsv_input_pass_with_various_delimiters(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    delimiter: str,
) -> None:
    """run_qsv_input honours whichever delimiter is configured on args, for every delimiter SEP_TO_EXT knows about."""
    content = build_xsv_content(VALID_ROWS, header=COLUMNS, delimiter=delimiter)
    source = write_source_file(content, f"data{SEP_TO_EXT[delimiter]}")
    args = make_args(source, delimiter=delimiter)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(content, input_file_name)

    result = run_qsv_input(args, input_file_name)
    assert result is not None

    expected_output = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    assert result == expected_output
    rows = parse_xsv((args.tmp_dir_path / result).read_text(), delimiter=delimiter)
    assert rows[0] == COLUMNS


def test_run_qsv_input_pass_header_only_file(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """A file containing only a header row (no data rows) is valid input to `qsv input`."""
    source = write_source_file(build_xsv_content([], header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(build_xsv_content([], header=COLUMNS), input_file_name)

    result: str | None = run_qsv_input(args, input_file_name)
    assert result

    expected_output = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    assert result == expected_output
    rows = parse_xsv((args.tmp_dir_path / result).read_text())
    assert rows == [COLUMNS]
    assert args.errors == []


def test_run_qsv_input_pass_strips_comment_lines(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """Comment lines are removed from the output file."""
    clean_content = build_xsv_content(VALID_ROWS, header=COLUMNS)
    commented_content = interleave_comments(clean_content, comment_char="#")

    source = write_source_file(commented_content, "data.tsv")
    args = make_args(source, comment_char="#")
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(commented_content, input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is not None
    output_content = (args.tmp_dir_path / result).read_text()
    rows = parse_xsv(output_content)

    # no comment lines survive
    assert not any(line.startswith("#") for line in output_content.splitlines())
    # exactly the original, uncommented data remains, in the same order
    assert rows[0] == COLUMNS
    assert rows[1:] == VALID_ROWS
    assert args.errors == []


@pytest.mark.parametrize("comment_char", ["#", ";", "!", "%"])
def test_run_qsv_input_pass_strips_comment_lines_with_various_comment_chars(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    comment_char: str,
) -> None:
    """Any single-character comment_char configured on args is honoured via QSV_COMMENT_CHAR."""
    clean_content = build_xsv_content(VALID_ROWS, header=COLUMNS)
    commented_content = interleave_comments(clean_content, comment_char=comment_char)

    source = write_source_file(commented_content, "data.tsv")
    args = make_args(source, comment_char=comment_char)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(commented_content, input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is not None
    rows = parse_xsv((args.tmp_dir_path / result).read_text())
    assert rows[0] == COLUMNS
    assert rows[1:] == VALID_ROWS


def test_run_qsv_input_pass_does_not_strip_non_comment_lines_resembling_comment_char(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """Only lines where the comment_char is the very first character are stripped."""
    rows_with_hash_in_field = [
        ["2", "2023-01-15", "3.14", "true", "key:value#1"],
        ["3", "2023-02-20", "1.11", "false", "not#a-comment"],
        ["", "", "#", "", "#"],
    ]
    content = build_xsv_content(rows_with_hash_in_field, header=COLUMNS)

    source = write_source_file(content, "data.tsv")
    args = make_args(source, comment_char="#")
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(content, input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is not None
    rows = parse_xsv((args.tmp_dir_path / result).read_text())
    assert rows[1:] == rows_with_hash_in_field


def test_run_qsv_input_pass_custom_quote_character_preserves_delimiter_within_field(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
) -> None:
    """A custom --quote character lets a field safely contain the delimiter without the row being ragged."""
    quote_char = "@"
    raw_content = (
        "number\tdate\tfloat\tboolean\tstring\n"
        f"2\t2023-01-15\t3.14\ttrue\t{quote_char}key:value1\twith-embedded-tab{quote_char}\n"
    )
    source = write_source_file(raw_content, "data.tsv")
    args = make_args(source)
    args = args.model_copy(update={"quote": quote_char})
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(raw_content, input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is not None
    assert args.errors == []
    output_rows = list(csv.reader((args.tmp_dir_path / result).read_text().splitlines(), delimiter="\t"))
    assert output_rows[0] == COLUMNS
    assert len(output_rows) == 2  # noqa: PLR2004
    assert output_rows[1][:4] == ["2", "2023-01-15", "3.14", "true"]
    # the embedded tab survived as part of a single field, rather than splitting the row
    assert output_rows[1][4] == "key:value1\twith-embedded-tab"


def test_run_qsv_input_pass_custom_escape_character_unescapes_embedded_quotes(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
) -> None:
    """A custom --escape character is honoured when unescaping quote characters embedded within a quoted field."""
    escape_char = "\\"
    raw_content = 'number\tdate\tfloat\tboolean\tstring\n2\t2023-01-15\t3.14\ttrue\t"contains \\"escaped\\" quotes"\n'
    source = write_source_file(raw_content, "data.tsv")
    args = make_args(source)
    args = args.model_copy(update={"escape": escape_char})
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(raw_content, input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is not None
    assert args.errors == []
    output_rows = list(csv.reader((args.tmp_dir_path / result).read_text().splitlines(), delimiter="\t"))
    assert output_rows[0] == COLUMNS
    assert len(output_rows) == 2  # noqa: PLR2004
    assert output_rows[1][:4] == ["2", "2023-01-15", "3.14", "true"]
    assert output_rows[1][4] == 'contains "escaped" quotes'


# Failure states


def test_run_qsv_input_fail_records_error_when_input_file_missing(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """If input_file_name doesn't exist in tmp_dir_path, an error is recorded before qsv is even invoked."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)

    result = run_qsv_input(args, "does-not-exist.tsv")

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.message.startswith(f"Input file not found at {args.tmp_dir_path / 'does-not-exist.tsv'!s}")


@pytest.mark.parametrize("ragged_row_fixture", [PARTIAL_RAGGED_ROWS, ALL_RAGGED_ROWS])
def test_run_qsv_input_fail_ragged_records_error_and_returns_none(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    ragged_row_fixture: list[list[str]],
) -> None:
    """A file with ragged rows triggers qsv to exit with return code 1.

    Note that qsv streams to the output file, so a partial output file will be produced if it exits partway through.
    """
    source = write_source_file(build_xsv_content(ragged_row_fixture, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(build_xsv_content(ragged_row_fixture, header=COLUMNS), input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == 1
    assert error.message.startswith("Invalid CSV. Last valid row")

    # partial output file is produced
    output_file_name = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    assert (args.tmp_dir_path / output_file_name).exists()


@pytest.mark.parametrize("utf8_lines", [INVALID_UTF8_WITH_NULLS, INVALID_UTF8_NO_NULLS, INVALID_UTF8_MIXED])
def test_run_qsv_input_fail_records_error_and_returns_none_on_invalid_utf8(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    utf8_lines: bytes,
) -> None:
    """--encoding-errors strict causes real qsv input failure on malformed UTF-8; recorded as a soft ErrorRecord."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(utf8_lines, input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == 5  # noqa: PLR2004
    assert error.message.startswith("encoding error: STRICT. Invalid UTF8 - ")

    # partial output file is produced
    output_file_name = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    assert (args.tmp_dir_path / output_file_name).exists()


@pytest.mark.usefixtures("with_invalid_flag_injected")
def test_run_qsv_input_fail_records_error_on_bad_cli_flag(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
) -> None:
    """A malformed CLI invocation (bad flag) is recorded as an ErrorRecord, and returns None."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(build_xsv_content(VALID_ROWS, header=COLUMNS), input_file_name)

    result = run_qsv_input(args, input_file_name)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode != 0
    assert error.message.startswith("Unknown flag: '--some-invalid-flag'")


"""run_qsv_null_replacement"""


# run_qsv_null_replacement -- mocked qsv
def test_run_qsv_null_replacement_mocked_qsv_pass_returns_norm_file_name(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """On a simulated exit code 0, the expected -norm file name is returned and populated."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, VALID_ROWS)
    mock_qsv_run(returncode=0, output_content="mocked normalised output\n")

    result = run_qsv_null_replacement(args, input_file_name)

    expected_output = f"{args.xsv_file_base_name}{NORM_SUFFIX}{args.ext}"
    assert result == expected_output
    assert (args.tmp_dir_path / expected_output).read_text() == "mocked normalised output\n"
    assert args.errors == []


def test_run_qsv_null_replacement_mocked_qsv_pass_invokes_expected_command_without_regex_cols(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """When null_regex_cols is unset, no -s flag/columns are added to the command."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, VALID_ROWS)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    run_qsv_null_replacement(args, input_file_name)

    output_file_name = f"{args.xsv_file_base_name}{NORM_SUFFIX}{args.ext}"
    mock.assert_called_once()
    (cmd,), kwargs = mock.call_args
    assert cmd == [
        FAKE_QSV_CMD,
        "replace",
        "--not-one",
        "--delimiter",
        args.delimiter,
        "--output",
        str(args.tmp_dir_path / output_file_name),
        NULL_REGEX,
        "",
        str(args.tmp_dir_path / input_file_name),
    ]
    assert kwargs["text"] is True
    assert kwargs["env"] == args.qsv_env


def test_run_qsv_null_replacement_mocked_qsv_pass_invokes_expected_command_with_regex_cols(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """Ensure that the `cols` to be regexed are correctly conveyed to qsv."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    args = args.model_copy(update={"null_regex_cols": ["float", "string"]})
    input_file_name = _write_cleaned_input_file(args, write_working_file, VALID_ROWS)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    run_qsv_null_replacement(args, input_file_name)

    output_file_name = f"{args.xsv_file_base_name}{NORM_SUFFIX}{args.ext}"
    mock.assert_called_once()
    (cmd,), _ = mock.call_args
    assert cmd == [
        FAKE_QSV_CMD,
        "replace",
        "--not-one",
        "--delimiter",
        args.delimiter,
        "--output",
        str(args.tmp_dir_path / output_file_name),
        "-s",
        "float,string",
        NULL_REGEX,
        "",
        str(args.tmp_dir_path / input_file_name),
    ]


@pytest.mark.parametrize("delimiter", DELIMITERS)
def test_run_qsv_null_replacement_mocked_qsv_pass_passes_correct_delimiter(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    delimiter: str,
) -> None:
    """Every delimiter in SEP_TO_EXT is passed through to qsv verbatim via --delimiter."""
    source = write_source_file(
        build_xsv_content(VALID_ROWS, header=COLUMNS, delimiter=delimiter), f"data{SEP_TO_EXT[delimiter]}"
    )
    args = make_mock_args(source, delimiter=delimiter, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, VALID_ROWS, delimiter=delimiter)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    run_qsv_null_replacement(args, input_file_name)

    (cmd,), _ = mock.call_args
    assert cmd[cmd.index("--delimiter") + 1] == delimiter


@pytest.mark.parametrize("returncode", [1, 2, 127, 255])
@pytest.mark.parametrize(
    "output_content",
    [None, "unexpected output written despite failure\n"],
    ids=["no-output-file-written", "output-file-written-anyway"],
)
def test_run_qsv_null_replacement_mocked_qsv_fail_records_error_and_returns_none_on_nonzero_exit(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    returncode: int,
    output_content: str | None,
) -> None:
    """Any non-zero exit code is handled uniformly with an ErrorRecord and return None."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, VALID_ROWS)
    original_input_content = (args.tmp_dir_path / input_file_name).read_text()
    mock_qsv_run(returncode=returncode, stderr="simulated failure", output_content=output_content)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == returncode
    assert error.message == "simulated failure"

    copied_path = args.qsv_output_dir_path / input_file_name
    assert copied_path.is_file()
    assert copied_path.read_text() == original_input_content

    # the qsv-produced (ignored) output file, if any, must not be what ended up in qsv_output_dir_path
    norm_output_name = f"{args.xsv_file_base_name}{NORM_SUFFIX}{args.ext}"
    if output_content is not None:
        assert not (args.qsv_output_dir_path / norm_output_name).exists()


def test_run_qsv_null_replacement_mocked_qsv_fail_records_both_qsv_and_copy_errors_when_copy_also_fails(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """A copy failure after a qsv failure generates a second error."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, VALID_ROWS)
    mock_qsv_run(returncode=1, stderr="simulated qsv failure", output_content=None)

    # remove the output directory out from under run_qsv_null_replacement so the subsequent
    # copy_safely call fails with a real OSError
    args.qsv_output_dir_path.rmdir()

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is None
    assert len(args.errors) == 2

    qsv_error, copy_error = args.errors
    assert qsv_error.file == args.file_name
    assert qsv_error.returncode == 1
    assert qsv_error.message == "simulated qsv failure"

    assert copy_error.file == args.file_name
    assert copy_error.returncode is None
    assert "No such file or directory" in copy_error.message


def test_run_qsv_null_replacement_mocked_qsv_fail_records_error_without_invoking_qsv_when_input_file_missing(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """A missing input file records an error and the mocked subprocess.run is never called."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    mock = mock_qsv_run(returncode=0, output_content="content\n")

    result = run_qsv_null_replacement(args, "does-not-exist.tsv")

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.message.startswith(f"Input file not found at {args.tmp_dir_path / 'does-not-exist.tsv'!s}")

    mock.assert_not_called()
    assert list(args.qsv_output_dir_path.iterdir()) == []


# real qsv!
def test_run_qsv_null_replacement_pass_replaces_matching_placeholder(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """On success, every exact match of null_regex is replaced with an empty string."""
    source = write_source_file(build_xsv_content(ROWS_WITH_NULLS, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, ROWS_WITH_NULLS)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    expected_output = f"{args.xsv_file_base_name}{NORM_SUFFIX}{args.ext}"
    assert result == expected_output
    rows = parse_xsv((args.tmp_dir_path / result).read_text())
    assert rows[0] == COLUMNS
    float_col_index = COLUMNS.index("float")
    for row in rows[1:]:
        assert row[float_col_index] == ""
    assert args.errors == []


@pytest.mark.parametrize("case", NULL_REGEX_CASES, ids=[c.placeholder for c in NULL_REGEX_CASES])
def test_run_qsv_null_replacement_pass_matches_various_placeholder_regex_pairs(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    case: NullValueCase,
) -> None:
    """Each documented (placeholder, regex) pair -- including one requiring escaping -- is matched exactly."""
    rows = rows_with_null_placeholder(case.placeholder)
    source = write_source_file(build_xsv_content(rows, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=case.regex)
    input_file_name = _write_cleaned_input_file(args, write_working_file, rows)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    parsed = parse_xsv((args.tmp_dir_path / result).read_text())
    float_col_index = COLUMNS.index("float")
    for row in parsed[1:]:
        assert row[float_col_index] == ""


def test_run_qsv_null_replacement_pass_case_insensitive_regex_matches_all_spellings(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """A single case-insensitive regex replaces differently-cased spellings across different rows."""
    source = write_source_file(build_xsv_content(ROWS_WITH_MIXED_CASE_NULLS, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=MIXED_CASE_NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, ROWS_WITH_MIXED_CASE_NULLS)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    parsed = parse_xsv((args.tmp_dir_path / result).read_text())
    float_col_index = COLUMNS.index("float")
    for row in parsed[1:]:
        assert row[float_col_index] == ""


def test_run_qsv_null_replacement_pass_alternation_regex_matches_multiple_spellings(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """A single regex with alternation replaces several distinct null-value spellings in one pass."""
    source = write_source_file(build_xsv_content(ROWS_WITH_MULTIPLE_NULL_SPELLINGS, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=MULTI_SPELLING_NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, ROWS_WITH_MULTIPLE_NULL_SPELLINGS)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    parsed = parse_xsv((args.tmp_dir_path / result).read_text())
    float_col_index = COLUMNS.index("float")
    for row in parsed[1:]:
        assert row[float_col_index] == ""


@pytest.mark.parametrize("cols", [["float"], ["float", "string"]])
def test_run_qsv_null_replacement_pass_restricts_replacement_to_selected_columns(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    cols: list[str],
) -> None:
    """When null_regex_cols is set, only the named column(s) are affected."""
    rows = [
        ["2", "2023-01-15", "NA", "true", "NA"],
        ["3", "2023-02-20", "NA", "false", "NA"],
    ]
    source = write_source_file(build_xsv_content(rows, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=NULL_REGEX)
    args = args.model_copy(update={"null_regex_cols": cols})
    input_file_name = _write_cleaned_input_file(args, write_working_file, rows)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    parsed = parse_xsv((args.tmp_dir_path / result).read_text())
    float_col_index = COLUMNS.index("float")
    string_col_index = COLUMNS.index("string")
    for row in parsed[1:]:
        assert row[float_col_index] == ""
        if "string" in cols:
            assert row[string_col_index] == ""
        else:
            assert row[string_col_index] == "NA"  # untouched: not in null_regex_cols


@pytest.mark.parametrize("delimiter", DELIMITERS)
def test_run_qsv_null_replacement_pass_with_various_delimiters(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    delimiter: str,
) -> None:
    """Every delimiter in SEP_TO_EXT is honoured, both for parsing input and writing output."""
    source = write_source_file(
        build_xsv_content(ROWS_WITH_NULLS, header=COLUMNS, delimiter=delimiter), f"data{SEP_TO_EXT[delimiter]}"
    )
    args = make_args(source, delimiter=delimiter, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, ROWS_WITH_NULLS, delimiter=delimiter)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    rows = parse_xsv((args.tmp_dir_path / result).read_text(), delimiter=delimiter)
    assert rows[0] == COLUMNS


def test_run_qsv_null_replacement_pass_returns_output_unchanged_when_no_matches_found(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """`--not-one` means qsv still exits 0 (success) even when the regex matches nothing at all."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=r"^this-will-never-match$")
    input_file_name = _write_cleaned_input_file(args, write_working_file, VALID_ROWS)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    assert args.errors == []
    rows = parse_xsv((args.tmp_dir_path / result).read_text())
    assert rows[1:] == VALID_ROWS


@pytest.mark.parametrize("utf8_lines", [INVALID_UTF8_NO_NULLS, INVALID_UTF8_WITH_NULLS, INVALID_UTF8_MIXED])
def test_run_qsv_null_replacement_pass_error_on_invalid_utf8_with_matches(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    utf8_lines: bytes,
) -> None:
    """Malformed (non-UTF-8) input bytes cause qsv replace to fail."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=NULL_REGEX)
    input_file_name = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    write_working_file(utf8_lines, input_file_name)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is not None
    assert args.errors == []

    # output file will contain unicode errors
    with pytest.raises(UnicodeDecodeError, match="codec can't decode byte"):
        (args.tmp_dir_path / result).read_text()

    if utf8_lines == INVALID_UTF8_NO_NULLS:
        # output file is the same as the input file
        assert (args.tmp_dir_path / result).read_bytes() == (args.tmp_dir_path / input_file_name).read_bytes()
    else:
        # output file has been changed
        assert (args.tmp_dir_path / result).read_bytes() != (args.tmp_dir_path / input_file_name).read_bytes()


def test_run_qsv_null_replacement_fail_records_error_when_input_file_missing(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """A missing input file records an error before qsv is ever invoked."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=NULL_REGEX)

    result = run_qsv_null_replacement(args, "does-not-exist.tsv")
    assert result is None

    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.message.startswith(f"Input file not found at {args.tmp_dir_path / 'does-not-exist.tsv'!s}")


@pytest.mark.parametrize("ragged_row_fixture", [PARTIAL_RAGGED_ROWS, ALL_RAGGED_ROWS])
def test_run_qsv_null_replacement_fail_ragged_rows_records_error_and_returns_none(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    ragged_row_fixture: list[list[str]],
) -> None:
    """When rows are ragged, qsv's CSV reader fails."""
    source = write_source_file(build_xsv_content(ragged_row_fixture, header=COLUMNS), "data.tsv")
    args = make_args(source, null_regex=NULL_REGEX)
    input_file_name = _write_cleaned_input_file(args, write_working_file, ragged_row_fixture)

    result = run_qsv_null_replacement(args, input_file_name)

    assert result is None
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode != 0
    assert error.message.startswith("csv error: CSV error: ")

    copied_path = args.qsv_output_dir_path / input_file_name
    assert copied_path.is_file()
    assert copied_path.read_text() == (args.tmp_dir_path / input_file_name).read_text()
    # confirm the -norm output name specifically was NOT what got copied
    norm_output_name = f"{args.xsv_file_base_name}{NORM_SUFFIX}{args.ext}"
    assert not (args.qsv_output_dir_path / norm_output_name).exists()


"""run_qsv_validate"""


# run_qsv_validate, qsv mocked
def test_run_qsv_validate_mocked_qsv_pass_returns_valid_output_when_produced(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """--valid-output returns a 1 error code and the valid output file path is returned."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)

    mock_qsv_run(returncode=1, output_content="mocked all-valid output\n", output_flag="--valid-output")

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result == files.valid_output.name
    assert files.valid_output.read_text() == "mocked all-valid output\n"
    assert args.errors == []


def test_run_qsv_validate_mocked_qsv_pass_invokes_expected_command_first_pass(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """First-pass invocations include --no-format-validation and every documented flag."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)
    mock = mock_qsv_run(returncode=1, output_content="content\n", output_flag="--valid-output")

    run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    mock.assert_called_once()
    (cmd,), kwargs = mock.call_args
    assert cmd == [
        FAKE_QSV_CMD,
        "validate",
        "--delimiter",
        args.delimiter,
        "--no-format-validation",
        "--split-ragged",
        "--trim",
        "--valid",
        args.valid_file_suffix,
        "--invalid",
        args.invalid_file_suffix,
        "--valid-output",
        str(files.valid_output),
        str(args.tmp_dir_path / input_file_name),
        str(files.schema),
    ]
    assert None not in cmd
    assert kwargs["text"] is True
    assert kwargs["env"] == args.qsv_env


def test_run_qsv_validate_mocked_qsv_pass_invokes_expected_command_second_pass_omits_no_format_validation(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """Second-pass invocations omit --no-format-validation entirely."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, NORM_SUFFIX)
    mock = mock_qsv_run(returncode=1, output_content="content\n", output_flag="--valid-output")

    run_qsv_validate(args, input_file_name, args.post_norm_schema, first_pass=False)

    (cmd,), _ = mock.call_args
    assert cmd == [
        FAKE_QSV_CMD,
        "validate",
        "--delimiter",
        args.delimiter,
        "--split-ragged",
        "--trim",
        "--valid",
        args.valid_file_suffix,
        "--invalid",
        args.invalid_file_suffix,
        "--valid-output",
        str(args.tmp_dir_path / f"data-normalised-validated{VALID_SUFFIX}.tsv"),
        str(args.tmp_dir_path / input_file_name),
        str(args.tmp_dir_path / "data-normalised-validated.jsonschema.json"),
    ]


@pytest.mark.parametrize("delimiter", DELIMITERS)
def test_run_qsv_validate_mocked_qsv_pass_passes_correct_delimiter(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    delimiter: str,
) -> None:
    """Every delimiter in SEP_TO_EXT is passed through to qsv verbatim via --delimiter."""
    source = write_source_file(
        build_xsv_content(VALID_ROWS, header=COLUMNS, delimiter=delimiter), f"data{SEP_TO_EXT[delimiter]}"
    )
    args = make_mock_args(source, delimiter=delimiter)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    mock = mock_qsv_run(returncode=1, output_content="content\n", output_flag="--valid-output")

    run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    (cmd,), _ = mock.call_args
    assert cmd[cmd.index("--delimiter") + 1] == delimiter


def test_run_qsv_validate_mocked_qsv_pass_recovers_valid_lines(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """Valid lines can be recovered if there are a non-zero number of non-header lines in the -valid file.

    When --valid-output isn't produced but qsv's default-named .valid file contains recoverable
    (non-header) rows, run_qsv_validate records the error, copies all present files to
    qsv_output_dir_path, and moves the valid subset into the --valid-output path.
    """
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    mock_qsv_run(
        returncode=1,
        stderr="1 out of 3 records invalid",
        extra_files={
            valid_lines_file: build_xsv_content(VALID_ROWS[:1], header=COLUMNS),
            invalid_lines_file: build_xsv_content(VALID_ROWS[1:], header=COLUMNS),
            errors_file: "row 2: field 'float': simulated error\n",
        },
    )

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result == output.name
    assert output.read_text() == build_xsv_content(VALID_ROWS[:1], header=COLUMNS)
    # the qsv-named valid file was moved (not copied) within tmp_dir_path, so it no longer exists there
    assert not valid_lines_file.exists()

    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == 1
    assert error.message == "1 out of 3 records invalid"

    # all three qsv-produced files should have been copied to qsv_output_dir_path under their own names
    assert (args.qsv_output_dir_path / valid_lines_file.name).read_text() == build_xsv_content(
        VALID_ROWS[:1], header=COLUMNS
    )
    assert (args.qsv_output_dir_path / invalid_lines_file.name).is_file()
    assert (args.qsv_output_dir_path / errors_file.name).is_file()


def test_run_qsv_validate_mocked_qsv_fail_returns_none_when_valid_file_only_contains_header(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """If the valid lines file has no non-header lines, run_qsv_validate returns None."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    mock_qsv_run(
        returncode=1,
        stderr="all records invalid",
        extra_files={
            valid_lines_file: build_xsv_content([], header=COLUMNS),  # header only
            invalid_lines_file: build_xsv_content(VALID_ROWS, header=COLUMNS),
            errors_file: "all rows failed\n",
        },
    )

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result is None
    assert not output.exists()
    # since it was only copied, not moved, the original (header-only) file remains in tmp_dir_path
    assert valid_lines_file.exists()
    assert not (args.qsv_output_dir_path / valid_lines_file.name).is_file()
    assert len(args.errors) == 1


def test_run_qsv_validate_mocked_qsv_fail_returns_none_when_no_output_files_produced_at_all(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
) -> None:
    """Simulate a total failure, with no output files produced.

    A total failure with no output files at all (e.g. a config-level problem) still results in
    a single recorded ErrorRecord and a None return, with nothing to copy to qsv_output_dir_path.
    """
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    mock_qsv_run(returncode=2, stderr="bad schema")

    output_dir_before = snapshot_dir(args.qsv_output_dir_path)
    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result is None
    assert snapshot_dir(args.qsv_output_dir_path) == output_dir_before
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.returncode == 2
    assert error.message == "bad schema"


def test_run_qsv_validate_mocked_qsv_fail_returns_none_when_move_of_valid_lines_fails(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If the valid_lines_file cannot be successfully renamed, run_qsv_validate returns None."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    mock_qsv_run(
        returncode=1,
        stderr="mixed results",
        extra_files={
            valid_lines_file: build_xsv_content(VALID_ROWS, header=COLUMNS),
            invalid_lines_file: build_xsv_content([], header=COLUMNS),
            errors_file: "some rows failed\n",
        },
    )
    monkeypatch.setattr(qsv, "move_safely", lambda *_a, **_k: False)

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result is None
    # no errors to check as the whole function was bypassed


def test_run_qsv_validate_mocked_qsv_fail_records_error_without_invoking_qsv_when_input_file_missing(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    mock_qsv_run: Callable[..., MagicMock],
    first_pass_validated_schema: ValidatedSchema,
) -> None:
    """A missing input file records an error and the mocked subprocess.run is never called."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    mock = mock_qsv_run(returncode=0, output_content="content\n", output_flag="--valid-output")

    result = run_qsv_validate(args, "does-not-exist.tsv", first_pass_validated_schema, first_pass=True)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.message.startswith(f"Input file not found at {args.tmp_dir_path / 'does-not-exist.tsv'!s}")

    mock.assert_not_called()


# qsv validate, real qsv binary used
def test_run_qsv_validate_pass_all_valid_first_pass_returns_valid_output_directly(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """All rows valid: qsv writes to --valid-output, return value is the --valid-output file name.

    When every row is structurally valid, qsv writes directly to --valid-output; run_qsv_validate
    returns that file name immediately, without recording an error, regardless of qsv's exit code.
    """
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)
    assert result is not None

    assert result == output.name
    assert args.tmp_dir_path / result == output
    rows = parse_xsv((args.tmp_dir_path / result).read_text())
    assert rows[0] == COLUMNS
    assert rows[1:] == VALID_ROWS
    assert args.errors == []


def test_run_qsv_validate_pass_all_valid_second_pass_returns_valid_output_directly(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """Same all-valid behaviour holds for the stricter, format-checking second pass."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, VALID_ROWS, NORM_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=False)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    result = run_qsv_validate(args, input_file_name, args.post_norm_schema, first_pass=False)
    assert result == output.name
    assert args.errors == []


def test_run_qsv_validate_pass_recovers_valid_lines_from_partially_ragged_input(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """Run qsv validate over a file with some ragged rows.

    With a mix of valid and ragged rows, --split-ragged routes the bad row to .invalid while
    the good rows land in the qsv-default-named .valid file; run_qsv_validate records the resulting
    non-zero-exit-code error, copies all produced files to qsv_output_dir_path, and promotes the
    recovered valid subset to the --valid-output path, returning its name.
    """
    source = write_source_file(build_xsv_content(PARTIAL_RAGGED_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, PARTIAL_RAGGED_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result == output.name
    assert output.is_file()
    assert non_header_lines_present(output)

    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == 1
    assert error.message == "1 out of 3 records invalid."

    # the qsv-default-named valid/invalid/errors files should have been copied to qsv_output_dir_path
    for f in [errors_file, valid_lines_file, invalid_lines_file]:
        assert (args.qsv_output_dir_path / f.name).is_file()
        if f == valid_lines_file:
            assert not f.exists()
        else:
            assert f.is_file()


def test_run_qsv_validate_fail_all_rows_ragged_returns_none(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """Run qsv validate over a file where all rows are ragged.

    When every row is ragged, the qsv-default .valid file ends up containing only the header
    (no recoverable data), so run_qsv_validate records the error but returns None -- there is nothing
    worth promoting to a --valid-output file.
    """
    source = write_source_file(build_xsv_content(ALL_RAGGED_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, ALL_RAGGED_ROWS, HEADER_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=True)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result is None
    assert not output.is_file()
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == 1
    assert error.message == "2 out of 2 records invalid."

    # the qsv-default-named valid/invalid/errors files should have been copied to qsv_output_dir_path
    for f in [errors_file, valid_lines_file, invalid_lines_file]:
        assert f.is_file()
        if f == valid_lines_file:
            assert not (args.qsv_output_dir_path / f.name).exists()
        else:
            assert (args.qsv_output_dir_path / f.name).is_file()


@pytest.mark.parametrize(
    "invalid_rows",
    [
        pytest.param(ROWS_INVALID_NUMBER, id="invalid-number"),
        pytest.param(ROWS_INVALID_DATE, id="invalid-date"),
        pytest.param(ROWS_INVALID_STRING, id="invalid-string"),
    ],
)
def test_run_qsv_validate_fail_second_pass_catches_field_format_errors(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    invalid_rows: list[list[str]],
) -> None:
    """Ensure that field format checks are executed when first_pass=False.

    Field-format checks (number/date/string pattern) are only enforced on the second pass
    (--no-format-validation is omitted); rows that would structurally pass the loose first-pass
    schema are correctly rejected here, with every row invalid and no valid rows to recover.
    """
    source = write_source_file(build_xsv_content(invalid_rows, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, invalid_rows, NORM_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=False)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    result = run_qsv_validate(args, input_file_name, args.post_norm_schema, first_pass=False)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == 1
    assert error.message == "2 out of 2 records invalid."

    # errors and invalid output files produced
    # no valid lines file as all lines were invalid
    assert set(args.qsv_output_dir_path.glob("*.tsv")) == {
        args.qsv_output_dir_path / f.name for f in [errors_file, invalid_lines_file]
    }
    assert valid_lines_file.exists()
    assert not non_header_lines_present(valid_lines_file)
    assert not output.exists()


def test_run_qsv_validate_pass_second_pass_recovers_valid_rows_mixed_with_invalid(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, write_working_file: WriteFile
) -> None:
    """Second pass validation: valid and invalid rows; files copied to output dir.

    A mix of valid and format-invalid rows on the second pass still recovers the valid subset,
    exactly as the first-pass ragged-row recovery case does.
    """
    mixed_rows = [*VALID_ROWS, *ROWS_INVALID_DATE]
    source = write_source_file(build_xsv_content(mixed_rows, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = _write_validator_input_file(args, write_working_file, mixed_rows, NORM_SUFFIX)
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=False)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    result = run_qsv_validate(args, input_file_name, args.post_norm_schema, first_pass=False)

    assert result == output.name
    recovered_rows = parse_xsv(output.read_text())
    assert recovered_rows[1:] == VALID_ROWS
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode == 1
    assert error.message == "2 out of 4 records invalid."

    # valid, errors, and invalid output files produced
    assert set(args.qsv_output_dir_path.glob("*.tsv")) == {
        args.qsv_output_dir_path / f.name for f in [errors_file, valid_lines_file, invalid_lines_file]
    }


@pytest.mark.parametrize("utf8_lines", [INVALID_UTF8_NO_NULLS, INVALID_UTF8_WITH_NULLS, INVALID_UTF8_MIXED])
def test_run_qsv_validate_fail_records_error_on_invalid_utf8_input(
    make_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    write_working_file: WriteFile,
    utf8_lines: bytes,
) -> None:
    """Invalid utf-8 ."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(utf8_lines, input_file_name)

    result = run_qsv_validate(args, input_file_name, args.first_pass_schema, first_pass=True)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode != 0
    assert error.message.startswith("encoding error: CSV header is not valid UTF-8")


"""Tests for clean_validate_file - the big cheese of the qsv functions!"""


# clean_validate_file, orchestration mocked
def test_clean_validate_file_raises_validation_error_for_invalid_args_type() -> None:
    """@validate_call enforces that clean_validate_file's sole argument is a CleanerValidatorArgs instance."""
    with pytest.raises(ValidationError, match="Input should be a valid dictionary or instance of CleanerValidatorArgs"):
        clean_validate_file(SimpleNamespace(errors=[]))  # type: ignore[arg-type]


def test_clean_validate_file_missing_header_true_prepends_header_and_calls_first_pass_validate_with_expected_arguments(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """missing_header=True prepends the header file to the source data before the first validation pass."""
    source = write_source_file(build_xsv_content(VALID_ROWS), "data.tsv")
    args = make_mock_args(source, missing_header=True)
    mock_validate = MagicMock(return_value=None)
    monkeypatch.setattr(qsv, "run_qsv_validate", mock_validate)

    result = clean_validate_file(args)

    assert result is None
    file_with_headers = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    written = (args.tmp_dir_path / file_with_headers).read_text()
    assert written == args.header_file_path.read_text() + build_xsv_content(VALID_ROWS)

    mock_validate.assert_called_once()
    call = mock_validate.call_args
    assert call.args == (args, file_with_headers)
    assert call.kwargs == {"schema": args.first_pass_schema, "first_pass": True}


def test_clean_validate_file_missing_header_false_copies_file_and_calls_first_pass_validate_with_expected_arguments(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """missing_header=False (the default) just copies the source file into tmp_dir_path under the -header name."""
    content = build_xsv_content(VALID_ROWS, header=COLUMNS)
    source = write_source_file(content, "data.tsv")
    args = make_mock_args(source)
    mock_validate = MagicMock(return_value=None)
    monkeypatch.setattr(qsv, "run_qsv_validate", mock_validate)

    result = clean_validate_file(args)

    assert result is None
    file_with_headers = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    assert (args.tmp_dir_path / file_with_headers).read_text() == content

    mock_validate.assert_called_once()
    call = mock_validate.call_args
    assert call.args == (args, file_with_headers)
    assert call.kwargs == {"schema": args.first_pass_schema, "first_pass": True}


def test_clean_validate_file_fail_missing_header_prepend_oserror_records_error_and_returns_none(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An OSError while prepending the header (e.g. the source file vanishing) is recorded and short-circuits."""
    source = write_source_file(build_xsv_content(VALID_ROWS), "data.tsv")
    args = make_mock_args(source, missing_header=True)
    mock_validate = MagicMock()
    monkeypatch.setattr(qsv, "run_qsv_validate", mock_validate)
    source.unlink()

    result = clean_validate_file(args)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert "No such file or directory" in error.message
    mock_validate.assert_not_called()


def test_clean_validate_file_fail_copy_safely_failure_records_error_and_returns_none(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When the header is already present, a copy_safely failure short-circuits before any validation runs."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    mock_validate = MagicMock()
    monkeypatch.setattr(qsv, "run_qsv_validate", mock_validate)
    source.unlink()

    result = clean_validate_file(args)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert "No such file or directory" in error.message
    mock_validate.assert_not_called()


def test_clean_validate_file_fail_first_pass_validate_returns_none_short_circuits(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If the first-pass validation step fails outright, run_qsv_input is never invoked."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    monkeypatch.setattr(qsv, "run_qsv_validate", MagicMock(return_value=None))
    mock_input = MagicMock()
    monkeypatch.setattr(qsv, "run_qsv_input", mock_input)

    result = clean_validate_file(args)

    assert result is None
    mock_input.assert_not_called()


def test_clean_validate_file_fail_run_qsv_input_returns_none_short_circuits(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If run_qsv_input fails, neither null replacement nor the second validation pass runs."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    monkeypatch.setattr(qsv, "run_qsv_validate", MagicMock(return_value="first-pass-output.tsv"))
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value=None))
    mock_null_replace = MagicMock()
    monkeypatch.setattr(qsv, "run_qsv_null_replacement", mock_null_replace)

    result = clean_validate_file(args)

    assert result is None
    mock_null_replace.assert_not_called()


def test_clean_validate_file_pass_skips_null_replacement_when_null_regex_unset(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When null_regex is unset, run_qsv_null_replacement is never called.

    The second validate pass is invoked directly on run_qsv_input's output file.
    """
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    mock_validate = MagicMock(side_effect=["first-pass-output.tsv", "second-pass-output.tsv"])
    monkeypatch.setattr(qsv, "run_qsv_validate", mock_validate)
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value="cleaned.tsv"))
    mock_null_replace = MagicMock()
    monkeypatch.setattr(qsv, "run_qsv_null_replacement", mock_null_replace)
    monkeypatch.setattr(qsv, "move_safely", MagicMock(return_value=True))

    result = clean_validate_file(args)

    assert result == "second-pass-output.tsv"
    mock_null_replace.assert_not_called()
    assert mock_validate.call_count == 2  # noqa: PLR2004
    second_call = mock_validate.call_args_list[1]
    assert second_call.args == (args, "cleaned.tsv", args.post_norm_schema)
    assert second_call.kwargs == {"first_pass": False}


def test_clean_validate_file_pass_runs_null_replacement_when_null_regex_set_and_feeds_its_output_to_second_validate(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When null_regex is set, run_qsv_null_replacement runs and its output feeds the second validate pass."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    mock_validate = MagicMock(side_effect=["first-pass-output.tsv", "final-output.tsv"])
    monkeypatch.setattr(qsv, "run_qsv_validate", mock_validate)
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value="cleaned.tsv"))
    mock_null_replace = MagicMock(return_value="normalised.tsv")
    monkeypatch.setattr(qsv, "run_qsv_null_replacement", mock_null_replace)
    monkeypatch.setattr(qsv, "move_safely", MagicMock(return_value=True))

    result = clean_validate_file(args)

    assert result == "final-output.tsv"
    mock_null_replace.assert_called_once_with(args, "cleaned.tsv")
    second_call = mock_validate.call_args_list[1]
    assert second_call.args == (args, "normalised.tsv", args.post_norm_schema)


def test_clean_validate_file_fail_null_replacement_returns_none_short_circuits(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If null replacement fails, the second (post-normalisation) validation pass never runs."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, null_regex=NULL_REGEX)
    mock_validate = MagicMock(return_value="first-pass-output.tsv")
    monkeypatch.setattr(qsv, "run_qsv_validate", mock_validate)
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value="cleaned.tsv"))
    monkeypatch.setattr(qsv, "run_qsv_null_replacement", MagicMock(return_value=None))

    result = clean_validate_file(args)

    assert result is None
    mock_validate.assert_called_once()


def test_clean_validate_file_fail_second_pass_validate_returns_none_skips_move_and_returns_none(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """If the second (post-normalisation) validation pass fails outright, nothing is moved and None is returned."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    monkeypatch.setattr(qsv, "run_qsv_validate", MagicMock(side_effect=["first-pass-output.tsv", None]))
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value="cleaned.tsv"))
    mock_move = MagicMock()
    monkeypatch.setattr(qsv, "move_safely", mock_move)

    result = clean_validate_file(args)

    assert result is None
    mock_move.assert_not_called()


def test_clean_validate_file_pass_moves_validated_output_to_validated_file_dir_path(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """On overall success, the validated file is moved from tmp_dir_path to validated_file_dir_path."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    monkeypatch.setattr(qsv, "run_qsv_validate", MagicMock(side_effect=["first-pass-output.tsv", "final-output.tsv"]))
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value="cleaned.tsv"))
    mock_move = MagicMock(return_value=True)
    monkeypatch.setattr(qsv, "move_safely", mock_move)

    result = clean_validate_file(args)

    assert result == "final-output.tsv"
    mock_move.assert_called_once_with(
        args.tmp_dir_path / "final-output.tsv", args.validated_file_dir_path / "final-output.tsv", args
    )


def test_clean_validate_file_fail_move_safely_failure_returns_none_despite_valid_output(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Even when both validation passes succeed, a failed move_safely means None is returned, not the file name."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    monkeypatch.setattr(qsv, "run_qsv_validate", MagicMock(side_effect=["first-pass-output.tsv", "final-output.tsv"]))
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value="cleaned.tsv"))
    mock_move = MagicMock(return_value=False)
    monkeypatch.setattr(qsv, "move_safely", mock_move)

    result = clean_validate_file(args)

    assert result is None
    mock_move.assert_called_once_with(
        args.tmp_dir_path / "final-output.tsv", args.validated_file_dir_path / "final-output.tsv", args
    )


def test_clean_validate_file_pass_move_safely_success_returns_validated_file_name(
    make_mock_args: Callable[..., CleanerValidatorArgs],
    write_source_file: WriteFile,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Conversely, when move_safely succeeds, the validated file name is returned as before."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)
    monkeypatch.setattr(qsv, "run_qsv_validate", MagicMock(side_effect=["first-pass-output.tsv", "final-output.tsv"]))
    monkeypatch.setattr(qsv, "run_qsv_input", MagicMock(return_value="cleaned.tsv"))
    mock_move = MagicMock(return_value=True)
    monkeypatch.setattr(qsv, "move_safely", mock_move)

    result = clean_validate_file(args)

    assert result == "final-output.tsv"
    mock_move.assert_called_once_with(
        args.tmp_dir_path / "final-output.tsv", args.validated_file_dir_path / "final-output.tsv", args
    )


# clean_validate_file, real qsv end-to-end
def test_clean_validate_file_pass_end_to_end_with_header_present(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """A well-formed file with a header already present passes through the whole pipeline unmolested."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)

    result = clean_validate_file(args)

    files = generate_qsv_validate_file_names(args, "ignored", first_pass=False)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    assert result == output.name

    validated_path = args.validated_file_dir_path / output.name
    assert validated_path.is_file()
    rows = parse_xsv(validated_path.read_text())
    assert rows[0] == COLUMNS
    assert rows[1:] == VALID_ROWS
    assert args.errors == []


def test_clean_validate_file_pass_end_to_end_with_missing_header(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """A headerless file has the header prepended before running through the rest of the pipeline."""
    source = write_source_file(build_xsv_content(VALID_ROWS), "data.tsv")
    args = make_args(source, missing_header=True)

    result = clean_validate_file(args)

    files = generate_qsv_validate_file_names(args, "ignored", first_pass=False)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema
    assert result == output.name

    validated_path = args.validated_file_dir_path / output.name
    assert validated_path.is_file()
    rows = parse_xsv(validated_path.read_text())
    assert rows[0] == COLUMNS
    assert rows[1:] == VALID_ROWS
    assert args.errors == []


def test_clean_validate_file_fail_end_to_end_move_failure_returns_none(
    make_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, monkeypatch: pytest.MonkeyPatch
) -> None:
    """clean_validate_file returns None if the final move to validated_file_dir_path fails."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_args(source)

    original_move_safely = qsv.move_safely

    def selective_move(src: Path, dst: Path, args: CleanerValidatorArgs) -> bool:
        """Throw an error if the move is to the validated_file_dir_path."""
        if src.name == "data-normalised-validated-valid.tsv":
            return original_move_safely(src, Path("/some/random/place"), args)
        return original_move_safely(src, dst, args)

    monkeypatch.setattr(qsv, "move_safely", selective_move)

    result = clean_validate_file(args)

    assert result is None
    assert len(args.errors) == 1
    error = args.errors[0]
    assert error.file == args.file_name
    assert error.returncode is None
    assert "No such file or directory: '/some/random/place'" in error.message

    # the validated content still exists in tmp_dir_path since the move never completed
    files = generate_qsv_validate_file_names(args, "ignored", first_pass=False)
    output = files.valid_output
    errors_file = files.errors
    valid_lines_file = files.valid_lines
    invalid_lines_file = files.invalid_lines
    schema_file = files.schema

    assert (args.tmp_dir_path / output.name).is_file()
