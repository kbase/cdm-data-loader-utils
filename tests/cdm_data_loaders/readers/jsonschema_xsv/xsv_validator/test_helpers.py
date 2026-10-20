"""Tests for the xsv validator helper (non-qsv-interacting) code."""

import os
import subprocess
from collections.abc import Callable
from pathlib import Path
from types import SimpleNamespace
from typing import Final
from unittest.mock import MagicMock

import pytest
from _pytest.mark.structures import ParameterSet

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.helpers import (
    SEP_TO_EXT,
    VALID_SUFFIX,
    VALIDATION_ERRORS,
    CleanerValidatorArgs,
    ErrorRecord,
    copy_safely,
    generate_qsv_validate_file_names,
    move_safely,
    non_header_lines_present,
    prepend_header,
)
from tests.cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.conftest import (
    COLUMNS,
    DELIMITERS,
    VALID_ROWS,
    WriteFile,
    _touch,
    build_xsv_content,
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

"""copy_safely and move_safely"""


@pytest.mark.parametrize("fn", [copy_safely, move_safely])
def test_copy_move_safely_pass_correct_content_at_destination(
    tmp_path: Path, fake_args: SimpleNamespace, fn: Callable
) -> None:
    """On success, the destination file is created with the source's content."""
    src = _touch(tmp_path / "src.txt", "hello world")
    dst = tmp_path / "dst.txt"

    assert fn(src, dst, fake_args) is True  # pyright: ignore[reportArgumentType]
    assert dst.read_text() == "hello world"
    assert fake_args.errors == []


@pytest.mark.parametrize("fn", [copy_safely, move_safely])
def test_copy_move_safely_pass_overwrites_existing_destination(
    tmp_path: Path, fake_args: SimpleNamespace, fn: Callable
) -> None:
    """An existing destination file is overwritten with the source's content."""
    src = _touch(tmp_path / "src.txt", "new content")
    dst = _touch(tmp_path / "dst.txt", "stale content")

    assert fn(src, dst, fake_args) is True  # pyright: ignore[reportArgumentType]
    assert dst.read_text() == "new content"


def test_copy_safely_pass_leaves_source_file_intact(tmp_path: Path, fake_args: SimpleNamespace) -> None:
    """Unlike move_safely, copy_safely must not remove the source file."""
    src = _touch(tmp_path / "src.txt", "hello")
    dst = tmp_path / "dst.txt"

    assert copy_safely(src, dst, fake_args) is True  # pyright: ignore[reportArgumentType]
    assert src.exists()
    assert src.read_text() == "hello"


def test_move_safely_pass_removes_source_file(tmp_path: Path, fake_args: SimpleNamespace) -> None:
    """Unlike copy_safely, move_safely must remove the source file on success."""
    src = _touch(tmp_path / "src.txt", "hello")
    dst = tmp_path / "dst.txt"

    assert move_safely(src, dst, fake_args) is True  # pyright: ignore[reportArgumentType]
    assert not src.exists()


@pytest.mark.parametrize(
    ("make_src", "make_dst"),
    [
        pytest.param(
            lambda tmp_path: tmp_path / "missing-src.txt",
            lambda tmp_path: tmp_path / "dst.txt",
            id="missing-source-file",
        ),
        pytest.param(
            lambda tmp_path: _touch(tmp_path / "src.txt"),
            lambda tmp_path: tmp_path / "no-such-dir" / "dst.txt",
            id="missing-destination-directory",
        ),
    ],
)
@pytest.mark.parametrize("fn", [copy_safely, move_safely])
def test_copy_safely_fail_records_error_and_returns_false(
    tmp_path: Path,
    fake_args: SimpleNamespace,
    make_src: Callable[[Path], Path],
    make_dst: Callable[[Path], Path],
    fn: Callable,
) -> None:
    """On an OSError, copy_safely and move_safely return False and records a matching ErrorRecord."""
    src = make_src(tmp_path)
    dst = make_dst(tmp_path)
    src_existed = src.exists()
    assert src_existed == (src.name == "src.txt")

    assert fn(src, dst, fake_args) is False  # pyright: ignore[reportArgumentType]
    assert not dst.exists()
    # if the source existed before the failed operation, it must not have been consumed
    assert src.exists() == src_existed

    assert len(fake_args.errors) == 1
    error = fake_args.errors[0]
    assert isinstance(error, ErrorRecord)
    assert error.file == fake_args.file_name
    assert error.returncode is None
    assert "No such file or directory" in error.message


"""non_header_lines_present"""


NON_HEADER_LINES_CASES: list[ParameterSet] = [
    pytest.param("", False, id="empty-file"),
    pytest.param("header\n", False, id="header-only"),
    pytest.param("header", False, id="header-only-no-trailing-newline"),
    pytest.param("header\n\n", False, id="header-plus-trailing-blank-line"),
    pytest.param("\n\n\n", False, id="only-blank-lines"),
    pytest.param(" \n\t\n", False, id="only-whitespace-lines"),
    pytest.param("header\ndata\n", True, id="header-plus-one-data-line"),
    pytest.param("header\ndata1\ndata2\n", True, id="header-plus-multiple-data-lines"),
    pytest.param("header\n\ndata\n", True, id="blank-line-between-header-and-data-is-ignored"),
    pytest.param("header\n \ndata\n", True, id="whitespace-only-line-is-ignored"),
]


@pytest.mark.parametrize(("content", "expected"), NON_HEADER_LINES_CASES)
def test_non_header_lines_present_pass_matches_expected_result(tmp_path: Path, content: str, expected: bool) -> None:
    """Blank/whitespace-only lines are ignored when deciding if non-header lines are present."""
    path = _touch(tmp_path / "file.txt", content)
    assert non_header_lines_present(path) is expected


def test_non_header_lines_present_fail_missing_file(tmp_path: Path) -> None:
    """A missing file raises FileNotFoundError rather than being silently treated as empty."""
    with pytest.raises(FileNotFoundError, match="No such file or directory:"):
        non_header_lines_present(tmp_path / "does-not-exist.txt")


"""prepend_header"""


def test_prepend_header_pass_concatenates_header_and_data(tmp_path: Path) -> None:
    """The destination file contains the header file's contents followed by the data file's contents."""
    header_path = _touch(tmp_path / "header.txt", "col1\tcol2\n")
    data_path = _touch(tmp_path / "data.txt", "1\t2\n3\t4\n")
    dest = tmp_path / "combined.txt"

    prepend_header(header_path, data_path, dest)

    assert dest.read_text() == "col1\tcol2\n1\t2\n3\t4\n"


def test_prepend_header_pass_writes_bytes_exactly(tmp_path: Path) -> None:
    """Content is copied byte-for-byte, so non-UTF-8/binary content is handled correctly."""
    header_path = tmp_path / "header.bin"
    header_path.write_bytes(b"col1,col2\n")
    data_path = tmp_path / "data.bin"
    data_path.write_bytes(b"\x00\x01,\xff\n")
    dest = tmp_path / "combined.bin"

    prepend_header(header_path, data_path, dest)

    assert dest.read_bytes() == b"col1,col2\n\x00\x01,\xff\n"


def test_prepend_header_pass_overwrites_existing_destination(tmp_path: Path) -> None:
    """A pre-existing destination file is fully overwritten, not appended to."""
    header_path = _touch(tmp_path / "header.txt", "h\n")
    data_path = _touch(tmp_path / "data.txt", "d\n")
    dest = _touch(tmp_path / "combined.txt", "stale content that should be gone\n")

    prepend_header(header_path, data_path, dest)

    assert dest.read_text() == "h\nd\n"


def test_prepend_header_fail_missing_header_file(tmp_path: Path) -> None:
    """A missing header file raises FileNotFoundError."""
    header_path = tmp_path / "missing-header.txt"
    data_path = _touch(tmp_path / "data.txt", "d\n")
    dest = tmp_path / "combined.txt"

    with pytest.raises(FileNotFoundError, match="No such file or directory:"):
        prepend_header(header_path, data_path, dest)


def test_prepend_header_fail_missing_data_file(tmp_path: Path) -> None:
    """A missing data file raises FileNotFoundError (no partial output file is left dangling)."""
    header_path = _touch(tmp_path / "header.txt", "h\n")
    data_path = tmp_path / "missing-data.txt"
    dest = tmp_path / "combined.txt"

    with pytest.raises(FileNotFoundError, match="No such file or directory:"):
        prepend_header(header_path, data_path, dest)


"""generate_qsv_validate_file_names"""


@pytest.mark.parametrize("first_pass_value", [None, True, False])
@pytest.mark.parametrize("input_file_name", ["some_file_name.tsv", "a_random_sequence_of_chars"])
@pytest.mark.parametrize("ext", [".tsv", ".csv", ".tar.gz", ".exe"])
def test_generate_qsv_validate_file_names(first_pass_value: None | bool, input_file_name: str, ext: str) -> None:
    """Ensure that the appropriate file names are generated for a variety of values.

    Parametrizes the file extension, the input file name, and whether or not this is the first pass naming.
    """
    tmp_path = Path("tmp")
    cv_args = MagicMock(
        ext=ext,
        tmp_dir_path=tmp_path,
        xsv_file_base_name="some_name",
        invalid_file_suffix=f"invalid{ext}",
        valid_file_suffix=f"valid{ext}",
    )

    args = (cv_args, input_file_name) if first_pass_value is None else (cv_args, input_file_name, first_pass_value)

    # file_names_tuple consists of
    # valid_output_file, _errors_file, _valid_lines_file, _invalid_lines_file
    file_names = generate_qsv_validate_file_names(*args)

    valid_output_file_name = f"some_name-first-pass{VALID_SUFFIX}{ext}"
    schema_file_name = "some_name-first-pass.jsonschema.json"
    if first_pass_value is not None and first_pass_value is False:
        valid_output_file_name = f"some_name-normalised-validated{VALID_SUFFIX}{ext}"
        schema_file_name = "some_name-normalised-validated.jsonschema.json"

    assert file_names.valid_output == tmp_path / valid_output_file_name
    assert file_names.errors == tmp_path / f"{input_file_name}{VALIDATION_ERRORS}"
    assert file_names.valid_lines == tmp_path / f"{input_file_name}.valid{ext}"
    assert file_names.invalid_lines == tmp_path / f"{input_file_name}.invalid{ext}"
    assert file_names.schema == tmp_path / schema_file_name

    names_list = [
        file_names.valid_output,
        file_names.errors,
        file_names.valid_lines,
        file_names.invalid_lines,
        file_names.schema,
    ]
    # all generated names are different
    assert len(set(names_list)) == len(names_list)


@pytest.mark.parametrize(
    "input_file_name",
    [
        pytest.param("data-header.tsv", id="simple-name"),
        pytest.param("data.with.dots-header.tsv", id="name-with-embedded-dots"),
        pytest.param("data-normalised.norm.tsv", id="name-resembling-a-previous-suffix"),
    ],
)
def test_generate_qsv_validate_file_names_pass_handles_various_input_file_names(input_file_name: str) -> None:
    """input_file_name is used verbatim as a prefix for errors_file/valid_lines_file/invalid_lines_file.

    Embedded dots or substrings that resemble other pipeline suffixes are ignored.
    """
    tmp_path = Path("tmp")
    ext = ".tar.gz"
    cv_args = MagicMock(
        ext=ext,
        tmp_dir_path=tmp_path,
        xsv_file_base_name="some_name",
        invalid_file_suffix=f"invalid{ext}",
        valid_file_suffix=f"valid{ext}",
    )

    file_names = generate_qsv_validate_file_names(cv_args, input_file_name)

    valid_output_file_name = f"some_name-first-pass{VALID_SUFFIX}{ext}"
    schema_file_name = "some_name-first-pass.jsonschema.json"

    assert file_names.valid_output == tmp_path / valid_output_file_name
    assert file_names.errors == tmp_path / f"{input_file_name}{VALIDATION_ERRORS}"
    assert file_names.valid_lines == tmp_path / f"{input_file_name}.valid{ext}"
    assert file_names.invalid_lines == tmp_path / f"{input_file_name}.invalid{ext}"
    assert file_names.schema == tmp_path / schema_file_name

    names_list = [
        file_names.valid_output,
        file_names.errors,
        file_names.valid_lines,
        file_names.invalid_lines,
        file_names.schema,
    ]
    # all generated names are different
    assert len(set(names_list)) == len(names_list)


"""ErrorRecord"""


def test_error_record_from_qsv_result_captures_file_message_and_returncode() -> None:
    """from_qsv_result pulls the file name, stripped stderr, and returncode straight from the completed process."""
    result = subprocess.CompletedProcess(args=["qsv", "validate"], returncode=1, stdout=None, stderr="boom\n")

    record = ErrorRecord.from_qsv_result("data.tsv", result)

    assert record.file == "data.tsv"
    assert record.message == "boom"
    assert record.returncode == 1


def test_error_record_from_qsv_result_strips_surrounding_whitespace_from_stderr() -> None:
    """Leading/trailing whitespace (including embedded newlines) around stderr is stripped."""
    result = subprocess.CompletedProcess(args=["qsv"], returncode=2, stdout=None, stderr="\n\n  multi-line error  \n\n")

    record = ErrorRecord.from_qsv_result("data.tsv", result)

    assert record.message == "multi-line error"


def test_error_record_from_qsv_result_handles_none_stderr_as_empty_message() -> None:
    """A None stderr (e.g. if it wasn't captured) results in an empty string message rather than raising."""
    result = subprocess.CompletedProcess(args=["qsv"], returncode=0, stdout=None, stderr=None)

    record = ErrorRecord.from_qsv_result("data.tsv", result)

    assert record.message == ""
    assert record.returncode == 0


def test_error_record_from_exception_captures_file_and_stringified_message() -> None:
    """from_exception stores the exception's str() as the message, leaving returncode unset (None)."""
    exc = OSError("No such file or directory: 'missing.tsv'")

    record = ErrorRecord.from_exception("data.tsv", exc)

    assert record.file == "data.tsv"
    assert record.message == "No such file or directory: 'missing.tsv'"
    assert record.returncode is None


def test_error_record_from_exception_works_with_arbitrary_exception_types() -> None:
    """from_exception isn't limited to OSError; any Exception subclass is supported."""
    exc = ValueError("some validation problem")

    record = ErrorRecord.from_exception("data.tsv", exc)

    assert record.message == "some validation problem"
    assert record.returncode is None


"""CleanerValidatorArgs -- defaults"""


def test_cleaner_validator_args_defaults_are_applied_when_omitted(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """Optional fields fall back to their documented defaults when not explicitly supplied."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)

    assert args.delimiter == "\t"
    assert args.comment_char == "#"
    assert args.quote is None
    assert args.escape is None
    assert args.null_regex is None
    assert args.null_regex_cols is None
    assert args.missing_header is False
    assert args.errors == []


def test_cleaner_validator_args_errors_default_factory_is_not_shared_between_instances(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """Each instance gets its own independent `errors` list; mutating one must not affect another.

    Pins down that the (unusually-shaped) `default_factory=lambda _: []` produces a fresh list per
    instance rather than a single shared mutable default.
    """
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    first = make_mock_args(source)
    second = make_mock_args(source)

    first.errors.append(ErrorRecord(file="a.tsv", message="oops"))

    assert first.errors != second.errors
    assert second.errors == []


"""CleanerValidatorArgs -- computed fields"""


def test_cleaner_validator_args_qsv_env_includes_comment_char_and_inherits_os_environ(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """qsv_env merges the process environment with QSV_COMMENT_CHAR set to comment_char."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, comment_char=";")

    env = args.qsv_env

    assert env["QSV_COMMENT_CHAR"] == ";"
    # every pre-existing environment variable is inherited too
    for key, value in os.environ.items():
        assert env[key] == value


def test_cleaner_validator_args_qsv_env_is_cached_after_first_access(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """qsv_env is a cached_property: mutating comment_char after first access does not change the cached value."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source, comment_char="#")

    first_env = args.qsv_env
    assert first_env["QSV_COMMENT_CHAR"] == "#"

    args.comment_char = ";"

    assert args.qsv_env is first_env
    assert args.qsv_env["QSV_COMMENT_CHAR"] == "#"


def test_cleaner_validator_args_xsv_file_base_name_and_file_name(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """xsv_file_base_name is the file stem; file_name is the full file name including extension."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "some-data.tsv")
    args = make_mock_args(source)

    assert args.xsv_file_base_name == "some-data"
    assert args.file_name == "some-data.tsv"


@pytest.mark.parametrize("delimiter", DELIMITERS)
def test_cleaner_validator_args_ext_uses_sep_to_ext_mapping_regardless_of_actual_file_suffix(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile, delimiter: str
) -> None:
    """For any delimiter known to SEP_TO_EXT, `ext` reflects the mapped extension.

    This is true even when the source file's suffix (`.weird`) is not the correct extension for the format.
    """
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS, delimiter=delimiter), "data.weird")
    args = make_mock_args(source, delimiter=delimiter)

    assert args.ext == SEP_TO_EXT[delimiter]


def test_cleaner_validator_args_ext_falls_back_to_file_suffix_for_unmapped_delimiter(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """A delimiter not present in SEP_TO_EXT (e.g. a pipe) falls back to the source file's actual suffix."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS, delimiter="|"), "data.psv")
    args = make_mock_args(source, delimiter="|")

    assert "|" not in SEP_TO_EXT
    assert args.ext == ".psv"


def test_cleaner_validator_args_valid_and_invalid_file_suffixes_are_derived_from_ext(
    make_mock_args: Callable[..., CleanerValidatorArgs], write_source_file: WriteFile
) -> None:
    """valid_file_suffix/invalid_file_suffix are literally 'valid'/'invalid' with `ext` appended."""
    source = write_source_file(build_xsv_content(VALID_ROWS, header=COLUMNS), "data.tsv")
    args = make_mock_args(source)

    assert args.valid_file_suffix == f"valid{args.ext}"
    assert args.invalid_file_suffix == f"invalid{args.ext}"
