"""Helpers for the qsv validator."""

import os
import subprocess
from dataclasses import dataclass
from functools import cached_property
from logging import Logger, getLogger
from pathlib import Path
from shutil import copy, move
from typing import Annotated, Final, Self

from pydantic import (
    BaseModel,
    DirectoryPath,
    Field,
    FilePath,
    StringConstraints,
    computed_field,
)

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import ValidatedSchema

NonEmptyStr = Annotated[str, StringConstraints(min_length=1)]
CharStr = Annotated[str, StringConstraints(min_length=1, max_length=1)]


VALIDATION_ERRORS: Final[str] = ".validation-errors.tsv"

SEP_TO_EXT: Final[dict[str, str]] = {"\t": ".tsv", ",": ".csv", ";": ".ssv"}

# Suffixes used when constructing intermediate/output file names for each pipeline stage.
HEADER_SUFFIX: Final[str] = "-header"
CLEANED_SUFFIX: Final[str] = "-cleaned"
NORM_SUFFIX: Final[str] = "-norm"
VALID_SUFFIX: Final[str] = "-valid"


logger: Logger = getLogger(__name__)


class ErrorRecord(BaseModel):
    """A single error encountered while cleaning/validating a file.

    Used so that ``errors`` always contains a uniform type, regardless of whether the
    error came from a failed qsv subprocess, an ``OSError`` raised while copying/moving files, or
    any other type of fun error that might occur.

    :param file: name of the xsv file being processed when the error occurred
    :type file: str
    :param message: human-readable description of the error
    :type message: str
    :param returncode: return code of the qsv subprocess, if the error originated from one
    :type returncode: int | None, optional
    """

    file: str
    message: str
    returncode: int | None = None

    @classmethod
    def from_qsv_result(cls, file_name: str, result: subprocess.CompletedProcess[str]) -> Self:
        """Build an `ErrorRecord` from a failed qsv subprocess result."""
        return cls(file=file_name, message=(result.stderr or "").strip(), returncode=result.returncode)

    @classmethod
    def from_exception(cls, file_name: str, exc: Exception) -> Self:
        """Build an `ErrorRecord` from a raised exception (e.g. `OSError`)."""
        return cls(file=file_name, message=str(exc))


class CleanerValidatorArgs(BaseModel):
    """Class containing required parameters for running the XSV file clean/validate functions.

    :param qsv_cmd: command for executing qsv
    :type qsv_cmd: str
    :param xsv_file_path: path to the data file to be validated
    :type xsv_file_path: FilePath
    :param header_file_path: path to the generated header file
    :type header_file_path: FilePath

    :param first_pass_schema: schema for initial validation of the data
    :type first_pass_schema: ValidatedSchema
    :param post_norm_schema: schema for validating data after cleaning and normalisation
    :type post_norm_schema: ValidatedSchema

    :param tmp_dir_path: temporary working directory
    :type tmp_dir_path: DirectoryPath
    :param qsv_output_dir_path: directory for qsv output from failed executions
    :type qsv_output_dir_path: DirectoryPath
    :param validated_file_dir_path: directory to put the cleaned, validated files in
    :type validated_file_dir_path: DirectoryPath
    :param errors: list of errors accumulated during processing
    :type errors: list[ErrorRecord]
    :param delimiter: delimiter used in the xsv file, defaults to tab
    :type delimiter: CharStr, optional
    :param comment_char: comment character used in the xsv file, defaults to `#`
    :type comment_char: CharStr, optional
    :param quote: quote character used in the xsv file; if unset, qsv will assume it is `"`
    :type quote: CharStr, optional
    :param escape: the escape character used in the xsv file. If not specified, qsv escapes quotes by doubling them. Default None
    :type escape: NonEmptyStr, optional
    :param null_regex: regex to replace with an empty string
    :type null_regex: NonEmptyStr, optional
    :param null_regex_cols: list of columns to apply the null regex to, optional; defaults to None (all cols)
    :type null_regex_cols: list[NonEmptyStr] | None
    :param missing_header: whether or not the xsv file has a missing header, defaults to False
    :type missing_header: bool, optional
    """

    qsv_cmd: NonEmptyStr

    xsv_file_path: FilePath
    header_file_path: FilePath

    first_pass_schema: ValidatedSchema
    post_norm_schema: ValidatedSchema

    tmp_dir_path: DirectoryPath
    qsv_output_dir_path: DirectoryPath
    validated_file_dir_path: DirectoryPath

    errors: list[ErrorRecord] = Field(default_factory=lambda _: [])

    delimiter: CharStr = "\t"
    comment_char: CharStr = "#"
    quote: CharStr | None = None
    escape: NonEmptyStr | None = None
    null_regex: NonEmptyStr | None = None
    null_regex_cols: list[NonEmptyStr] | None = None
    missing_header: bool = False

    @computed_field
    @cached_property
    def qsv_env(self: Self) -> dict[str, str]:
        """Environment for qsv commands to run in."""
        return {**os.environ, "QSV_COMMENT_CHAR": self.comment_char}

    @computed_field
    @property
    def xsv_file_base_name(self: Self) -> str:
        """Base name of the input xsv file, without the final suffix."""
        return self.xsv_file_path.stem

    @computed_field
    @property
    def file_name(self) -> str:
        """Name of the input xsv file."""
        return self.xsv_file_path.name

    @computed_field
    @property
    def ext(self) -> str:
        """Extension for the xsv files; includes the `.`."""
        return SEP_TO_EXT.get(self.delimiter) or self.xsv_file_path.suffix

    @computed_field
    @property
    def valid_file_suffix(self) -> str:
        """Suffix added to validator output file containing valid lines."""
        return f"valid{self.ext}"

    @computed_field
    @property
    def invalid_file_suffix(self) -> str:
        """Suffix added to validator output file containing valid lines."""
        return f"invalid{self.ext}"


def copy_safely(src: Path, dst: Path, args: "CleanerValidatorArgs") -> bool:
    """Copy `src` to `dst`, recording an `ErrorRecord` and returning False on failure.

    :param src: source path
    :type src: Path
    :param dst: destination path
    :type dst: Path
    :param args: args for the cleaner validator
    :type args: CleanerValidatorArgs
    :return: True if the copy succeeded, False otherwise
    :rtype: bool
    """
    try:
        copy(src, dst)
    except OSError as exc:
        args.errors.append(ErrorRecord.from_exception(args.file_name, exc))
        return False
    return True


def move_safely(src: Path, dst: Path, args: "CleanerValidatorArgs") -> bool:
    """Move `src` to `dst`, recording an `ErrorRecord` and returning False on failure.

    :param src: source path
    :type src: Path
    :param dst: destination path
    :type dst: Path
    :param args: args for the cleaner validator
    :type args: CleanerValidatorArgs
    :return: True if the move succeeded, False otherwise
    :rtype: bool
    """
    try:
        move(src, dst)
    except OSError as exc:
        args.errors.append(ErrorRecord.from_exception(args.file_name, exc))
        return False
    return True


def non_header_lines_present(path: Path) -> bool:
    """Return True if `path` has any non-blank lines beyond the header line.

    Stops reading as soon as a second non-blank line is found.
    """
    non_blank_lines_seen = 0
    with path.open() as f:
        for line in f:
            if line.strip():
                non_blank_lines_seen += 1
                if non_blank_lines_seen > 1:
                    return True
    return False


def prepend_header(header_file_path: Path, xsv_file_path: Path, dest: Path) -> None:
    """Write a new file at `dest` consisting of the header file's contents followed by the data file's contents."""
    with dest.open("wb") as out:
        out.write(header_file_path.read_bytes())
        out.write(xsv_file_path.read_bytes())


@dataclass
class FileNames:
    """Dataclass for storing derived file names."""

    valid_output: Path
    errors: Path
    valid_lines: Path
    invalid_lines: Path
    schema: Path


def generate_qsv_validate_file_names(
    args: CleanerValidatorArgs,
    input_file_name: str,
    first_pass: bool = True,  # noqa: FBT001, FBT002
) -> FileNames:
    """Compute the file names run_qsv_validate generates from a file with validation errors.

    :param args: arguments
    :type args: CleanerValidatorArgs
    :param input_file_name: input file name
    :type input_file_name: str
    :param first_pass: whether or not this is the first validation pass, defaults to True
    :type first_pass: bool, optional
    :return: dataclass containing derived paths (valid_output, errors, valid_lines, invalid_lines, schema)
    :rtype: FileNames
    """
    # file name stem for validated lines and schema: {file_base_name}-header-valid.ext
    output_file_stem = (
        f"{args.xsv_file_base_name}-first-pass" if first_pass else f"{args.xsv_file_base_name}-normalised-validated"
    )

    return FileNames(
        valid_output=args.tmp_dir_path / f"{output_file_stem}{VALID_SUFFIX}{args.ext}",
        errors=args.tmp_dir_path / f"{input_file_name}{VALIDATION_ERRORS}",
        valid_lines=args.tmp_dir_path / f"{input_file_name}.{args.valid_file_suffix}",
        invalid_lines=args.tmp_dir_path / f"{input_file_name}.{args.invalid_file_suffix}",
        schema=args.tmp_dir_path / f"{output_file_stem}.jsonschema.json",
    )
