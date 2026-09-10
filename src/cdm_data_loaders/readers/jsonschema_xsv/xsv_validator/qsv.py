"""Core qsv-interacting code for the xsv validator."""

import json
import os
import shutil
import subprocess

from pydantic import validate_call

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.helpers import (
    CLEANED_SUFFIX,
    HEADER_SUFFIX,
    NORM_SUFFIX,
    CleanerValidatorArgs,
    ErrorRecord,
    copy_safely,
    generate_qsv_validate_file_names,
    move_safely,
    non_header_lines_present,
    prepend_header,
)
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import ValidatedSchema


def qsv_check() -> str | None:
    """Check whether the qsv binary is available and whether it is functioning as expected."""
    # check for qsv
    qsv_cmd = shutil.which("qsv") or os.environ.get("QSV_BIN")
    if qsv_cmd:
        err_msg = None
        try:
            # check qsv is functional
            result = subprocess.run([qsv_cmd, "--version"], text=True, stderr=subprocess.PIPE, check=True)
            if result.returncode == 0:
                return qsv_cmd
            err_msg = f"`qsv --version` exited with code {result.returncode}"
            if result.stderr:
                err_msg += f"; STDERR: {result.stderr.strip()}"
        except Exception as e:
            err_msg = f"Cannot perform validation with qsv: {e!s}"
            raise RuntimeError(err_msg) from e
    else:
        err_msg = "Could not locate the qsv binary"
    raise RuntimeError(err_msg)


@validate_call
def clean_validate_file(args: CleanerValidatorArgs) -> None | str:
    """Clean and validate an xSV file.

    Files undergo the following steps:

    - headers are added if missing_header=True
    - preliminary validation removes ragged rows
    - files are normalised using xsv input
    - nulls are replaced if null_regex != None
    - second validation occurs with stricter criteria to ensure fields match

    Files are validated in two stages: the first pass schema is used for very rough validation and basically
    ensures that there are the correct number of columns in each row.

    The second validation occurs after the files are tidied up, nulls removed, and so on.

    Note: there is currently no way to set custom escape characters or quotes in `qsv validate`, so lines with
    a delimiter character as part of the field may fail the first pass validation for having too many cols.

    :param args: object with all the useful params required for running the clean/validate workflow
    :type args: CleanerValidatorArgs
    """
    file_with_headers = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"

    if args.missing_header:
        # add header and save in tmp dir
        try:
            prepend_header(args.header_file_path, args.xsv_file_path, args.tmp_dir_path / file_with_headers)
        except OSError as exc:
            args.errors.append(ErrorRecord.from_exception(args.file_name, exc))
            return None
    # otherwise, just copy the file into the temp dir with the appropriate extension
    elif not copy_safely(args.xsv_file_path, args.tmp_dir_path / file_with_headers, args):
        return None

    file_for_input_cmd = run_qsv_validate(args, file_with_headers, schema=args.first_pass_schema, first_pass=True)

    if not file_for_input_cmd:
        return None

    file_for_rplc_cmd = run_qsv_input(args, file_for_input_cmd)
    if not file_for_rplc_cmd:
        return None

    # run qsv replace
    # if there is nothing specified as null_regex, do not run the command
    if args.null_regex:
        file_for_validate_cmd = run_qsv_null_replacement(args, file_for_rplc_cmd)
        if not file_for_validate_cmd:
            return None
    else:
        file_for_validate_cmd = file_for_rplc_cmd

    validated_file = run_qsv_validate(args, file_for_validate_cmd, args.post_norm_schema, first_pass=False)

    # save to output dir and return
    if validated_file and move_safely(
        args.tmp_dir_path / validated_file, args.validated_file_dir_path / validated_file, args
    ):
        return validated_file
    return None


def _run_qsv_step(cmd: list[str], args: "CleanerValidatorArgs") -> subprocess.CompletedProcess[str] | None:
    """Run a qsv subcommand, adding in the qsv environment.

    If an error occurs during execution of the command, it is stored in the `args.errors` array.

    :param cmd: qsv command and arguments to run
    :type cmd: list[str]
    :param args: args for the cleaner validator
    :type args: CleanerValidatorArgs
    :return: the completed process or None if some sort of error occurred during execution
    :rtype: subprocess.CompletedProcess[str] | None
    """
    try:
        return subprocess.run(cmd, text=True, stderr=subprocess.PIPE, env=args.qsv_env)
    except Exception as exc:  # noqa: BLE001
        args.errors.append(ErrorRecord.from_exception(args.file_name, exc))
    return None


def run_qsv_input(args: CleanerValidatorArgs, input_file_name: str) -> str | None:
    """Run qsv input to normalise the file input_file_name.

    See https://github.com/dathere/qsv/blob/master/docs/help/input.md for more information.

    Does not currently take account of quote or escape characters.

    :param args: args for the cleaner validator
    :type args: CleanerValidatorArgs
    :param input_file_name: name of the input file
    :type input_file_name: str
    :return: name of the output file if all went ok; otherwise None
    :rtype: str | None
    """
    if not (args.tmp_dir_path / input_file_name).is_file():
        err_msg = f"Input file not found at {(args.tmp_dir_path / input_file_name)!s}"
        args.errors.append(ErrorRecord.from_exception(args.file_name, FileNotFoundError(err_msg)))
        return None

    output_file_name = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    output_path = args.tmp_dir_path / output_file_name

    xsv_file_config = ["--delimiter", args.delimiter]
    for conf in ["quote", "escape"]:
        if value := getattr(args, conf):
            xsv_file_config.extend([f"--{conf}", value])

    # comment character is in the environment
    cmd = [
        args.qsv_cmd,
        "input",
        "--trim-headers",
        "--trim-fields",
        "--encoding-errors",
        "strict",
        *xsv_file_config,
        "--output",
        str(output_path),
        # input file
        str(args.tmp_dir_path / input_file_name),
    ]
    result = _run_qsv_step(cmd, args)
    if result is None:
        return None

    if result.returncode == 0:
        return output_file_name

    # any other return code: record the error
    # note: the output file will still exist, but it is ignored
    args.errors.append(ErrorRecord.from_qsv_result(args.file_name, result))
    return None


def run_qsv_null_replacement(args: CleanerValidatorArgs, input_file_name: str) -> str | None:
    """Run qsv replace to replace any exact matches of `args.null_regex` in input_file_name with "".

    :param args: args for the cleaner validator
    :type args: CleanerValidatorArgs
    :param input_file_name: name of the input file
    :type input_file_name: str
    :return: name of the output file if all went ok; otherwise None
    :rtype: str | None
    """
    if not (args.tmp_dir_path / input_file_name).is_file():
        err_msg = f"Input file not found at {(args.tmp_dir_path / input_file_name)!s}"
        args.errors.append(ErrorRecord.from_exception(args.file_name, FileNotFoundError(err_msg)))
        return None

    output_file_name = f"{args.xsv_file_base_name}{NORM_SUFFIX}{args.ext}"
    regex_cols = ["-s", f"{','.join(args.null_regex_cols)}"] if args.null_regex_cols else []

    cmd = [
        args.qsv_cmd,
        "replace",
        # qsv returns 0 even if there are no matches
        "--not-one",
        "--delimiter",
        args.delimiter,
        "--output",
        str(args.tmp_dir_path / output_file_name),
        *regex_cols,
        # pattern
        args.null_regex,
        # replacement
        "",
        # input file
        str(args.tmp_dir_path / input_file_name),
    ]
    result = _run_qsv_step(cmd, args)
    if result is None:
        return None

    if result.returncode == 0:
        return output_file_name
    # any other return code: record the error
    # note: the output file may still exist, but it is ignored
    args.errors.append(ErrorRecord.from_qsv_result(args.file_name, result))
    # copy the input file to the output directory for later perusal
    copy_safely(args.tmp_dir_path / input_file_name, args.qsv_output_dir_path / input_file_name, args)
    return None


def run_qsv_validate(
    args: CleanerValidatorArgs,
    input_file_name: str,
    schema: ValidatedSchema,
    first_pass: bool = True,  # noqa: FBT001, FBT002
) -> str | None:
    """Run the validator over the file input_file_name.

    :param args: args for the cleaner validator
    :type args: CleanerValidatorArgs
    :param input_file_name: name of the input file
    :type input_file_name: str
    :param schema: ValidatedSchema object to use for validation
    :type schema: ValidatedSchema
    :param first_pass: if True, this runs the first pass validator, which does not check format
    :type first_pass: bool, defaults to True
    :return: name of the output file if all went ok; otherwise None
    :rtype: str | None
    """
    if not (args.tmp_dir_path / input_file_name).is_file():
        err_msg = f"Input file not found at {(args.tmp_dir_path / input_file_name)!s}"
        args.errors.append(ErrorRecord.from_exception(args.file_name, FileNotFoundError(err_msg)))
        return None

    # generate the full paths for the derived files
    files = generate_qsv_validate_file_names(args, input_file_name, first_pass=first_pass)

    # save the schema to a temp file
    files.schema.write_text(json.dumps(schema.jsonschema))

    cmd = [
        args.qsv_cmd,
        "validate",
        "--delimiter",
        args.delimiter,
        # only validate field format when running the post-normalisation validation
        "--no-format-validation" if first_pass else None,
        # ragged rows are put in the invalid lines file
        "--split-ragged",
        # trim whitespace from fields
        "--trim",
        "--valid",
        args.valid_file_suffix,
        "--invalid",
        args.invalid_file_suffix,
        "--valid-output",
        str(files.valid_output),
        # file to validate
        str(args.tmp_dir_path / input_file_name),
        # the schema
        str(files.schema),
    ]
    result = _run_qsv_step([c for c in cmd if c], args)
    if result is None:
        return None

    # check whether the output was all valid -- if so, --valid-output will have been produced
    if files.valid_output.is_file():
        return files.valid_output.name

    # otherwise, the return code will have useful info, so save STDERR as an error message
    if result.returncode != 0:
        args.errors.append(ErrorRecord.from_qsv_result(args.file_name, result))

    # if errors are found, produces three files:
    # {input_file_name}{VALIDATION_ERRORS} -- list of all the errors found
    # {input_file_name}{args.valid_file_suffix} -- valid lines
    # {input_file_name}{args.invalid_file_suffix}-- lines that fail validation
    # comments are automatically removed and do not appear in .valid or .invalid files
    # headers appear in the valid_lines_file

    valid_lines_found = files.valid_lines.exists() and non_header_lines_present(files.valid_lines)

    for f in [files.errors, files.invalid_lines]:
        if f.is_file():
            # copy to the output directory for later consumption
            copy_safely(f, args.qsv_output_dir_path / f.name, args)

    # if the valid lines file exists, has more than just a header line, and can be successfully renamed,
    # return the new file name
    if valid_lines_found:
        # copy to the output folder
        copy_safely(files.valid_lines, args.qsv_output_dir_path / files.valid_lines.name, args)
        # rename the existing file for use as the output file
        if move_safely(files.valid_lines, files.valid_output, args):
            return files.valid_output.name

    return None
