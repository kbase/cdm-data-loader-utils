"""Fixtures and helpers for the xsv validator module parts."""

import json
import os
import shutil
import subprocess
from collections.abc import Callable, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Final
from unittest.mock import MagicMock

import pytest

from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator import qsv
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.helpers import (
    CLEANED_SUFFIX,
    HEADER_SUFFIX,
    SEP_TO_EXT,
    CleanerValidatorArgs,
)
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.schema_utils import (
    ValidatedSchema,
    generate_header,
)

FAKE_QSV_CMD: Final[str] = "mock-qsv"

SCHEMAS_DIR: Final[Path] = Path("tests") / "data" / "qsv" / "schemas"

Row = Sequence[str]
WriteFile = Callable[[str | bytes, str], Path]

# Canonical column order shared by both schema fixtures.
COLUMNS: Final[list[str]] = ["number", "date", "float", "boolean", "string"]

# valid delimiters from the extension <=> delimiter mapping
DELIMITERS: Final[list[str]] = list(SEP_TO_EXT)

VALID_ROWS: Final[list[list[str]]] = [
    ["2", "2023-01-15", "3.14", "true", "key:value1"],
    ["3", "2023-02-20", "1.11", "false", "key:value2"],
]


@pytest.fixture(scope="session")
def qsv_cmd() -> str:
    """Provide the path/name of the qsv binary; override via the QSV_BIN env var if not on PATH.

    If the qsv binary is not available, all tests that rely on this fixture are skipped.
    """
    cmd = shutil.which("qsv") or os.environ.get("QSV_BIN")
    if not cmd or not Path(cmd).is_file():
        pytest.skip(f"qsv binary not found at {cmd!r}; set QSV_BIN to point at a real binary")
    return cmd


@pytest.fixture(scope="session")
def first_pass_schema_path() -> Path:
    """Path to the first-pass JSON Schema fixture."""
    return SCHEMAS_DIR / "first_pass_schema.json"


@pytest.fixture(scope="session")
def first_pass_validated_schema(first_pass_schema_path: Path) -> ValidatedSchema:
    """ValidatedSchema object for the first pass schema."""
    return ValidatedSchema(jsonschema=json.loads(first_pass_schema_path.read_bytes()))


@pytest.fixture(scope="session")
def post_norm_schema_path() -> Path:
    """Path to the post-normalisation JSON Schema fixture."""
    return SCHEMAS_DIR / "post_norm_schema.json"


@pytest.fixture(scope="session")
def post_norm_validated_schema(post_norm_schema_path: Path) -> ValidatedSchema:
    """ValidatedSchema object for the post-norm schema."""
    return ValidatedSchema(jsonschema=json.loads(post_norm_schema_path.read_bytes()))


@pytest.fixture(scope="session")
def derived_first_pass_schema_path() -> Path:
    """First pass schema derived from post_norm_schema. Derivation performed manually."""
    return SCHEMAS_DIR / "derived_first_pass_schema.json"


@pytest.fixture(scope="session")
def derived_first_pass_schema(derived_first_pass_schema_path: Path) -> dict[str, Any]:
    """Content of the schema at derived_first_pass_schema_path."""
    return json.loads(derived_first_pass_schema_path.read_bytes())


@pytest.fixture
def make_schema_file(tmp_path: Path) -> Callable[[Any, str], Path]:
    """Factory fixture that dumps an arbitrary JSON-serialisable object to a temp schema file."""

    def _make(schema: list | dict, file_name: str = "schema.json") -> Path:
        path = tmp_path / file_name
        path.write_text(json.dumps(schema))
        return path

    return _make


@pytest.fixture
def tmp_dir_path(tmp_path: Path) -> Path:
    """Temporary working directory used as CleanerValidatorArgs.tmp_dir_path."""
    d = tmp_path / "tmp"
    d.mkdir()
    return d


@pytest.fixture
def output_dir_path(tmp_path: Path) -> Path:
    """Temporary output directory used as a base for qsv_output_dir_path and validated_file_dir_path."""
    d = tmp_path / "output"
    d.mkdir()
    return d


@pytest.fixture
def output_dir_path_derivatives(output_dir_path: Path) -> dict[str, Path]:
    """Qsv_output and validated directories."""
    (output_dir_path / "qsv_output").mkdir(exist_ok=True, parents=True)
    (output_dir_path / "validated").mkdir(exist_ok=True, parents=True)
    return {
        "qsv_output_dir_path": output_dir_path / "qsv_output",
        "validated_file_dir_path": output_dir_path / "validated",
    }


@pytest.fixture
def source_dir_path(tmp_path: Path) -> Path:
    """Directory representing the 'incoming' location of source xsv files, separate from tmp_dir_path."""
    d = tmp_path / "source"
    d.mkdir()
    return d


@pytest.fixture
def write_source_file(source_dir_path: Path) -> WriteFile:
    """Factory fixture that writes content to a new file in source_dir_path (an xsv "input" location)."""
    return lambda content, file_name: _write_file(source_dir_path, content, file_name)


@pytest.fixture
def write_working_file(tmp_dir_path: Path) -> WriteFile:
    """Factory fixture that writes content directly into tmp_dir_path (simulating an already-prepped file)."""
    return lambda content, file_name: _write_file(tmp_dir_path, content, file_name)


@pytest.fixture
def with_invalid_flag_injected(monkeypatch: pytest.MonkeyPatch) -> None:
    """Patch `_run_qsv_step` to insert an unrecognised CLI flag into `cmd` at position 2 before running it.

    The corrupted command is still executed against the real qsv binary via the original
    `_run_qsv_step`, so this exercises qsv's actual argument-parsing failure path (expected exit
    code > 1) end-to-end.
    """
    real_run_qsv_step = qsv._run_qsv_step  # noqa: SLF001

    def _patched(cmd: list[str], args: "CleanerValidatorArgs") -> subprocess.CompletedProcess[str] | None:
        corrupted_cmd = [*cmd[:2], "--some-invalid-flag", *cmd[2:]]
        return real_run_qsv_step(corrupted_cmd, args)

    monkeypatch.setattr(qsv, "_run_qsv_step", _patched)


@pytest.fixture
def write_raw_file(tmp_path: Path) -> Callable[[str, str], Path]:
    """Factory fixture that writes raw, unvalidated text content directly to a file (bypassing JSON encoding)."""

    def _write(content: str, file_name: str) -> Path:
        path = tmp_path / file_name
        path.write_text(content)
        return path

    return _write


@pytest.fixture
def mock_qsv_run(monkeypatch: pytest.MonkeyPatch) -> Callable[..., MagicMock]:
    """Factory fixture that replaces `qsv.subprocess.run` with a configurable fake.

    Returns a factory function accepting:

    - `returncode`: exit code the fake process should report
    - `stderr`: stderr text the fake process should report
    - `output_content`: if not None, content to write to the path following `output_flag` in the
      invoked command (simulating qsv actually producing an output file); if None, no file is
      written, simulating qsv failing before producing any output
    - `extra_files`: allows writing several files at once (e.g. for qsv validate, which produces three files)

    The installed mock itself is returned, so callers can assert on how it was invoked (e.g. the
    constructed command, or that it was never called at all).
    """

    def _install(
        returncode: int = 0,
        stderr: str = "",
        output_content: str | None = None,
        output_flag: str = "--output",
        extra_files: dict[Path, str] | None = None,
    ) -> MagicMock:
        def _fake_run(cmd: list[str], **_kwargs: Any) -> subprocess.CompletedProcess[str]:  # noqa: ANN401
            if output_content is not None and output_flag in cmd:
                output_path = Path(cmd[cmd.index(output_flag) + 1])
                output_path.write_text(output_content)
            for path, content in (extra_files or {}).items():
                path.write_text(content)
            return subprocess.CompletedProcess(cmd, returncode=returncode, stdout=None, stderr=stderr)

        mock = MagicMock(side_effect=_fake_run)
        monkeypatch.setattr(qsv.subprocess, "run", mock)
        return mock

    return _install


@pytest.fixture
def fake_args() -> SimpleNamespace:
    """Minimal stand-in for CleanerValidatorArgs, providing only what copy_safely/move_safely need."""
    return SimpleNamespace(errors=[], file_name="source.txt")


def _make_args_factory(
    qsv_cmd: str,
    tmp_dir_path: Path,
    output_dir_path_derivatives: dict[str, Path],
    first_pass_schema: ValidatedSchema,
    post_norm_schema: ValidatedSchema,
) -> Callable[..., CleanerValidatorArgs]:
    """Build a `_make_args` factory closed over a fixed `qsv_cmd`."""

    def _make_args(
        xsv_file_path: Path,
        delimiter: str = "\t",
        comment_char: str = "#",
        null_regex: str | None = None,
        missing_header: bool = False,
        first_pass_schema_override: ValidatedSchema | None = None,
        post_norm_schema_override: ValidatedSchema | None = None,
    ) -> CleanerValidatorArgs:
        # generate_header expects a (first pass) ValidatedSchema object
        header_file_path = generate_header(
            first_pass_schema_override or first_pass_schema,
            tmp_dir_path,
            header_file_name=f"header-{ord(delimiter)}.txt",
            delimiter=delimiter,
        )
        return CleanerValidatorArgs(
            errors=[],
            qsv_cmd=qsv_cmd,
            xsv_file_path=xsv_file_path,
            header_file_path=header_file_path,
            first_pass_schema=first_pass_schema_override or first_pass_schema,
            post_norm_schema=post_norm_schema_override or post_norm_schema,
            tmp_dir_path=tmp_dir_path,
            delimiter=delimiter,
            comment_char=comment_char,
            null_regex=null_regex,
            missing_header=missing_header,
            **output_dir_path_derivatives,
        )

    return _make_args


@pytest.fixture
def make_args(
    tmp_dir_path: Path,
    output_dir_path_derivatives: dict[str, Path],
    first_pass_validated_schema: ValidatedSchema,
    post_norm_validated_schema: ValidatedSchema,
    qsv_cmd: str,
) -> Callable[..., CleanerValidatorArgs]:
    """Factory fixture that builds a CleanerValidatorArgs, generating a matching header on the fly.

    `first_pass_schema`/`post_norm_schema` overrides only affect the schema path used by qsv
    validate commands; the header file is always generated from the known-good default schema.
    """
    return _make_args_factory(
        qsv_cmd, tmp_dir_path, output_dir_path_derivatives, first_pass_validated_schema, post_norm_validated_schema
    )


@pytest.fixture
def make_mock_args(
    tmp_dir_path: Path,
    output_dir_path_derivatives: dict[str, Path],
    first_pass_validated_schema: ValidatedSchema,
    post_norm_validated_schema: ValidatedSchema,
) -> Callable[..., CleanerValidatorArgs]:
    """Generate a set of CleanerValidatorArgs with an invalid qsv path.

    Use with `mock_qsv_run` for testing qsv commands with the qsv binary mocked out.
    """
    return _make_args_factory(
        FAKE_QSV_CMD, tmp_dir_path, output_dir_path_derivatives, first_pass_validated_schema, post_norm_validated_schema
    )


def _touch(path: Path, content: str = "hello") -> Path:
    """Write `content` to `path` and return it."""
    path.write_text(content)
    return path


def build_xsv_content(
    rows: Sequence[Row],
    *,
    delimiter: str = "\t",
    header: Sequence[str] | None = None,
) -> str:
    """Build delimited text content from a list of rows, optionally prefixed with a header row.

    :param rows: sequence of rows, each a sequence of field values
    :param delimiter: field delimiter used to join each row
    :param header: header columns to prepend; omitted entirely if None
    :return: the fully assembled file content, newline-terminated
    """
    lines: list[str] = []
    if header is not None:
        lines.append(delimiter.join(header))
    lines.extend(delimiter.join(row) for row in rows)
    return "\n".join(lines) + "\n"


def interleave_comments(content: str, comment_char: str) -> str:
    """Insert a comment line before every line in `content`, using the given comment character."""
    lines = content.splitlines()
    commented: list[str] = []
    for i, line in enumerate(lines):
        commented.append(f"{comment_char} comment line {i}")
        commented.append(line)
    # ensure that we have at least one comment line!
    comment_lines = [line for line in commented if line.startswith(f"{comment_char} comment line")]
    assert len(comment_lines) >= 1
    return "\n".join(commented) + "\n"


def _write_file(dir_path: Path, content: str | bytes, file_name: str) -> Path:
    """Write `content` (str or bytes) to `file_name` inside `dir_path`, returning the resulting path."""
    path = dir_path / file_name
    if isinstance(content, bytes):
        path.write_bytes(content)
    else:
        path.write_text(content)
    return path


def _write_header_input_file(args: CleanerValidatorArgs, write_working_file: WriteFile) -> str:
    """Write a `-header`-suffixed working file for `args` and return its name."""
    input_file_name = f"{args.xsv_file_base_name}{HEADER_SUFFIX}{args.ext}"
    write_working_file(build_xsv_content(VALID_ROWS, header=COLUMNS), input_file_name)
    return input_file_name


def _write_cleaned_input_file(
    args: CleanerValidatorArgs, write_working_file: WriteFile, rows: Sequence[Row], delimiter: str | None = None
) -> str:
    """Write a `-cleaned`-suffixed working file for `args` (mimicking run_qsv_input's output) and return its name."""
    input_file_name = f"{args.xsv_file_base_name}{CLEANED_SUFFIX}{args.ext}"
    content = build_xsv_content(rows, header=COLUMNS, delimiter=delimiter or args.delimiter)
    write_working_file(content, input_file_name)
    return input_file_name


def _write_validator_input_file(
    args: CleanerValidatorArgs, write_working_file: WriteFile, rows: Sequence[Row], suffix: str
) -> str:
    """Write a working file with the given suffix (e.g. HEADER_SUFFIX, NORM_SUFFIX) for use as run_qsv_validate input."""
    input_file_name = f"{args.xsv_file_base_name}{suffix}{args.ext}"
    write_working_file(build_xsv_content(rows, header=COLUMNS), input_file_name)
    return input_file_name


def snapshot_dir(dir_path: Path) -> set[str]:
    """Return the set of file names currently present in `dir_path`, for before/after diffing."""
    return {p.name for p in dir_path.iterdir()}


def parse_xsv(content: str, delimiter: str = "\t") -> list[list[str]]:
    """Parse simple xsv content (no embedded delimiters/quotes) into a list of rows of fields."""
    return [line.split(delimiter) for line in content.splitlines() if line]
