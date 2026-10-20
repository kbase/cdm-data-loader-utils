"""Tests for the local path validators and the LocalInputDir / LocalOutputDir types."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from pydantic import BaseModel, ValidationError
from pydantic_core import ErrorDetails

from cdm_data_loaders.core.path_fields import (
    DIR_NOT_CREATED,
    PATH_EMPTY,
    PATH_NOT_LOCAL,
    VALID_PATH_ONE_CHAR,
    LocalInputDir,
    LocalOutputDir,
    check_local_path,
    ensure_dir,
    is_local,
    normalise_dir,
    url_protocol,
)


class InputDirModel(BaseModel):
    """Model with one LocalInputDir field."""

    path: LocalInputDir


class OutputDirModel(BaseModel):
    """Model with one LocalOutputDir field."""

    path: LocalOutputDir


BOTH_MODELS = [
    pytest.param(InputDirModel, id="LocalInputDir"),
    pytest.param(OutputDirModel, id="LocalOutputDir"),
]


def get_single_error(exc: ValidationError) -> ErrorDetails:
    """Check that a ValidationError holds exactly one error on the `path` field, and return it."""
    errors = exc.errors()
    assert len(errors) == 1
    assert errors[0]["loc"] == ("path",)
    return errors[0]


URL_PROTOCOL_CASES = [
    pytest.param("/tmp/data", "file", id="absolute_posix_path"),
    pytest.param("relative/dir", "file", id="relative_path"),
    pytest.param("", "file", id="empty_string"),
    pytest.param("C:\\data", "file", id="windows_drive_backslash"),
    pytest.param("C://data", "file", id="windows_drive_single_char_protocol"),
    pytest.param("file:///tmp/data", "file", id="file_protocol"),
    pytest.param("FILE:///tmp/data", "file", id="file_protocol_upper_case"),
    pytest.param("file://", "file", id="file_protocol_only"),
    pytest.param("s3://bucket/key", "s3", id="s3"),
    pytest.param("S3://bucket/key", "s3", id="s3_upper_case"),
    pytest.param("s3a://bucket/key", "s3a", id="s3a_not_aliased"),
    pytest.param("gs://bucket/key", "gs", id="gcs"),
    pytest.param("https://example.com/data", "https", id="https"),
    pytest.param("data:text/plain,hello", "data", id="data_uri"),
]


@pytest.mark.parametrize(("url", "expected"), URL_PROTOCOL_CASES)
def test_url_protocol_pass(url: str, expected: str) -> None:
    """Test that url_protocol returns the lower-case protocol, defaulting to 'file'."""
    assert url_protocol(url) == expected


@pytest.mark.parametrize(("url", "protocol"), URL_PROTOCOL_CASES)
def test_is_local_pass(url: str, protocol: str) -> None:
    """Test that is_local is True only for paths with the file protocol."""
    assert is_local(url) is (protocol == "file")


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param("/tmp/data", "/tmp/data", id="absolute_posix_path"),
        pytest.param("relative/dir", "relative/dir", id="relative_path"),
        pytest.param(".", ".", id="current_dir"),
        pytest.param("C:\\data", "C:\\data", id="windows_drive_path"),
        pytest.param("  /tmp/data \n", "/tmp/data", id="surrounding_whitespace"),
        pytest.param("file:///tmp/data", "/tmp/data", id="file_protocol_absolute"),
        pytest.param("FILE:///tmp/data", "/tmp/data", id="file_protocol_upper_case"),
        pytest.param("file://relative/dir", "relative/dir", id="file_protocol_relative"),
        pytest.param("  file:///tmp/data  ", "/tmp/data", id="file_protocol_with_whitespace"),
    ],
)
def test_check_local_path_pass_strings(value: str, expected: str) -> None:
    """Test that local string paths are stripped of whitespace and the file protocol."""
    assert check_local_path(value) == expected


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(Path("/tmp/data"), id="path_object"),
        pytest.param(None, id="none"),
        pytest.param(123, id="int"),
        pytest.param(b"/tmp/data", id="bytes"),
    ],
)
def test_check_local_path_pass_non_string_unchanged(value: Any) -> None:  # noqa: ANN401
    """Test that non-string values are returned as the same object."""
    assert check_local_path(value) is value


@pytest.mark.parametrize(
    "value",
    [
        pytest.param("https://example.com/data", id="https"),
        pytest.param("s3://bucket/key", id="s3"),
        pytest.param("s3a://bucket/key", id="s3a"),
        pytest.param("data:text/plain,hello", id="data_uri"),
        pytest.param("  gs://bucket/key  ", id="gcs_with_whitespace"),
    ],
)
def test_check_local_path_fail_not_local(value: str) -> None:
    """Test that non-local URLs raise a ValueError that quotes the stripped value."""
    with pytest.raises(ValueError, match="Path must be on the local filesystem") as exc_info:
        check_local_path(value)
    assert str(exc_info.value) == f"{PATH_NOT_LOCAL}: {value.strip()!r}"


@pytest.mark.parametrize(
    "value",
    [
        pytest.param("", id="empty_string"),
        pytest.param("   ", id="spaces_only"),
        pytest.param("\t\n", id="tab_newline_only"),
        pytest.param("file://", id="file_protocol_only"),
        pytest.param("FILE://", id="file_protocol_only_upper_case"),
        pytest.param("  file://  ", id="file_protocol_only_with_whitespace"),
    ],
)
def test_check_local_path_fail_empty(value: str) -> None:
    """Test that empty paths, including a bare file protocol, raise a ValueError."""
    with pytest.raises(ValueError, match="Path must not be empty") as exc_info:
        check_local_path(value)
    assert str(exc_info.value) == PATH_EMPTY


@pytest.mark.parametrize(
    "parts",
    [
        pytest.param((), id="existing_dir"),
        pytest.param(("new",), id="new_dir"),
        pytest.param(("a", "b", "c"), id="new_nested_dirs"),
    ],
)
def test_ensure_dir_pass(tmp_path: Path, parts: tuple[str, ...]) -> None:
    """Test that ensure_dir creates missing directories and returns the same Path object."""
    target = tmp_path.joinpath(*parts)
    result = ensure_dir(target)
    assert result is target
    assert target.is_dir()


@pytest.mark.parametrize(
    "parts",
    [
        pytest.param((), id="target_is_file"),
        pytest.param(("child",), id="parent_is_file"),
    ],
)
def test_ensure_dir_fail_os_error(tmp_path: Path, parts: tuple[str, ...]) -> None:
    """Test that OS errors from mkdir are re-raised as ValueError, with the OSError as the cause."""
    blocker = tmp_path / "blocker"
    blocker.write_text("not a directory")
    target = blocker.joinpath(*parts)

    with pytest.raises(ValueError, match="Cannot create directory") as exc_info:
        ensure_dir(target)

    assert str(exc_info.value).startswith(f"{DIR_NOT_CREATED}{target!r}: ")
    assert isinstance(exc_info.value.__cause__, OSError)
    assert blocker.is_file()


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        pytest.param("/", "/", id="root"),
        pytest.param("///", "/", id="multiple_slashes_only"),
        pytest.param("/a/b", "/a/b", id="no_trailing_slash"),
        pytest.param("/a/b/", "/a/b", id="one_trailing_slash"),
        pytest.param("/a/b///", "/a/b", id="many_trailing_slashes"),
        pytest.param("a/b/", "a/b", id="relative_path"),
        pytest.param(".", ".", id="current_dir"),
        pytest.param("/my dir/", "/my dir", id="inner_whitespace_kept"),
        pytest.param("file:///", "file:///", id="file_root"),
        pytest.param("file:///a/b/", "file:///a/b", id="file_protocol_trailing_slash"),
        pytest.param("file://a/b", "file://a/b", id="file_protocol_relative"),
        pytest.param("s3://bucket/", "s3://bucket", id="s3_bucket_trailing_slash"),
        pytest.param("s3://bucket/prefix//", "s3://bucket/prefix", id="s3_prefix_trailing_slashes"),
        pytest.param("s3:///", "s3:///", id="s3_root"),
        pytest.param("HTTPS://Example.com/", "HTTPS://Example.com", id="protocol_case_preserved"),
        pytest.param("data:abc/", "data:abc", id="data_uri"),
        pytest.param("data:/", "data:/", id="data_uri_root"),
        pytest.param("  /a/b/  ", "/a/b", id="whitespace_path_trailing_slash"),
        pytest.param("/a/b/\n", "/a/b", id="trailing_newline"),
        pytest.param("\t/\t", "/", id="whitespace_root"),
        pytest.param("  relative/dir/ ", "relative/dir", id="whitespace_relative_path"),
        pytest.param("  /my dir/  ", "/my dir", id="whitespace_inner_whitespace_kept"),
        pytest.param("  file:///a/b/  ", "file:///a/b", id="whitespace_file_protocol"),
        pytest.param(" file:/// ", "file:///", id="whitespace_file_root"),
        pytest.param("  s3://bucket/  ", "s3://bucket", id="whitespace_s3"),
        pytest.param("  data:abc/  ", "data:abc", id="whitespace_data_uri"),
    ],
)
def test_normalise_dir_pass(value: str, expected: str) -> None:
    """Test that whitespace and trailing slashes are removed while roots and protocols are kept."""
    assert normalise_dir(value) == expected


@pytest.mark.parametrize(
    "value",
    [
        pytest.param(None, id="none"),
        pytest.param("", id="empty_string"),
        pytest.param("   ", id="spaces_only"),
        pytest.param("\t\n", id="tab_newline_only"),
    ],
)
def test_normalise_dir_fail_no_content(value: str | None) -> None:
    """Test that None, empty and whitespace-only input raise a ValueError."""
    with pytest.raises(ValueError, match="at least one character") as exc_info:
        normalise_dir(value)
    assert str(exc_info.value) == VALID_PATH_ONE_CHAR


@pytest.mark.parametrize(
    "value",
    [
        pytest.param("file://", id="file"),
        pytest.param("https://", id="https"),
        pytest.param("s3://", id="s3"),
        pytest.param("data:", id="data"),
        pytest.param("  file://  ", id="file_with_whitespace"),
        pytest.param(" s3://\n", id="s3_with_whitespace"),
        pytest.param("  data:  ", id="data_with_whitespace"),
    ],
)
def test_normalise_dir_fail_protocol_only(value: str) -> None:
    """Test that a protocol with no path, with or without surrounding whitespace, raises a ValueError."""
    with pytest.raises(ValueError, match="Path must not be empty") as exc_info:
        normalise_dir(value)
    assert str(exc_info.value) == PATH_EMPTY


PATH_FORMS = [
    pytest.param(str, id="str"),
    pytest.param(lambda p: p, id="path_object"),
    pytest.param(lambda p: f"file://{p}", id="file_protocol"),
    pytest.param(lambda p: f"  {p}\n", id="str_with_whitespace"),
]


@pytest.mark.parametrize("make_value", PATH_FORMS)
def test_local_input_dir_pass_existing_dir(tmp_path: Path, make_value: Callable[[Path], Any]) -> None:
    """Test that an existing directory is accepted in each form and returned as an absolute Path."""
    model = InputDirModel.model_validate({"path": make_value(tmp_path)})
    assert model.path == tmp_path.resolve()
    assert model.path.is_absolute()


@pytest.mark.parametrize(
    ("value", "parts"),
    [
        pytest.param(".", (), id="current_dir"),
        pytest.param("sub", ("sub",), id="sub_dir"),
        pytest.param("file://sub", ("sub",), id="file_protocol_sub_dir"),
    ],
)
def test_local_input_dir_pass_relative_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, value: str, parts: tuple[str, ...]
) -> None:
    """Test that relative paths are resolved against the current working directory."""
    (tmp_path / "sub").mkdir()
    monkeypatch.chdir(tmp_path)
    model = InputDirModel.model_validate({"path": value})
    assert model.path == tmp_path.joinpath(*parts).resolve()
    assert model.path.is_absolute()


@pytest.mark.parametrize(
    "create_file",
    [
        pytest.param(False, id="does_not_exist"),
        pytest.param(True, id="is_a_file"),
    ],
)
def test_local_input_dir_fail_not_a_directory(tmp_path: Path, create_file: bool) -> None:
    """Test that a missing path or a file is rejected by LocalInputDir."""
    target = tmp_path / "target"
    if create_file:
        target.write_text("content")

    with pytest.raises(ValidationError) as exc_info:
        InputDirModel.model_validate({"path": str(target)})

    assert get_single_error(exc_info.value)["type"] == "path_not_directory"
    assert target.is_file() is create_file


@pytest.mark.parametrize("make_value", PATH_FORMS)
def test_local_output_dir_pass_creates_dir(tmp_path: Path, make_value: Callable[[Path], Any]) -> None:
    """Test that a missing nested directory is created in each input form and returned as an absolute Path."""
    target = tmp_path / "out" / "nested"
    model = OutputDirModel.model_validate({"path": make_value(target)})
    assert model.path == target.resolve()
    assert model.path.is_absolute()
    assert target.is_dir()


def test_local_output_dir_pass_existing_dir(tmp_path: Path) -> None:
    """Test that an existing directory, with its contents, is accepted as is."""
    existing_file = tmp_path / "keep.txt"
    existing_file.write_text("keep")
    model = OutputDirModel.model_validate({"path": str(tmp_path)})
    assert model.path == tmp_path.resolve()
    assert existing_file.read_text() == "keep"


def test_local_output_dir_pass_relative_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test that a relative path is created under the current working directory."""
    monkeypatch.chdir(tmp_path)
    model = OutputDirModel.model_validate({"path": "out/nested"})
    expected = (tmp_path / "out" / "nested").resolve()
    assert model.path == expected
    assert expected.is_dir()


def test_local_output_dir_fail_cannot_create(tmp_path: Path) -> None:
    """Test that a failure to create the directory is reported as a pydantic value error."""
    blocker = tmp_path / "blocker"
    blocker.write_text("not a directory")

    with pytest.raises(ValidationError) as exc_info:
        OutputDirModel.model_validate({"path": str(blocker)})

    error = get_single_error(exc_info.value)
    assert error["type"] == "value_error"
    assert error["msg"].startswith(f"Value error, {DIR_NOT_CREATED}{blocker.resolve()!r}: ")
    assert blocker.is_file()


@pytest.mark.parametrize("model_cls", BOTH_MODELS)
@pytest.mark.parametrize(
    "value",
    [
        pytest.param("https://example.com/data", id="https"),
        pytest.param("s3://bucket/key", id="s3"),
        pytest.param("data:text/plain,hello", id="data_uri"),
        pytest.param("  gs://bucket  ", id="gcs_with_whitespace"),
    ],
)
def test_local_dir_models_fail_not_local(model_cls: type[BaseModel], value: str) -> None:
    """Test that both local dir types reject remote URLs."""
    with pytest.raises(ValidationError) as exc_info:
        model_cls.model_validate({"path": value})

    error = get_single_error(exc_info.value)
    assert error["type"] == "value_error"
    assert error["msg"] == f"Value error, {PATH_NOT_LOCAL}: {value.strip()!r}"


@pytest.mark.parametrize("model_cls", BOTH_MODELS)
@pytest.mark.parametrize(
    "value",
    [
        pytest.param("", id="empty_string"),
        pytest.param("   ", id="whitespace_only"),
        pytest.param("file://", id="file_protocol_only"),
    ],
)
def test_local_dir_models_fail_empty(model_cls: type[BaseModel], value: str) -> None:
    """Test that both local dir types reject empty paths instead of treating them as the current directory."""
    with pytest.raises(ValidationError) as exc_info:
        model_cls.model_validate({"path": value})

    error = get_single_error(exc_info.value)
    assert error["type"] == "value_error"
    assert error["msg"] == f"Value error, {PATH_EMPTY}"


@pytest.mark.parametrize("model_cls", BOTH_MODELS)
@pytest.mark.parametrize(
    "value",
    [
        pytest.param(None, id="none"),
        pytest.param(123, id="int"),
        pytest.param(["/tmp"], id="list"),
    ],
)
def test_local_dir_models_fail_invalid_type(model_cls: type[BaseModel], value: Any) -> None:  # noqa: ANN401
    """Test that values that are not paths pass the before validator and are rejected by pydantic."""
    with pytest.raises(ValidationError) as exc_info:
        model_cls.model_validate({"path": value})

    assert get_single_error(exc_info.value)["type"] == "path_type"
