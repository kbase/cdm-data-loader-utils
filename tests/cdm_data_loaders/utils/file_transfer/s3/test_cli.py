"""Tests for s3_utils.py using moto to mock AWS S3."""

import json
import re
from pathlib import Path
from unittest.mock import call, patch

import pytest
from types_boto3_s3.client import S3Client

from cdm_data_loaders.utils.file_transfer.s3.cli import (
    cmd_cp,
    cmd_head,
    cmd_ls,
    cmd_mb,
)
from tests.cdm_data_loaders.utils.file_transfer.s3.conftest import (
    COULD_NOT_PARSE,
    NO_PATH_FOUND,
    SAMPLE_FILES,
    START_WITH_BUCKET_NAME,
)

# CLI helper function tests

TEST_PROTOCOLS = [
    "",
    "s3://",
    "s3a://",
]

TEST_MB_PARAMS = [
    pytest.param("one-bucket", "one-bucket", id="name only"),
    pytest.param("two-bucket/", "two-bucket", id="trailing slash"),
    pytest.param("red-bucket/key", "red-bucket", id="red-bucket"),
    pytest.param("blue-bucket/foo/bar.txt", "blue-bucket", id="blue-bucket"),
]

ERROR_MB_USAGE = "Usage: s3_local.py mb s3://BUCKET"

TEST_MB_ARGS_ERROR = [
    pytest.param([], ERROR_MB_USAGE, id="no args"),
    pytest.param(["foo", "bar"], ERROR_MB_USAGE, id="two args"),
    pytest.param(["foo", "bar", "baz"], ERROR_MB_USAGE, id="three args"),
    pytest.param([""], NO_PATH_FOUND, id="missing bucket"),
    pytest.param(["/bucket"], START_WITH_BUCKET_NAME, id="preceding slash"),
    pytest.param(["/"], START_WITH_BUCKET_NAME, id="root"),
]


@pytest.mark.s3
@pytest.mark.parametrize("protocol", TEST_PROTOCOLS)
@pytest.mark.parametrize(("path", "bucket"), TEST_MB_PARAMS)
def test_cmd_mb_creates_bucket(mock_s3_client: S3Client, protocol: str, path: str, bucket: str) -> None:
    """CLI mb helper creates bucket for valid paths."""
    with patch("builtins.print") as mock_print:
        cmd_mb([f"{protocol}{path}"])
    mock_print.assert_called_once_with(f"Created bucket: {bucket}")
    mock_s3_client.head_bucket(Bucket=bucket)


@pytest.mark.s3
@pytest.mark.parametrize("protocol", TEST_PROTOCOLS)
@pytest.mark.parametrize(("path", "bucket"), TEST_MB_PARAMS)
def test_cmd_mb_handles_existing_bucket(mock_s3_client: S3Client, protocol: str, path: str, bucket: str) -> None:
    """CLI mb helper prints message when bucket already exists."""
    args = [f"{protocol}{path}"]
    with patch("builtins.print") as mock_print:
        cmd_mb(args)
        cmd_mb(args)
    mock_print.assert_has_calls(
        [
            call(f"Created bucket: {bucket}"),
            call(f"Bucket already exists: {bucket}"),
        ]
    )
    mock_s3_client.head_bucket(Bucket=bucket)


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize(("args", "err_msg"), TEST_MB_ARGS_ERROR)
def test_cmd_mb_prints_error(args: list[str], err_msg: str) -> None:
    """CLI mb helper prints usage on invalid argument list."""
    with pytest.raises(SystemExit, match=err_msg):
        cmd_mb(args)


TEST_CP_PATHS = [
    pytest.param("bucket/path/to/file.txt", "bucket", "path/to/file.txt", id="no prefix"),
    pytest.param("s3://foo/bar/baz.qux", "foo", "bar/baz.qux", id="with prefix"),
    pytest.param("bucket/root.txt", "bucket", "root.txt", id="no prefix root file"),
    pytest.param("s3://foo/bar.baz", "foo", "bar.baz", id="with prefix root file"),
]


@pytest.mark.s3
@pytest.mark.parametrize(("path", "bucket", "key"), TEST_CP_PATHS)
def test_cmd_cp_copies_file(mock_s3_client: S3Client, sample_file: Path, path: str, bucket: str, key: str) -> None:
    """CLI cp helper copies local file to store."""
    mock_s3_client.create_bucket(Bucket=bucket)
    with patch("builtins.print") as mock_print:
        cmd_cp([str(sample_file), path])
    mock_print.assert_has_calls(
        [
            call(f"  {key}"),
            call(f"Uploaded 1 files to s3://{bucket}/{key}"),
        ]
    )
    obj = mock_s3_client.get_object(Bucket=bucket, Key=key)
    assert obj["Body"].read() == b"hello s3"


TEST_CP_KEY_PREFIXES = [
    pytest.param("bucket/path/to", "bucket", "path/to/", id="no prefix"),
    pytest.param("bucket/path/to/", "bucket", "path/to/", id="no prefix trailing slash"),
    pytest.param("s3://foo/bar", "foo", "bar/", id="with prefix"),
    pytest.param("s3://foo/bar/", "foo", "bar/", id="with prefix trailing slash"),
    pytest.param("bucket", "bucket", "", id="no prefix root"),
    pytest.param("bucket/", "bucket", "", id="no prefix root trailing slash"),
    pytest.param("s3://foo", "foo", "", id="with prefix root"),
    pytest.param("s3://foo/", "foo", "", id="with prefix root trailing slash"),
]


@pytest.mark.s3
@pytest.mark.parametrize(("path", "bucket", "prefix"), TEST_CP_KEY_PREFIXES)
def test_cmd_cp_copies_dir(mock_s3_client: S3Client, sample_dir: Path, path: str, bucket: str, prefix: str) -> None:
    """CLI cp helper copies local folder to store."""
    mock_s3_client.create_bucket(Bucket=bucket)
    with patch("builtins.print") as mock_print:
        cmd_cp([str(sample_dir), path])
    keys = {obj["Key"] for obj in mock_s3_client.list_objects_v2(Bucket=bucket)["Contents"]}
    expected_paths = [Path(prefix) / rel_path if prefix else Path(rel_path) for rel_path in SAMPLE_FILES]
    mock_print.assert_has_calls(
        [
            *[call(f"  {key}") for key in expected_paths],
            call(f"Uploaded {len(keys)} files to s3://{Path(bucket) / prefix}/"),
        ]
    )
    assert keys == {f"{Path(prefix) / f}" for f in SAMPLE_FILES}


ERROR_CP_USAGE = "Usage: s3_local.py cp [LOCAL_DIR | LOCAL_FILE] s3://BUCKET[/PREFIX/]"

TEST_CP_ARGS_ERROR = [
    pytest.param([], ERROR_CP_USAGE, id="no args"),
    pytest.param(["foo"], ERROR_CP_USAGE, id="one arg"),
    pytest.param(["foo", "bar", "baz"], ERROR_CP_USAGE, id="three args"),
    pytest.param(["local/file.txt", ""], NO_PATH_FOUND, id="local file; missing bucket"),
    pytest.param(["local/folder", ""], NO_PATH_FOUND, id="local folder; missing bucket"),
    pytest.param(["local/folder/", ""], NO_PATH_FOUND, id="local folder/;missing bucket"),
    pytest.param(["local/file.txt", "/bucket"], START_WITH_BUCKET_NAME, id="local file; preceding slash"),
    pytest.param(["local/folder", "/bucket"], START_WITH_BUCKET_NAME, id="local folder; preceding slash"),
    pytest.param(["local/folder/", "/bucket"], START_WITH_BUCKET_NAME, id="local folder/; preceding slash"),
    pytest.param(["local/file.txt", "/"], START_WITH_BUCKET_NAME, id="local file; root"),
    pytest.param(["local/folder", "/"], START_WITH_BUCKET_NAME, id="local folder; root"),
    pytest.param(["local/folder/", "/"], START_WITH_BUCKET_NAME, id="local folder/; root"),
]


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize(("args", "err_msg"), TEST_CP_ARGS_ERROR)
def test_cmd_cp_prints_error(args: list[str], err_msg: str) -> None:
    """CLI cp helper prints usage on invalid argument list."""
    with pytest.raises(SystemExit, match=re.escape(err_msg)):
        cmd_cp(args)


LS_FILES = [
    "key/prefix/file1.txt",
    "key/prefix/file2.txt",
    "key/prefix/sub/file3.txt",
]

TEST_LS_KEY_PREFIXES = [
    pytest.param("some-bucket/key/prefix/", "some-bucket", id="trailing slash"),
    pytest.param("other-bucket/key/prefix", "other-bucket", id="no trailing slash"),
    pytest.param("this-bucket/", "this-bucket", id="bucket trailing slash"),
    pytest.param("that-bucket", "that-bucket", id="bucket only"),
    pytest.param("s3://s3-bucket/key/prefix/", "s3-bucket", id="s3 protocol trailing slash"),
    pytest.param("s3a://s3a-bucket/key/prefix/", "s3a-bucket", id="s3a protocol trailing slash"),
]

ERROR_LS_USAGE = "Usage: s3_local.py ls s3://BUCKET[/PREFIX/] [--limit N]"

TEST_LS_ARGS_ERROR = [
    pytest.param([], ERROR_LS_USAGE, id="no args"),
    pytest.param([""], NO_PATH_FOUND, id="empty path"),
    pytest.param(["/bucket"], START_WITH_BUCKET_NAME, id="preceding slash"),
    pytest.param(["/"], START_WITH_BUCKET_NAME, id="root"),
]


@pytest.mark.s3
@pytest.mark.parametrize(("prefix", "bucket"), TEST_LS_KEY_PREFIXES)
def test_cmd_ls_lists_keys(mock_s3_client: S3Client, prefix: str, bucket: str) -> None:
    """CLI ls helper lists keys under the given prefix."""
    mock_s3_client.create_bucket(Bucket=bucket)
    for key in LS_FILES:
        mock_s3_client.put_object(Bucket=bucket, Key=key, Body=b"a")
    with patch("builtins.print") as mock_print:
        cmd_ls([prefix])
    printed_keys = {call.args[0].split()[-1] for call in mock_print.call_args_list}
    assert printed_keys == set(LS_FILES)


@pytest.mark.s3
def test_cmd_ls_respects_limit(mock_s3_client: S3Client) -> None:
    """CLI ls helper stops after --limit objects."""
    mock_s3_client.create_bucket(Bucket="some-bucket")
    for i in range(5):
        mock_s3_client.put_object(Bucket="some-bucket", Key=f"key/file{i}.txt", Body=b"x")

    with patch("builtins.print") as mock_print:
        cmd_ls(["some-bucket/key", "--limit", "2"])

    assert mock_print.call_count == 2  # noqa: PLR2004


@pytest.mark.s3
def test_cmd_ls_default_limit_is_20(mock_s3_client: S3Client) -> None:
    """CLI ls helper defaults to a limit of 20 objects."""
    bucket = "large-bucket"
    mock_s3_client.create_bucket(Bucket=bucket)
    for i in range(25):
        mock_s3_client.put_object(Bucket=bucket, Key=f"key/file{i:02d}.txt", Body=b"x")

    with patch("builtins.print") as mock_print:
        cmd_ls([f"{bucket}/key"])

    assert mock_print.call_count == 20  # noqa: PLR2004


@pytest.mark.s3
def test_cmd_ls_empty_prefix_returns_nothing(mock_s3_client: S3Client) -> None:
    """CLI ls helper prints nothing when no objects match the prefix."""
    bucket = "foo-bucket"
    mock_s3_client.create_bucket(Bucket=bucket)
    with patch("builtins.print") as mock_print:
        cmd_ls([f"{bucket}/nonexistent/"])

    mock_print.assert_not_called()


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize(("args", "err_msg"), TEST_LS_ARGS_ERROR)
def test_cmd_ls_prints_error(args: list[str], err_msg: str) -> None:
    """CLI ls helper raises SystemExit with a useful message on invalid arguments."""
    with pytest.raises(SystemExit, match=re.escape(err_msg)):
        cmd_ls(args)


# cmd_head

ERROR_HEAD_USAGE = "Usage: s3_local.py head s3://BUCKET/KEY"

TEST_HEAD_ARGS_ERROR = [
    pytest.param([], ERROR_HEAD_USAGE, id="no args"),
    pytest.param([""], NO_PATH_FOUND, id="empty path"),
    pytest.param(["/bucket/key"], START_WITH_BUCKET_NAME, id="preceding slash"),
    pytest.param(["bucket"], COULD_NOT_PARSE, id="bucket only"),
    pytest.param(["bucket/"], COULD_NOT_PARSE, id="bucket trailing slash"),
]


@pytest.mark.s3
def test_cmd_head_prints_metadata(mock_s3_client: S3Client) -> None:
    """CLI head helper prints metadata for an existing object."""
    bucket = "bar-bucket"
    mock_s3_client.create_bucket(Bucket=bucket)
    mock_s3_client.put_object(
        Bucket=bucket,
        Key="some/file.txt",
        Body=b"data",
        Metadata={"md5": "abc123", "source": "ncbi"},
    )
    with patch("builtins.print") as mock_print:
        cmd_head([f"{bucket}/some/file.txt"])
    mock_print.assert_any_call(f"Metadata for {bucket}/some/file.txt:")
    # second call is the JSON dump — check both keys are present
    json_output = mock_print.call_args_list[1].args[0]
    parsed = json.loads(json_output)
    assert parsed["md5"] == "abc123"
    assert parsed["source"] == "ncbi"


@pytest.mark.s3
def test_cmd_head_prints_empty_metadata(mock_s3_client: S3Client) -> None:
    """CLI head helper prints empty JSON object when no metadata is set."""
    bucket = "baz-bucket"
    mock_s3_client.create_bucket(Bucket=bucket)
    mock_s3_client.put_object(Bucket=bucket, Key="bare/file.txt", Body=b"data")
    with patch("builtins.print") as mock_print:
        cmd_head([f"{bucket}/bare/file.txt"])

    mock_print.assert_any_call(f"Metadata for {bucket}/bare/file.txt:")
    json_output = mock_print.call_args_list[1].args[0]
    assert json.loads(json_output) == {}


@pytest.mark.s3
def test_cmd_head_prints_not_found_for_missing_key(mock_s3_client: S3Client) -> None:
    """CLI head helper prints a not-found message instead of raising for a 404."""
    bucket = "empty-bucket"
    mock_s3_client.create_bucket(Bucket=bucket)
    with patch("builtins.print") as mock_print:
        cmd_head([f"{bucket}/does/not/exist.txt"])

    mock_print.assert_called_once_with(f"File not found in store: {bucket}/does/not/exist.txt")


@pytest.mark.s3
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_cmd_head_handles_protocols(mock_s3_client: S3Client, protocol: str) -> None:
    """CLI head helper accepts all valid S3 protocol prefixes."""
    bucket = "proto-bucket"
    mock_s3_client.create_bucket(Bucket=bucket)
    mock_s3_client.put_object(Bucket=bucket, Key="proto/file.txt", Body=b"x")

    with patch("builtins.print") as mock_print:
        cmd_head([f"{protocol}{bucket}/proto/file.txt"])

    mock_print.assert_any_call(f"Metadata for {bucket}/proto/file.txt:")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize(("args", "err_msg"), TEST_HEAD_ARGS_ERROR)
def test_cmd_head_prints_error(args: list[str], err_msg: str) -> None:
    """CLI head helper raises SystemExit with a useful message on invalid arguments."""
    with pytest.raises(SystemExit, match=err_msg):
        cmd_head(args)
