"""Tests for s3 utils.py using moto to mock AWS S3."""

from pathlib import Path
from typing import Final

import pytest
from types_boto3_s3.client import S3Client

import cdm_data_loaders.utils.file_transfer.s3.client as s3_client
from cdm_data_loaders.utils.file_transfer.s3.client import reset_s3_client
from tests.conftest import ALT_BUCKET, TEST_BUCKET

HTTP_200: Final[int] = 200  # Status OK
HTTP_204: Final[int] = 204  # Status no content

SIZE_HELLO: Final[int] = 5
SIZE_DATA: Final[int] = 4

SAMPLE_FILES = [
    "dir_one/file1.txt",
    "dir_one/file2.txt",
    "dir_one/sub_dir/file3.txt",
    "dir_one/sub_dir/under_dir/file4.txt",
]

FILES_IN_BUCKETS = {
    TEST_BUCKET: SAMPLE_FILES,
    ALT_BUCKET: ["dir_one/file1.txt"],
}


# CLI helper function tests

NO_PATH_FOUND: Final[str] = "No path found"
START_WITH_BUCKET_NAME: Final[str] = "s3 paths must start with the bucket name"
COULD_NOT_PARSE: Final[str] = "Could not parse out bucket and key"

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

ERROR_MB_USAGE: Final[str] = "Usage: s3_local.py mb s3://BUCKET"

TEST_MB_ARGS_ERROR = [
    pytest.param([], ERROR_MB_USAGE, id="no args"),
    pytest.param(["foo", "bar"], ERROR_MB_USAGE, id="two args"),
    pytest.param(["foo", "bar", "baz"], ERROR_MB_USAGE, id="three args"),
    pytest.param([""], NO_PATH_FOUND, id="missing bucket"),
    pytest.param(["/bucket"], START_WITH_BUCKET_NAME, id="preceding slash"),
    pytest.param(["/"], START_WITH_BUCKET_NAME, id="root"),
]


@pytest.fixture
def sample_file(tmp_path: Path) -> Path:
    """Create a small temporary file for upload tests."""
    f = tmp_path / "sample.txt"
    f.write_text("hello s3")
    return f


@pytest.fixture
def sample_dir(tmp_path: Path) -> Path:
    """Create a small temporary directory tree for upload_dir tests.

    Structure (same as TEST_BUCKET files)

    dir_one/file1.txt
    dir_one/file2.txt
    dir_one/sub_dir/file3.txt
    dir_one/sub_dir/under_dir/file4.txt
    """
    sample_dir = tmp_path / "sample_dir"
    for f in SAMPLE_FILES:
        new_file = sample_dir / f
        # ensure the parent dir exists
        new_file.parent.mkdir(parents=True, exist_ok=True)
        # add the (relative) path as the content
        new_file.write_text(f)
    return sample_dir


def populate_mock_s3(client: S3Client, file_list_by_bucket: dict[str, list[str]]) -> None:
    """Populate buckets with a list of files.

    File names should be a list, indexed by bucket.

    Files will be populated with the file name as bytes if the top level directory is `dir_one`;
    otherwise, the content will just be `x`.

    :param client: s3 client
    :type client: S3Client
    :param file_list: list of files, indexed by bucket
    :type file_list: dict[str, list[str]]
    """
    for bucket, file_list in file_list_by_bucket.items():
        for file in file_list:
            full_path = f"{bucket}/{file}"
            if file.startswith("dir_one"):
                client.put_object(Bucket=bucket, Key=file, Body=full_path.encode("utf-8"))
            else:
                client.put_object(Bucket=bucket, Key=file, Body=b"x")
            # if this errors, the transfer was not successful
            client.head_object(Bucket=bucket, Key=file)


def prep_client_init(monkeypatch: pytest.MonkeyPatch) -> None:
    """Set up environment variables to allow get_s3_client to initialize without error."""
    reset_s3_client()
    assert s3_client._s3_client is None  # noqa: SLF001
    # set up env vars to ensure that the argument takes precedence
    monkeypatch.setenv("AWS_ENDPOINT_URL", "http://env-endpoint.com")
    monkeypatch.delenv("AWS_ENDPOINT_URL_S3", raising=False)
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "aws_access_key_id_env_var")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "aws_secret_access_key_env_var")
