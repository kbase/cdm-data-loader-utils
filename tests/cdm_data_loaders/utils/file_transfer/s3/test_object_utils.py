"""Tests for s3_object utils.py using moto to mock AWS S3."""

import io
import logging
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch

import pytest
from boto3.s3.transfer import TransferConfig
from botocore.exceptions import ClientError
from types_boto3_s3.client import S3Client

from cdm_data_loaders.utils.file_transfer.s3 import object_utils
from cdm_data_loaders.utils.file_transfer.s3.object_utils import (
    DEFAULT_EXTRA_ARGS,
    copy_directory,
    copy_object,
    delete_object,
    delete_objects,
    download_file,
    get_existing_object_info,
    head_object,
    list_objects,
    object_exists,
    split_s3_path,
    upload_dir,
    upload_file,
    upload_fileobj,
)
from tests.cdm_data_loaders.utils.file_transfer.s3.conftest import (
    ALT_BUCKET,
    FILES_IN_BUCKETS,
    HTTP_200,
    HTTP_204,
    SAMPLE_FILES,
    SIZE_DATA,
    SIZE_HELLO,
    TEST_BUCKET,
    populate_mock_s3,
)
from tests.conftest import BUCKETS

# split_s3_path
PATH = "path"
TO = "to"
TO_FILE = "to/file.txt"
PATH_TO_FILE = f"{PATH}/{TO_FILE}"

EXPECTED = {
    "path/to": (PATH, TO),
    "path/to/": (PATH, "to/"),
    "path/to/file.txt": (PATH, TO_FILE),
    "s3://path/to": (PATH, TO),
    "s3://path/to/": (PATH, "to/"),
    "s3://path/to/file.txt": (PATH, TO_FILE),
    "s3a://path/to": (PATH, TO),
    "s3a://path/to/": (PATH, "to/"),
    "s3a://path/to/file.txt": (PATH, TO_FILE),
}

NO_PATH_FOUND = "No path found"
START_WITH_BUCKET_NAME = "s3 paths must start with the bucket name"
COULD_NOT_PARSE = "Could not parse out bucket and key"

INVALID_PATH_ERRORS = {
    "": NO_PATH_FOUND,
    "/": START_WITH_BUCKET_NAME,
    "/path": START_WITH_BUCKET_NAME,
    "/path/to/file.txt": START_WITH_BUCKET_NAME,
    "path": COULD_NOT_PARSE,
    "path/": COULD_NOT_PARSE,
    "s3://": NO_PATH_FOUND,
    "s3:///": START_WITH_BUCKET_NAME,
    "s3://path": COULD_NOT_PARSE,
    "s3://path/": COULD_NOT_PARSE,
    "s3a://": NO_PATH_FOUND,
    "s3a://path": COULD_NOT_PARSE,
    "s3a://path/": COULD_NOT_PARSE,
}


def make_mock_requests(
    content: bytes = b"hello world",
    status_code: int = 200,
    content_type: str = "application/octet-stream",
) -> tuple[MagicMock, MagicMock]:
    """Build a mock requests module whose .get() returns a mock response."""
    mock_response = MagicMock()
    mock_response.status_code = status_code
    mock_response.raw = io.BytesIO(content)
    mock_response.raw.decode_content = True
    mock_response.headers = {
        "content-type": content_type,
    }
    mock_response.raise_for_status = MagicMock()
    mock_response.__enter__ = lambda s: s
    mock_response.__exit__ = MagicMock(return_value=False)

    mock_requests = MagicMock()
    mock_requests.get.return_value = mock_response

    return mock_requests, mock_response


@pytest.mark.parametrize("invalid_path", list(INVALID_PATH_ERRORS.keys()))
def test_split_s3_path_errors(invalid_path: str) -> None:
    """Ensure that an error is thrown if an invalid s3 path is passed in."""
    with pytest.raises(ValueError, match=INVALID_PATH_ERRORS[invalid_path]):
        split_s3_path(invalid_path)


@pytest.mark.parametrize("valid_path", list(EXPECTED.keys()))
def test_split_s3_path_success(valid_path: str) -> None:
    """Verify that a valid path is correctly split into bucket and key."""
    (bucket, path) = split_s3_path(valid_path)
    assert (bucket, path) == EXPECTED[valid_path]


EXPECTED_ALLOW_BUCKET_ONLY = {
    "path/to": (PATH, TO),
    "path/to/": (PATH, "to/"),
    "path/to/file.txt": (PATH, TO_FILE),
    "path": (PATH, ""),
    "path/": (PATH, ""),
    "s3://path/to": (PATH, TO),
    "s3://path/to/": (PATH, "to/"),
    "s3://path/to/file.txt": (PATH, TO_FILE),
    "s3://path": (PATH, ""),
    "s3://path/": (PATH, ""),
    "s3a://path/to": (PATH, TO),
    "s3a://path/to/": (PATH, "to/"),
    "s3a://path/to/file.txt": (PATH, TO_FILE),
    "s3a://path": (PATH, ""),
    "s3a://path/": (PATH, ""),
}

INVALID_PATH_ALLOW_BUCKET_ONLY_ERRORS = {
    "": NO_PATH_FOUND,
    "/": START_WITH_BUCKET_NAME,
    "/path": START_WITH_BUCKET_NAME,
    "/path/to/file.txt": START_WITH_BUCKET_NAME,
    "s3://": NO_PATH_FOUND,
    "s3:///": START_WITH_BUCKET_NAME,
    "s3a://": NO_PATH_FOUND,
}


@pytest.mark.parametrize("invalid_path", list(INVALID_PATH_ALLOW_BUCKET_ONLY_ERRORS.keys()))
def test_split_s3_path_allow_bucket_only_errors(invalid_path: str) -> None:
    """Ensure that an error is thrown if an invalid s3 path is passed in."""
    with pytest.raises(ValueError, match=INVALID_PATH_ALLOW_BUCKET_ONLY_ERRORS[invalid_path]):
        split_s3_path(invalid_path, allow_bucket_only=True)


@pytest.mark.parametrize("valid_path", list(EXPECTED_ALLOW_BUCKET_ONLY.keys()))
def test_split_s3_path_allow_bucket_only_success(valid_path: str) -> None:
    """Verify that a valid path is correctly split into bucket and key."""
    (bucket, path) = split_s3_path(valid_path, allow_bucket_only=True)
    assert (bucket, path) == EXPECTED_ALLOW_BUCKET_ONLY[valid_path]


# head_object
@pytest.mark.s3
def test_head_object_returns_info(mock_s3_client: S3Client) -> None:
    """Verify that head_object returns size, metadata, and checksum fields."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key="info/file.txt", Body=b"hello", Metadata={"md5": "abc123"})
    result = head_object(f"{TEST_BUCKET}/info/file.txt")
    assert result is not None
    assert result["ContentLength"] == SIZE_HELLO
    assert result["Metadata"]["md5"] == "abc123"


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_head_object_raises_for_missing() -> None:
    """Verify that head_object returns None for a non-existent object."""
    with pytest.raises(ClientError, match="404"):
        head_object(f"{TEST_BUCKET}/does/not/exist.txt")


@pytest.mark.s3
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_head_object_with_protocols(mock_s3_client: S3Client, protocol: str) -> None:
    """Verify that head_object handles all valid protocol prefixes."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key="proto/file.txt", Body=b"data")
    result = head_object(f"{protocol}{TEST_BUCKET}/proto/file.txt")
    assert result is not None
    assert result["ContentLength"] == SIZE_DATA


# get_existing_object_info
@pytest.mark.s3
@pytest.mark.parametrize("metadata", [{"md5": "abc123"}, {}, None])
def test_get_existing_object_info_pass_returns_info_for_existing_object(
    mock_s3_client: S3Client, metadata: dict[str, str] | None
) -> None:
    """Verify size, etag, and lowercased metadata are all populated for an existing object."""
    kwargs = {"Bucket": TEST_BUCKET, "Key": "info/file.txt", "Body": b"hello"}
    if metadata is not None:
        kwargs["Metadata"] = metadata

    mock_s3_client.put_object(**kwargs)
    info = get_existing_object_info(f"{TEST_BUCKET}/info/file.txt")
    assert info is not None
    assert info["ContentLength"] == SIZE_HELLO
    assert info["ETag"]
    assert info["Metadata"] == (metadata or {})


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_get_existing_object_info_pass_returns_none_for_missing_object() -> None:
    """Verify None is returned (not an exception) for a missing object."""
    assert get_existing_object_info(f"{TEST_BUCKET}does/not/exist.txt") is None


@pytest.mark.s3
@pytest.mark.parametrize("error_code", ["404", "NoSuchKey", "NotFound", "NoSuchBucket"])
def test_get_existing_object_info_pass_recognizes_all_not_found_error_codes(
    mock_s3_client: S3Client, monkeypatch: pytest.MonkeyPatch, error_code: str
) -> None:
    """NOT_FOUND_ERROR_CODES are treated as 'not found', even with a non-404 HTTP status."""

    def raise_not_found(**kwargs: Any) -> None:
        raise ClientError(
            {"Error": {"Code": error_code}, "ResponseMetadata": {"HTTPStatusCode": 400}},
            "HeadObject",
        )

    monkeypatch.setattr(mock_s3_client, "head_object", raise_not_found)
    assert get_existing_object_info(f"{TEST_BUCKET}some/key.txt") is None


@pytest.mark.s3
def test_get_existing_object_info_pass_recognizes_404_status_with_unknown_error_code(
    mock_s3_client: S3Client, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A 404 HTTP status is treated as 'not found' even if the error code itself isn't recognized."""

    def raise_404_status(**kwargs: Any) -> None:
        raise ClientError(
            {"Error": {"Code": "SomeWeirdCode"}, "ResponseMetadata": {"HTTPStatusCode": 404}},
            "HeadObject",
        )

    monkeypatch.setattr(mock_s3_client, "head_object", raise_404_status)
    assert get_existing_object_info(f"{TEST_BUCKET}some/key.txt") is None


@pytest.mark.s3
def test_get_existing_object_info_fail_non_not_found_error_propagates(
    mock_s3_client: S3Client, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A ClientError that is neither a recognized not-found code nor a 404 status is re-raised."""

    def raise_access_denied(**kwargs: Any) -> None:
        raise ClientError(
            {"Error": {"Code": "AccessDenied"}, "ResponseMetadata": {"HTTPStatusCode": 403}},
            "HeadObject",
        )

    monkeypatch.setattr(mock_s3_client, "head_object", raise_access_denied)
    with pytest.raises(ClientError, match="AccessDenied"):
        get_existing_object_info(f"{TEST_BUCKET}some/key.txt")


# object_exists
@pytest.mark.s3
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_head_object_and_object_exists_true_and_false(mock_s3_client: S3Client, protocol: str) -> None:
    """Verify that object_exists returns True for an object that exists in the bucket."""
    populate_mock_s3(mock_s3_client, FILES_IN_BUCKETS)
    for bucket, file_list in FILES_IN_BUCKETS.items():
        for f in file_list:
            output = head_object(f"{protocol}{bucket}/{f}")
            assert output is not None
            assert isinstance(output["ContentLength"], int)
            assert object_exists(f"{protocol}{bucket}/{f}") is True

        nonexistent_file = f"{protocol}{bucket}/a-file-i-just-made-up.txt"
        assert object_exists(nonexistent_file) is False
        with pytest.raises(ClientError, match="404"):
            head_object(nonexistent_file)


@pytest.mark.s3
@pytest.mark.parametrize("s3_path", ["absent", "dir_one", "dir_one/", "dir_one/file1.tnt"])
@pytest.mark.parametrize("bucket", BUCKETS)
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_object_exists_returns_false_when_absent(
    mock_s3_client: S3Client, s3_path: str, protocol: str, bucket: str
) -> None:
    """Verify that object_exists returns False for an object that does not exist."""
    populate_mock_s3(mock_s3_client, FILES_IN_BUCKETS)
    assert object_exists(f"{protocol}{bucket}/{s3_path}") is False


# list_objects
@pytest.mark.s3
@pytest.mark.parametrize("bucket", BUCKETS)
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_list_objects_lists_objects(
    mock_s3_client: S3Client,
    bucket: str,
    protocol: str,
) -> None:
    """Verify that all objects under a given prefix are returned, regardless of whether the protocol is supplied."""
    populate_mock_s3(mock_s3_client, FILES_IN_BUCKETS)
    contents = list_objects(f"{protocol}{bucket}/dir_one")
    keys = {obj["Key"] for obj in contents}
    assert keys == {f for f in FILES_IN_BUCKETS[bucket] if f.startswith("dir_one")}


@pytest.mark.s3
@pytest.mark.parametrize("dir_path", ["dir_one/sub_dir", "dir_one/sub_dir/", "dir_one/sub_dir/und"])
def test_list_objects_filters_by_prefix(
    mock_s3_client: S3Client,
    dir_path: str,
) -> None:
    """Check that more specific queries, including those that have 'incomplete' dir/file names, return correct results."""
    bucket = TEST_BUCKET
    populate_mock_s3(mock_s3_client, {bucket: FILES_IN_BUCKETS[bucket]})
    contents = list_objects(f"{bucket}/{dir_path}")
    keys = {obj["Key"] for obj in contents}
    # make sure this is a subset of all the files in the bucket
    assert len(keys) < len(FILES_IN_BUCKETS[bucket])
    assert keys == {f for f in FILES_IN_BUCKETS[bucket] if f.startswith(dir_path)}


@pytest.mark.s3
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_list_objects_empty_for_missing_prefix(
    mock_s3_client: S3Client,
    protocol: str,
) -> None:
    """Verify that an empty list is returned when no objects match the given prefix."""
    populate_mock_s3(mock_s3_client, FILES_IN_BUCKETS)
    for bucket in FILES_IN_BUCKETS:
        contents = list_objects(f"{protocol}{bucket}/nonexistent/")
        assert contents == []


N_FILES = 1005
DIR_TWO_FILES = [f"dir_two/file_{i:04d}.txt" for i in range(N_FILES)]
DIRTY_DATA = [f"dirty_data/file_{i:04d}.txt" for i in range(N_FILES)]
# pagination tests (1005 objects each, to exceed the 1000-item S3 page limit)
LOTS_OF_FILES = {
    TEST_BUCKET: [
        *DIR_TWO_FILES,
        *DIRTY_DATA,
    ]
}

EXPECTED_FILE_LIST = {
    "di": [*FILES_IN_BUCKETS[TEST_BUCKET], *LOTS_OF_FILES[TEST_BUCKET]],
    "dir": [*FILES_IN_BUCKETS[TEST_BUCKET], *LOTS_OF_FILES[TEST_BUCKET]],
    "dir_": [*FILES_IN_BUCKETS[TEST_BUCKET], *DIR_TWO_FILES],
    "dirty_data": DIRTY_DATA,
}


# NOTE: These tests currently compose multiple fixtures explicitly for readability.
@pytest.mark.s3
@pytest.mark.parametrize("dir_path", EXPECTED_FILE_LIST.keys())
def test_list_objects_returns_more_than_1000_entries(
    mock_s3_client: S3Client,
    dir_path: str,
) -> None:
    """Verify that pagination is followed so that more than 1000 objects are returned."""
    populate_mock_s3(mock_s3_client, FILES_IN_BUCKETS)
    # this adds two extra dirs to TEST_BUCKET with 1005 files in each
    populate_mock_s3(mock_s3_client, LOTS_OF_FILES)

    contents = list_objects(f"{TEST_BUCKET}/{dir_path}")
    keys = {obj["Key"] for obj in contents}
    assert keys == set(EXPECTED_FILE_LIST[dir_path])


@pytest.mark.s3
def test_list_objects_pass_respects_custom_max_keys(mock_s3_client: S3Client) -> None:
    """A small custom max_keys value still results in complete results via pagination."""
    populate_mock_s3(mock_s3_client, {TEST_BUCKET: FILES_IN_BUCKETS[TEST_BUCKET]})
    contents = list_objects(f"{TEST_BUCKET}/dir_one", max_keys=1)
    keys = {obj["Key"] for obj in contents}
    assert keys == {f for f in FILES_IN_BUCKETS[TEST_BUCKET] if f.startswith("dir_one")}


# upload_file
@pytest.mark.s3
@pytest.mark.parametrize("destination_dir", ["uploads", "uploads/", "some/uploads"])
@pytest.mark.parametrize("bucket", BUCKETS)
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
@pytest.mark.parametrize("file_name", [None, "", "custom.txt"])
def test_upload_file_succeeds(
    mock_s3_client: S3Client, sample_file: Path, protocol: str, bucket: str, destination_dir: str, file_name: str | None
) -> None:
    """Verify that a file is uploaded to the correct key in the specified bucket."""
    upload_args = {
        "local_file_path": sample_file,
        "destination_dir": f"{protocol}{bucket}/{destination_dir}",
    }
    if file_name is not None:
        upload_args["object_name"] = file_name

    result = upload_file(**upload_args)
    assert result is True
    expected_key = f"{destination_dir.removesuffix('/')}/{file_name or sample_file.name}"
    obj = mock_s3_client.get_object(Bucket=bucket, Key=expected_key)
    assert obj["Body"].read() == b"hello s3"


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize("path_type", [str, Path])
def test_upload_file_accepts_str_and_path(sample_file: Path, path_type: type[str] | type[Path]) -> None:
    """Verify that upload_file accepts both str and Path objects for the local file path."""
    result = upload_file(path_type(sample_file), f"{TEST_BUCKET}/uploads")
    assert result is True


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize(
    ("file_name", "error"),
    [
        (None, "No local file path specified"),
        ("", "No local file path specified"),
        # current working dir -- not a file
        (".", "No object_name"),
        # all the following are equivalent
        (Path(), "No object_name"),
        (Path(""), "No object_name"),  # noqa: PTH201
        (Path("."), "No object_name"),  # noqa: PTH201
    ],
)
def test_upload_file_fail_empty_object_name(file_name: str | Path | None, error: str) -> None:
    """A local_file_path whose .name is empty (e.g. '.') and no explicit object_name raises ValueError."""
    with pytest.raises(ValueError, match=error):
        upload_file(file_name, f"{TEST_BUCKET}/uploads")  # pyright: ignore[reportArgumentType]


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_file_error(sample_file: Path) -> None:
    """Verify that upload_file raises ValueError when no destination directory is provided."""
    with pytest.raises(ValueError, match="No destination directory"):
        upload_file(sample_file, "")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_file_fail_upload_error_returns_false(sample_file: Path, caplog: pytest.LogCaptureFixture) -> None:
    """A failure during the actual upload (e.g. destination bucket does not exist) returns False, not an exception."""
    with caplog.at_level(logging.ERROR):
        result = upload_file(sample_file, "nonexistent-bucket/uploads")
    assert result is False
    assert any("Error uploading to s3" in r.message for r in caplog.records)


@pytest.mark.s3
def test_upload_file_with_metadata_still_applies_default_extra_args(
    mock_s3_client: S3Client, sample_file: Path
) -> None:
    """DEFAULT_EXTRA_ARGS (e.g. ChecksumAlgorithm) is preserved even when user_metadata is also supplied."""
    upload_file(sample_file, f"{TEST_BUCKET}/uploads", user_metadata={"k": "v"})
    resp = mock_s3_client.head_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}", ChecksumMode="ENABLED")
    assert "ChecksumCRC64NVME" in resp
    assert resp["Metadata"]["k"] == "v"


# transfer_config_kwargs
@pytest.mark.s3
def test_upload_file_pass_transfer_config_kwargs_forwarded(
    mock_s3_client: S3Client, sample_file: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify that transfer_config_kwargs passed to upload_file reach build_transfer_config."""
    captured: dict[str, Any] = {}

    def fake_build_transfer_config(file_size: int, **kwargs: Any) -> TransferConfig:
        captured["file_size"] = file_size
        captured["kwargs"] = kwargs
        return TransferConfig()

    monkeypatch.setattr(object_utils, "build_transfer_config", fake_build_transfer_config)

    result = upload_file(
        sample_file,
        f"{TEST_BUCKET}/uploads",
        transfer_config_kwargs={"max_concurrency": 5, "use_threads": False},
    )

    assert result is True
    assert captured["kwargs"] == {"max_concurrency": 5, "use_threads": False}
    assert captured["file_size"] == sample_file.stat().st_size

    obj = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}")
    assert obj["Body"].read() == b"hello s3"


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_file_pass_no_transfer_config_kwargs_uses_empty_dict(
    sample_file: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """When transfer_config_kwargs is not supplied, build_transfer_config still receives no extra kwargs."""
    captured: dict[str, Any] = {}

    def fake_build_transfer_config(_: int, **kwargs: Any) -> TransferConfig:
        captured["kwargs"] = kwargs
        return TransferConfig()

    monkeypatch.setattr(object_utils, "build_transfer_config", fake_build_transfer_config)

    result = upload_file(sample_file, f"{TEST_BUCKET}/uploads")

    assert result is True
    assert captured["kwargs"] == {}


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_file_pass_show_progress_false_disables_progress_bar(sample_file: Path) -> None:
    """Verify that show_progress=False disables the tqdm progress bar for uploads."""
    captured_disable: list[bool] = []
    original_make_progress_bar = object_utils.make_progress_bar

    def spy_make_progress_bar(*, disable: bool, **kwargs: Any) -> Any:
        captured_disable.append(disable)
        return original_make_progress_bar(disable=disable, **kwargs)

    with patch.object(object_utils, "make_progress_bar", side_effect=spy_make_progress_bar):
        result = upload_file(sample_file, f"{TEST_BUCKET}/uploads", show_progress=False)

    assert result is True
    assert captured_disable == [True]


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_file_pass_show_progress_true_enables_progress_bar(sample_file: Path) -> None:
    """Verify that show_progress=True (the default) enables the tqdm progress bar for uploads."""
    captured_disable: list[bool] = []
    original_make_progress_bar = object_utils.make_progress_bar

    def spy_make_progress_bar(*, disable: bool, **kwargs: Any) -> Any:
        captured_disable.append(disable)
        return original_make_progress_bar(disable=disable, **kwargs)

    with patch.object(object_utils, "make_progress_bar", side_effect=spy_make_progress_bar):
        result = upload_file(sample_file, f"{TEST_BUCKET}/uploads", show_progress=True)

    assert result is True
    assert captured_disable == [False]


@pytest.mark.s3
def test_upload_file_overwrites_existing_object_by_default(mock_s3_client: S3Client, sample_file: Path) -> None:
    """Since skip_if_exists is not yet implemented, upload_file always overwrites an existing object."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}", Body=b"old")
    result = upload_file(sample_file, f"{TEST_BUCKET}/uploads")
    assert result is True
    obj = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}")
    assert obj["Body"].read() == b"hello s3"


@pytest.mark.s3
@pytest.mark.xfail(reason="skip_if_exists is not yet implemented (TODO in upload_fileobj)", strict=True)
def test_upload_file_skip_if_exists_skips_when_already_present(mock_s3_client: S3Client, sample_file: Path) -> None:
    """Once implemented, skip_if_exists=True should skip the upload, leaving the existing object untouched.

    Marked xfail(strict=True) so that this test starts *failing* (unexpected pass) the moment
    skip_if_exists is actually implemented — a forcing function to come remove this marker.
    """
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}", Body=b"old")
    result = upload_file(sample_file, f"{TEST_BUCKET}/uploads", skip_if_exists=True)
    assert result is True
    obj = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}")
    assert obj["Body"].read() == b"old"


# upload_file => upload_fileobj
@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_file_forwards_file_path_and_size_to_upload_fileobj(
    sample_file: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify upload_file passes file_path and file_size through to upload_fileobj correctly."""
    captured: dict[str, Any] = {}
    original_upload_fileobj = object_utils.upload_fileobj

    def spy_upload_fileobj(fileobj: Any, s3_path: str, **kwargs: Any) -> bool:
        captured["s3_path"] = s3_path
        captured["file_path"] = kwargs.get("file_path")
        captured["file_size"] = kwargs.get("file_size")
        return original_upload_fileobj(fileobj, s3_path, **kwargs)

    monkeypatch.setattr(object_utils, "upload_fileobj", spy_upload_fileobj)
    result = upload_file(sample_file, f"{TEST_BUCKET}/uploads")

    assert result is True
    assert captured["s3_path"] == f"{TEST_BUCKET}/uploads/{sample_file.name}"
    assert captured["file_path"] == sample_file
    assert captured["file_size"] == sample_file.stat().st_size


# upload_fileobj (direct)
@pytest.mark.s3
@pytest.mark.parametrize("destination_dir", ["uploads", "uploads/", "some/uploads"])
@pytest.mark.parametrize("bucket", BUCKETS)
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_upload_fileobj_succeeds(mock_s3_client: S3Client, protocol: str, bucket: str, destination_dir: str) -> None:
    """Verify upload_fileobj uploads stream content to the correct key, with any valid protocol prefix."""
    data = b"streamed bytes content"
    key = f"{destination_dir.removesuffix('/')}/stream.bin"
    s3_path = f"{protocol}{bucket}/{key}"
    result = upload_fileobj(io.BytesIO(data), s3_path)
    assert result is True
    obj = mock_s3_client.get_object(Bucket=bucket, Key=key)
    assert obj["Body"].read() == data


@pytest.mark.s3
def test_upload_fileobj_default_file_size_is_zero_when_unspecified(
    mock_s3_client: S3Client, monkeypatch: pytest.MonkeyPatch
) -> None:
    """When file_size is not supplied, build_transfer_config receives 0 (the declared default)."""
    captured: dict[str, Any] = {}

    def fake_build_transfer_config(file_size: int, **kwargs: Any) -> TransferConfig:
        captured["file_size"] = file_size
        return TransferConfig()

    monkeypatch.setattr(object_utils, "build_transfer_config", fake_build_transfer_config)
    result = upload_fileobj(io.BytesIO(b"data"), f"{TEST_BUCKET}/uploads/stream.bin")
    assert result is True
    assert captured["file_size"] == 0


@pytest.mark.s3
def test_upload_fileobj_explicit_file_size_used_for_progress_bar(mock_s3_client: S3Client) -> None:
    """Passing file_size explicitly sizes the progress bar accordingly."""
    data = b"12345"
    captured_totals: list[int | None] = []
    original_make_progress_bar = object_utils.make_progress_bar

    def spy_make_progress_bar(*, total: int | None, **kwargs: Any) -> Any:
        captured_totals.append(total)
        return original_make_progress_bar(total=total, **kwargs)

    with patch.object(object_utils, "make_progress_bar", side_effect=spy_make_progress_bar):
        upload_fileobj(io.BytesIO(data), f"{TEST_BUCKET}/uploads/sized.bin", file_size=len(data))

    assert captured_totals == [len(data)]


@pytest.mark.s3
def test_upload_fileobj_applies_default_extra_args_and_metadata(mock_s3_client: S3Client) -> None:
    """DEFAULT_EXTRA_ARGS and user_metadata are both applied when calling upload_fileobj directly."""
    upload_fileobj(io.BytesIO(b"payload"), f"{TEST_BUCKET}/uploads/meta.bin", user_metadata={"k": "v"})
    resp = mock_s3_client.head_object(Bucket=TEST_BUCKET, Key="uploads/meta.bin", ChecksumMode="ENABLED")
    assert "ChecksumCRC64NVME" in resp
    assert resp["Metadata"]["k"] == "v"


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_fileobj_fail_upload_error_returns_false(caplog: pytest.LogCaptureFixture) -> None:
    """A failure during the actual upload (e.g. destination bucket does not exist) returns False, not an exception."""
    with caplog.at_level(logging.ERROR):
        result = upload_fileobj(io.BytesIO(b"data"), "nonexistent-bucket/uploads/file.bin")
    assert result is False
    assert any("Error uploading to s3" in r.message for r in caplog.records)


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_fileobj_raises_for_invalid_s3_path() -> None:
    """An s3_path that can't be split into bucket/key propagates the ValueError from split_s3_path."""
    with pytest.raises(ValueError, match="Could not parse"):
        upload_fileobj(io.BytesIO(b"data"), "no-slash-path")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_fileobj_pass_transfer_config_kwargs_forwarded(monkeypatch: pytest.MonkeyPatch) -> None:
    """Verify transfer_config_kwargs passed directly to upload_fileobj reach build_transfer_config."""
    captured: dict[str, Any] = {}

    def fake_build_transfer_config(file_size: int, **kwargs: Any) -> TransferConfig:
        captured["file_size"] = file_size
        captured["kwargs"] = kwargs
        return TransferConfig()

    monkeypatch.setattr(object_utils, "build_transfer_config", fake_build_transfer_config)
    data = b"transfer config test"
    result = upload_fileobj(
        io.BytesIO(data),
        f"{TEST_BUCKET}/uploads/tc.bin",
        file_size=len(data),
        transfer_config_kwargs={"max_concurrency": 2},
    )
    assert result is True
    assert captured["kwargs"] == {"max_concurrency": 2}
    assert captured["file_size"] == len(data)


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_fileobj_pass_show_progress_false_disables_progress_bar() -> None:
    """Verify that show_progress=False disables the tqdm progress bar when calling upload_fileobj directly."""
    captured_disable: list[bool] = []
    original_make_progress_bar = object_utils.make_progress_bar

    def spy_make_progress_bar(*, disable: bool, **kwargs: Any) -> Any:
        captured_disable.append(disable)
        return original_make_progress_bar(disable=disable, **kwargs)

    with patch.object(object_utils, "make_progress_bar", side_effect=spy_make_progress_bar):
        result = upload_fileobj(io.BytesIO(b"data"), f"{TEST_BUCKET}/uploads/np.bin", show_progress=False)

    assert result is True
    assert captured_disable == [True]


@pytest.mark.s3
def test_upload_fileobj_skip_if_exists_flag_currently_logs_placeholder_only(
    mock_s3_client: S3Client, caplog: pytest.LogCaptureFixture
) -> None:
    """Documents current stub behavior: skip_if_exists=True logs a TODO but does not skip the upload."""
    with caplog.at_level(logging.DEBUG):
        result = upload_fileobj(io.BytesIO(b"data"), f"{TEST_BUCKET}/uploads/skip.bin", skip_if_exists=True)
    assert result is True
    assert any("To be implemented" in r.message for r in caplog.records)
    obj = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key="uploads/skip.bin")
    assert obj["Body"].read() == b"data"


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_fileobj_logs_file_path_when_provided(sample_file: Path, caplog: pytest.LogCaptureFixture) -> None:
    """The debug log includes the local file path when file_path is supplied."""
    with caplog.at_level(logging.DEBUG), sample_file.open("rb") as fh:
        upload_fileobj(fh, f"{TEST_BUCKET}/uploads/{sample_file.name}", file_path=sample_file)
    assert any(f"uploading {sample_file} to" in r.message for r in caplog.records)


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_fileobj_logs_generic_message_without_file_path(caplog: pytest.LogCaptureFixture) -> None:
    """The debug log uses the generic 'fileobj' message when file_path is not supplied."""
    with caplog.at_level(logging.DEBUG):
        upload_fileobj(io.BytesIO(b"data"), f"{TEST_BUCKET}/uploads/stream.bin")
    assert any("uploading fileobj to" in r.message for r in caplog.records)


# upload_dir
@pytest.mark.s3
@pytest.mark.parametrize("bucket", [TEST_BUCKET, ALT_BUCKET])
def test_upload_dir_uploads_recursively(mock_s3_client: S3Client, bucket: str, sample_dir: Path) -> None:
    """Verify that upload_dir recurses into subdirectories and uploads nested files."""
    result = upload_dir(sample_dir, f"{bucket}/remote")
    assert result is True
    keys = {obj["Key"] for obj in mock_s3_client.list_objects_v2(Bucket=bucket)["Contents"]}
    assert keys == {f"remote/{f}" for f in SAMPLE_FILES}


@pytest.mark.s3
@pytest.mark.parametrize("path_type", [str, Path])
def test_upload_dir_accepts_str_and_path(
    mock_s3_client: S3Client, sample_dir: Path, path_type: type[str] | type[Path]
) -> None:
    """Verify that upload_dir accepts both str and Path objects for the local directory path."""
    bucket = TEST_BUCKET
    result = upload_dir(path_type(sample_dir), f"{bucket}/remote")
    assert result is True
    keys = {obj["Key"] for obj in mock_s3_client.list_objects_v2(Bucket=bucket)["Contents"]}
    assert keys == {f"remote/{f}" for f in SAMPLE_FILES}


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_dir_raises_on_empty_source() -> None:
    """Verify that upload_dir raises ValueError when no source directory is provided."""
    with pytest.raises(ValueError, match="No source directory"):
        upload_dir("", f"{TEST_BUCKET}/remote")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_dir_raises_on_empty_destination(sample_dir: Path) -> None:
    """Verify that upload_dir raises ValueError when no destination directory is provided."""
    with pytest.raises(ValueError, match="No destination directory"):
        upload_dir(sample_dir, "")


@pytest.mark.s3
def test_upload_dir_pass_respects_custom_file_glob(mock_s3_client: S3Client, sample_dir: Path) -> None:
    """A custom file_glob restricts which files get uploaded."""
    result = upload_dir(sample_dir, f"{TEST_BUCKET}/remote", file_glob="dir_one/*.txt")
    assert result is True
    keys = {obj["Key"] for obj in mock_s3_client.list_objects_v2(Bucket=TEST_BUCKET)["Contents"]}
    assert keys == {"remote/dir_one/file1.txt", "remote/dir_one/file2.txt"}


@pytest.mark.s3
def test_upload_dir_pass_empty_directory_returns_true_and_uploads_nothing(
    mock_s3_client: S3Client, tmp_path: Path
) -> None:
    """An empty source directory (no matching files) returns True without uploading anything."""
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    result = upload_dir(empty_dir, f"{TEST_BUCKET}/remote")
    assert result is True
    assert "Contents" not in mock_s3_client.list_objects_v2(Bucket=TEST_BUCKET)


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_dir_fail_partial_failure_returns_false(sample_dir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """If any individual file upload fails, upload_dir returns False even though other files succeeded."""
    original_upload_file = object_utils.upload_file
    call_count = {"n": 0}

    def flaky_upload_file(*args: Any, **kwargs: Any) -> bool:
        call_count["n"] += 1
        if call_count["n"] == 1:
            return False
        return original_upload_file(*args, **kwargs)

    monkeypatch.setattr(object_utils, "upload_file", flaky_upload_file)
    result = object_utils.upload_dir(sample_dir, f"{TEST_BUCKET}/remote")
    assert result is False
    assert call_count["n"] == len(SAMPLE_FILES)


# download_file
@pytest.mark.s3
@pytest.mark.parametrize("bucket", BUCKETS)
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_download_file_retrieves_correct_content(
    mock_s3_client: S3Client, protocol: str, bucket: str, tmp_path: Path
) -> None:
    """Verify that download_file writes the correct file content to disk for each valid bucket."""
    content = b"some important content"
    mock_s3_client.put_object(Bucket=bucket, Key="remote/data.txt", Body=content, **DEFAULT_EXTRA_ARGS)
    local_file = str(tmp_path / "data.txt")
    download_file(f"{protocol}{bucket}/remote/data.txt", local_file)
    assert Path(local_file).read_bytes() == content


@pytest.mark.s3
@pytest.mark.parametrize("path_type", [str, Path])
def test_download_file_use_str_or_path_for_local_file(
    mock_s3_client: S3Client, tmp_path: Path, path_type: type[str] | type[Path]
) -> None:
    """Verify that download_file can create a new directory if need be."""
    content = b"some cool file stuff"
    bucket = BUCKETS[0]
    key = "to/the/door.txt"
    mock_s3_client.put_object(Bucket=bucket, Key=key, Body=content, **DEFAULT_EXTRA_ARGS)
    assert object_exists(f"{bucket}/{key}")
    local_file = tmp_path / "file.txt"
    assert not local_file.exists()
    download_file(f"{bucket}/{key}", path_type(local_file))
    assert local_file.read_bytes() == content


@pytest.mark.s3
def test_download_file_save_to_new_dir(mock_s3_client: S3Client, tmp_path: Path) -> None:
    """Verify that download_file can create a new directory if need be."""
    content = b"some cool file stuff"
    bucket = BUCKETS[0]
    key = "to/the/door.txt"
    mock_s3_client.put_object(Bucket=bucket, Key=key, Body=content, **DEFAULT_EXTRA_ARGS)
    assert object_exists(f"{bucket}/{key}")
    local_file = tmp_path / "some" / "convoluted" / "path" / "to" / "file.txt"
    assert not local_file.exists()
    assert not local_file.parents[2].exists()
    download_file(f"{bucket}/{key}", local_file)
    assert local_file.read_bytes() == content


@pytest.mark.s3
def test_download_file_clobbers_existing_file(mock_s3_client: S3Client, tmp_path: Path) -> None:
    """Verify that download_file can create a new directory if need be."""
    bucket = BUCKETS[0]
    key = "to/the/door.txt"
    local_file_content = b"some old crap"
    remote_file_content = b"some remote crap"
    local_file = tmp_path / "file.txt"
    local_file.write_bytes(local_file_content)

    mock_s3_client.put_object(Bucket=bucket, Key=key, Body=remote_file_content, **DEFAULT_EXTRA_ARGS)
    assert object_exists(f"{bucket}/{key}")

    assert local_file.exists()
    assert local_file.read_bytes() == local_file_content

    download_file(f"{bucket}/{key}", local_file)
    assert local_file.read_bytes() == remote_file_content


@pytest.mark.s3
def test_download_file_does_not_clobber_existing_file_to_mkdir(mock_s3_client: S3Client, tmp_path: Path) -> None:
    """Verify that download_file will not overwrite an existing file whilst trying to make a directory."""
    bucket = BUCKETS[0]
    key = "to/the/door.txt"
    mock_s3_client.put_object(Bucket=bucket, Key=key, Body=b"some crappy nonsense", **DEFAULT_EXTRA_ARGS)
    assert object_exists(f"{bucket}/{key}")
    local_file = tmp_path / "to"
    local_file.touch()

    with pytest.raises(FileExistsError, match=f"File exists: '{local_file!s}'"):
        download_file(f"{bucket}/{key}", local_file / "file.txt")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_download_file_does_not_exist(tmp_path: Path) -> None:
    """Ensure that attempting to download a file that does not exist raises an error."""
    bucket = BUCKETS[0]
    key = "to/the/door.txt"
    assert not object_exists(f"{bucket}/{key}")

    with pytest.raises(
        FileNotFoundError,
        match=f"File not found: {bucket}/{key}",
    ):
        download_file(f"{bucket}/{key}", tmp_path / "file.txt")


@pytest.mark.s3
def test_download_file_fail_non_404_error_propagates(
    mock_s3_client: S3Client, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A non-404 ClientError while checking object existence propagates rather than being treated as 'not found'."""

    def raise_access_denied(**kwargs: Any) -> None:
        raise ClientError(
            {"Error": {"Code": "AccessDenied"}, "ResponseMetadata": {"HTTPStatusCode": 403}},
            "HeadObject",
        )

    monkeypatch.setattr(mock_s3_client, "head_object", raise_access_denied)
    with pytest.raises(ClientError, match="AccessDenied"):
        download_file(f"{TEST_BUCKET}/some/key.txt", tmp_path / "file.txt")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_download_file_fail_directory_creation_error_other_than_exists_propagates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    """A directory-creation failure other than 'already exists' (e.g. permission denied) is logged and re-raised."""
    err_msg = "Permission denied"

    def raise_permission_error(self: Path, *args: Any, **kwargs: Any) -> None:
        raise PermissionError(err_msg)

    monkeypatch.setattr(Path, "mkdir", raise_permission_error)

    with pytest.raises(PermissionError, match=err_msg):
        download_file(f"{TEST_BUCKET}/some/key.txt", tmp_path / "newdir" / "file.txt")

    assert any("Could not save s3 file" in r.message for r in caplog.records)


@pytest.mark.s3
def test_download_file_pass_transfer_config_kwargs_forwarded(
    mock_s3_client: S3Client, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Verify that transfer_config_kwargs passed to download_file reach build_transfer_config."""
    bucket = TEST_BUCKET
    key = "transfer/config.txt"
    content = b"config test content"
    mock_s3_client.put_object(Bucket=bucket, Key=key, Body=content, **DEFAULT_EXTRA_ARGS)

    captured: dict[str, Any] = {}

    def fake_build_transfer_config(file_size: int, **kwargs: Any) -> TransferConfig:
        captured["file_size"] = file_size
        captured["kwargs"] = kwargs
        return TransferConfig()

    monkeypatch.setattr(object_utils, "build_transfer_config", fake_build_transfer_config)

    local_file = tmp_path / "config.txt"
    download_file(
        f"{bucket}/{key}",
        local_file,
        transfer_config_kwargs={"max_concurrency": 3},
    )

    assert local_file.read_bytes() == content
    assert captured["kwargs"] == {"max_concurrency": 3}
    assert captured["file_size"] == len(content)


@pytest.mark.s3
def test_download_file_pass_no_transfer_config_kwargs_uses_empty_dict(
    mock_s3_client: S3Client, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """When transfer_config_kwargs is not supplied, build_transfer_config still receives no extra kwargs."""
    bucket = TEST_BUCKET
    key = "transfer/config2.txt"
    mock_s3_client.put_object(Bucket=bucket, Key=key, Body=b"data", **DEFAULT_EXTRA_ARGS)

    captured: dict[str, Any] = {}

    def fake_build_transfer_config(_: int, **kwargs: Any) -> TransferConfig:
        captured["kwargs"] = kwargs
        return TransferConfig()

    monkeypatch.setattr(object_utils, "build_transfer_config", fake_build_transfer_config)

    download_file(f"{bucket}/{key}", tmp_path / "config2.txt")

    assert captured["kwargs"] == {}


@pytest.mark.s3
def test_download_file_pass_show_progress_false_disables_progress_bar(mock_s3_client: S3Client, tmp_path: Path) -> None:
    """Verify that show_progress=False disables the tqdm progress bar for downloads."""
    bucket = TEST_BUCKET
    key = "progress/no_bar.txt"
    mock_s3_client.put_object(Bucket=bucket, Key=key, Body=b"content", **DEFAULT_EXTRA_ARGS)

    captured_disable: list[bool] = []
    original_make_progress_bar = object_utils.make_progress_bar

    def spy_make_progress_bar(*, disable: bool, **kwargs: Any) -> Any:
        captured_disable.append(disable)
        return original_make_progress_bar(disable=disable, **kwargs)

    with patch.object(object_utils, "make_progress_bar", side_effect=spy_make_progress_bar):
        download_file(f"{bucket}/{key}", tmp_path / "no_bar.txt", show_progress=False)

    assert captured_disable == [True]


# copy_object
@pytest.mark.s3
@pytest.mark.parametrize("destination", BUCKETS)
def test_copy_object(mock_s3_client: S3Client, destination: str) -> None:
    """Verify that copy_object copies an object to a new key within the same bucket."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key="src/file.txt", Body=b"copy me")
    assert object_exists(f"{TEST_BUCKET}/src/file.txt")
    response = copy_object(f"{TEST_BUCKET}/src/file.txt", f"{destination}/dst/path/to/file.txt")

    # check both objects exist
    assert object_exists(f"{TEST_BUCKET}/src/file.txt")
    assert object_exists(f"{destination}/dst/path/to/file.txt")

    obj = mock_s3_client.get_object(Bucket=destination, Key="dst/path/to/file.txt")
    assert obj["Body"].read() == b"copy me"
    assert response["ResponseMetadata"]["HTTPStatusCode"] == HTTP_200


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_copy_object_source_object_nonexistent() -> None:
    """Ensure that the code throws an error if the source object does not exist."""
    s3_path = f"{TEST_BUCKET}/some/path/to/file"
    assert object_exists(s3_path) is False
    with pytest.raises(Exception, match="The specified key does not exist"):
        copy_object(s3_path, f"{TEST_BUCKET}/a/different/path/to/file")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_copy_object_source_bucket_nonexistent() -> None:
    """Ensure that the code throws an error if the bucket does not exist."""
    s3_path = "some-bucket/some/path/to/file"
    assert object_exists(s3_path) is False
    with pytest.raises(Exception, match="The specified bucket does not exist"):
        copy_object(s3_path, f"{TEST_BUCKET}/a/different/path/to/file")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_copy_file_source_object_nonexistent() -> None:
    """Ensure that the code throws an error if the source object does not exist."""
    s3_path = f"{TEST_BUCKET}/some/path/to/file"
    assert object_exists(s3_path) is False
    with pytest.raises(Exception, match="The specified key does not exist"):
        copy_object(s3_path, f"{TEST_BUCKET}/a/different/path/to/file")


# copy_directory tests
def put_objects(mock_s3_client: S3Client, bucket: str, keys: list[str], body: bytes = b"data") -> None:
    """Helper to seed objects into a bucket."""
    for key in keys:
        mock_s3_client.put_object(Bucket=bucket, Key=key, Body=body)


def list_keys(mock_s3_client: S3Client, bucket: str, prefix: str = "") -> set[str]:
    """Helper to list all keys in a bucket under a prefix."""
    paginator = mock_s3_client.get_paginator("list_objects_v2")
    keys = []
    for page in paginator.paginate(Bucket=bucket, Prefix=prefix):
        keys.extend(obj["Key"] for obj in page.get("Contents", []))
    return set(keys)


@pytest.mark.s3
@pytest.mark.parametrize("source_suffix", ["", "/"])
@pytest.mark.parametrize("dest_suffix", ["", "/"])
def test_copy_directory_copies_all_objects_to_dest(
    mock_s3_client: S3Client, source_suffix: str, dest_suffix: str
) -> None:
    """Verify that all objects under the source prefix are present in the successes dict.

    Ensure that copy works correctly with or without a slash at the end of the directory name.
    """
    populate_mock_s3(mock_s3_client, {TEST_BUCKET: FILES_IN_BUCKETS[TEST_BUCKET]})
    source_bucket_files = list_keys(mock_s3_client, TEST_BUCKET)
    assert set(source_bucket_files) == set(SAMPLE_FILES)
    dest_bucket_files = list_keys(mock_s3_client, ALT_BUCKET)
    assert dest_bucket_files == set()

    successes, errors = copy_directory(
        f"s3://{TEST_BUCKET}/dir_one{source_suffix}", f"s3://{ALT_BUCKET}/some/destination/dir{dest_suffix}"
    )

    assert errors == {}
    expected_files = {
        f"{TEST_BUCKET}/{f}": f"{ALT_BUCKET}/some/destination/dir{f.replace('dir_one', '')}" for f in SAMPLE_FILES
    }
    assert successes == expected_files
    # ensure that the original files are still in place
    assert set(list_keys(mock_s3_client, TEST_BUCKET)) == set(SAMPLE_FILES)
    # destination should have new files from the source
    assert set(list_keys(mock_s3_client, ALT_BUCKET)) == {
        f.removeprefix(f"{ALT_BUCKET}/") for f in expected_files.values()
    }

    # check the content
    for src, dest in expected_files.items():
        src_resp = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key=src.removeprefix(f"{TEST_BUCKET}/"))
        assert src_resp["Body"].read() == src.encode()
        dest_resp = mock_s3_client.get_object(Bucket=ALT_BUCKET, Key=dest.removeprefix(f"{ALT_BUCKET}/"))
        assert dest_resp["Body"].read() == src.encode()


@pytest.mark.s3
def test_copy_directory_copy_within_same_bucket(mock_s3_client: S3Client) -> None:
    """Verify that copying between two prefixes within the same bucket works correctly."""
    populate_mock_s3(mock_s3_client, {TEST_BUCKET: ["foo/a.txt", "foo/b.txt"]})

    successes, errors = copy_directory(f"s3://{TEST_BUCKET}/foo", f"s3://{TEST_BUCKET}/bar")

    assert successes == {
        f"{TEST_BUCKET}/foo/a.txt": f"{TEST_BUCKET}/bar/a.txt",
        f"{TEST_BUCKET}/foo/b.txt": f"{TEST_BUCKET}/bar/b.txt",
    }
    assert errors == {}
    assert list_keys(mock_s3_client, TEST_BUCKET, prefix="bar") == {"bar/a.txt", "bar/b.txt"}
    assert list_keys(mock_s3_client, TEST_BUCKET) == {"foo/a.txt", "foo/b.txt", "bar/a.txt", "bar/b.txt"}


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_copy_directory_empty_directory_returns_empty_dicts() -> None:
    """Verify that when the source prefix matches no objects, both the successes and errors dictionaries are returned empty."""
    successes, errors = copy_directory(f"s3://{TEST_BUCKET}/nonexistent/", f"s3://{ALT_BUCKET}/bar/")

    assert successes == {}
    assert errors == {}


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_copy_directory_does_not_copy_objects_outside_prefix(mock_s3_client: S3Client) -> None:
    """Verify that objects whose keys share a prefix string but are not under the source directory.

    Example: 'foobar/' when copying 'foo/'.
    """
    populate_mock_s3(
        mock_s3_client,
        {
            TEST_BUCKET: [
                "foo/a.txt",
                "foobar/should-not-be-copied.txt",
            ]
        },
    )
    copy_directory(f"s3://{TEST_BUCKET}/foo", f"s3://{ALT_BUCKET}/bar")
    assert list_keys(mock_s3_client, ALT_BUCKET) == {"bar/a.txt"}


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_copy_directory_missing_source_bucket_returns_error() -> None:
    """Verify that when the source bucket does not exist, botocore throws an error."""
    # FIXME: throws a s3.Client.exceptions.NoSuchBucket
    with pytest.raises(Exception, match="The specified bucket does not exist"):
        copy_directory("s3://nonexistent-bucket/bar", f"s3://{TEST_BUCKET}/foo")


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_copy_directory_missing_dest_bucket_records_errors(mock_s3_client: S3Client) -> None:
    """Verify that when the destination bucket does not exist, the errors dict contains all objects under the original dir."""
    # FIXME: throw a bucket not exists error?
    populate_mock_s3(mock_s3_client, {TEST_BUCKET: ["foo/a.txt", "foo/b.txt"]})

    # with pytest.raises(FileNotFoundError, match="The specified bucket does not exist"):
    successes, errors = copy_directory(f"s3://{TEST_BUCKET}/foo", "s3://nonexistent-bucket/bar")

    assert successes == {}
    assert f"{TEST_BUCKET}/foo/a.txt" in errors
    assert f"{TEST_BUCKET}/foo/b.txt" in errors
    assert isinstance(errors[f"{TEST_BUCKET}/foo/a.txt"], Exception)


@pytest.mark.s3
def test_copy_directory_records_error_for_non_200_response(
    mock_s3_client: S3Client, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A copy_object call that returns (rather than raises) a non-200 status is recorded as an error."""
    populate_mock_s3(mock_s3_client, {TEST_BUCKET: ["foo/a.txt"]})

    fake_response = {"ResponseMetadata": {"HTTPStatusCode": 500}}
    monkeypatch.setattr(mock_s3_client, "copy_object", lambda **kwargs: fake_response)  # noqa: ARG005

    successes, errors = copy_directory(f"s3://{TEST_BUCKET}/foo", f"s3://{ALT_BUCKET}/bar")

    assert successes == {}
    assert errors == {f"{TEST_BUCKET}/foo/a.txt": fake_response}


# upload_file with metadata
@pytest.mark.s3
@pytest.mark.parametrize("bucket", BUCKETS)
def test_upload_file_with_metadata_attaches_metadata(mock_s3_client: S3Client, sample_file: Path, bucket: str) -> None:
    """Verify that upload_file with metadata stores user metadata on the uploaded object."""
    metadata = {"md5": "abc123", "source": "ncbi"}
    result = upload_file(sample_file, f"{bucket}/uploads", user_metadata=metadata)
    assert result is True

    resp = mock_s3_client.head_object(Bucket=bucket, Key=f"uploads/{sample_file.name}")
    assert resp["Metadata"]["md5"] == "abc123"
    assert resp["Metadata"]["source"] == "ncbi"


@pytest.mark.s3
def test_upload_file_with_metadata_custom_object_name(mock_s3_client: S3Client, sample_file: Path) -> None:
    """Verify that the object_name parameter overrides the filename."""
    result = upload_file(sample_file, f"{TEST_BUCKET}/uploads", user_metadata={"k": "v"}, object_name="renamed.txt")
    assert result is True
    obj = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key="uploads/renamed.txt")
    assert obj["Body"].read() == b"hello s3"


@pytest.mark.s3
def test_upload_file_with_metadata_overwrites_existing(mock_s3_client: S3Client, sample_file: Path) -> None:
    """Verify that upload_file with metadata uploads even when the object already exists."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}", Body=b"old")
    result = upload_file(sample_file, f"{TEST_BUCKET}/uploads", user_metadata={"new": "true"})
    assert result is True
    obj = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key=f"uploads/{sample_file.name}")
    assert obj["Body"].read() == b"hello s3"


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_upload_file_with_metadata_raises_on_empty_destination(sample_file: Path) -> None:
    """Verify ValueError when destination_dir is empty."""
    with pytest.raises(ValueError, match="No destination directory"):
        upload_file(sample_file, "", user_metadata={"k": "v"})


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
@pytest.mark.parametrize("path_type", [str, Path])
def test_upload_file_with_metadata_accepts_str_and_path(sample_file: Path, path_type: type[str] | type[Path]) -> None:
    """Verify that upload_file with metadata accepts both str and Path."""
    result = upload_file(path_type(sample_file), f"{TEST_BUCKET}/uploads", user_metadata={})
    assert result is True


# copy_object
@pytest.mark.s3
@pytest.mark.parametrize("destination", BUCKETS)
def test_copy_object_preserves_user_metadata(mock_s3_client: S3Client, destination: str) -> None:
    """copy_object preserves source user metadata (MetadataDirective=COPY default)."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key="src/file.txt", Body=b"archive me", Metadata={"md5": "abc123"})
    response = copy_object(
        f"{TEST_BUCKET}/src/file.txt",
        f"{destination}/archive/file.txt",
    )
    assert response["ResponseMetadata"]["HTTPStatusCode"] == HTTP_200

    # source user metadata is preserved (MetadataDirective=COPY)
    resp = mock_s3_client.head_object(Bucket=destination, Key="archive/file.txt")
    assert resp["Metadata"].get("md5") == "abc123"

    # verify source still exists
    assert object_exists(f"{TEST_BUCKET}/src/file.txt")


@pytest.mark.s3
def test_copy_object_preserves_content(mock_s3_client: S3Client) -> None:
    """Verify that the content of the copied object matches the original."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key="src/data.bin", Body=b"binary data")
    copy_object(
        f"{TEST_BUCKET}/src/data.bin",
        f"{TEST_BUCKET}/dst/data.bin",
    )
    obj = mock_s3_client.get_object(Bucket=TEST_BUCKET, Key="dst/data.bin")
    assert obj["Body"].read() == b"binary data"


# delete_object
@pytest.mark.s3
@pytest.mark.parametrize("bucket", BUCKETS)
@pytest.mark.parametrize("protocol", ["", "s3://", "s3a://"])
def test_delete_object_removes_object(mock_s3_client: S3Client, bucket: str, protocol: str) -> None:
    """Verify that delete_object removes the object from the specified bucket."""
    mock_s3_client.put_object(Bucket=bucket, Key="to/delete.txt", Body=b"bye")
    s3_path = f"{protocol}{bucket}/to/delete.txt"
    assert object_exists(s3_path) is True

    resp = delete_object(s3_path)
    assert object_exists(s3_path) is False
    assert resp.get("ResponseMetadata", {}).get("HTTPStatusCode") == HTTP_204

    # retry the deletion
    resp = delete_object(s3_path)
    assert object_exists(s3_path) is False
    assert resp.get("ResponseMetadata", {}).get("HTTPStatusCode") == HTTP_204


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_delete_object_pass_non_existent_file() -> None:
    """Verify that delete_object does not raise an error when the file does not exist."""
    test_file = f"{TEST_BUCKET}/nonexistent-file.txt"
    assert object_exists(test_file) is False
    resp = delete_object(test_file)
    assert resp.get("ResponseMetadata", {}).get("HTTPStatusCode") == HTTP_204


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_delete_object_fail_no_such_bucket() -> None:
    """Verify that delete_object removes the object from the specified bucket."""
    test_file = "made-up-bucket/nonexistent-file.txt"
    assert object_exists(test_file) is False
    with pytest.raises(ClientError, match="The specified bucket does not exist"):
        delete_object(test_file)


# delete_objects
@pytest.mark.s3
def test_delete_objects_pass_splits_into_batches_of_1000(mock_s3_client: S3Client) -> None:
    """More than 1000 keys are deleted via multiple batched API calls, each capped at 1000 objects."""
    total_keys = 1005
    keys = [f"batch/file_{i:04d}.txt" for i in range(total_keys)]
    for k in keys:
        mock_s3_client.put_object(Bucket=TEST_BUCKET, Key=k, Body=b"x")

    original_delete_objects = mock_s3_client.delete_objects
    call_sizes: list[int] = []

    def spy_delete_objects(**kwargs: Any) -> Any:
        call_sizes.append(len(kwargs["Delete"]["Objects"]))
        return original_delete_objects(**kwargs)

    mock_s3_client.delete_objects = spy_delete_objects

    errors = delete_objects(TEST_BUCKET, keys)

    assert errors == []
    assert call_sizes == [1000, 5]
    for k in keys:
        assert object_exists(f"{TEST_BUCKET}/{k}") is False


@pytest.mark.s3
def test_delete_objects_pass_exact_multiple_of_batch_size(mock_s3_client: S3Client) -> None:
    """Exactly 1000 keys results in a single batched delete_objects call."""
    total_keys = 1000
    keys = [f"exact_batch/file_{i:04d}.txt" for i in range(total_keys)]
    for k in keys:
        mock_s3_client.put_object(Bucket=TEST_BUCKET, Key=k, Body=b"x")

    original_delete_objects = mock_s3_client.delete_objects
    call_sizes: list[int] = []

    def spy_delete_objects(**kwargs: Any) -> Any:
        call_sizes.append(len(kwargs["Delete"]["Objects"]))
        return original_delete_objects(**kwargs)

    mock_s3_client.delete_objects = spy_delete_objects

    errors = delete_objects(TEST_BUCKET, keys)

    assert errors == []
    assert call_sizes == [1000]
    for k in keys:
        assert object_exists(f"{TEST_BUCKET}/{k}") is False


@pytest.mark.s3
def test_delete_objects_pass_propagates_per_key_errors_from_response(mock_s3_client: S3Client) -> None:
    """Per-key errors returned in the S3 response's Errors field are surfaced in the returned list."""
    fake_error = {"Key": "bad/key.txt", "Code": "AccessDenied", "Message": "Access Denied"}
    original_delete_objects = mock_s3_client.delete_objects

    def mock_delete_objects(**kwargs: Any) -> Any:
        resp = original_delete_objects(**kwargs)
        resp["Errors"] = [fake_error]
        return resp

    mock_s3_client.delete_objects = mock_delete_objects
    errors = delete_objects(TEST_BUCKET, ["some/key.txt"])
    assert errors == [fake_error]


@pytest.mark.s3
def test_delete_objects_removes_all(mock_s3_client: S3Client) -> None:
    """delete_objects removes every listed key in a single call."""
    keys = ["bulk/a.txt", "bulk/b.txt", "bulk/c.txt"]
    for k in keys:
        mock_s3_client.put_object(Bucket=TEST_BUCKET, Key=k, Body=b"data")

    errors = delete_objects(TEST_BUCKET, keys)

    assert errors == []
    for k in keys:
        assert object_exists(f"{TEST_BUCKET}/{k}") is False


@pytest.mark.s3
def test_delete_objects_empty_list_is_noop(mock_s3_client: S3Client) -> None:
    """delete_objects with an empty list makes no API call and returns no errors."""
    mock_s3_client.put_object(Bucket=TEST_BUCKET, Key="keep/me.txt", Body=b"safe")
    errors = delete_objects(TEST_BUCKET, [])
    assert errors == []
    assert object_exists(f"{TEST_BUCKET}/keep/me.txt") is True


@pytest.mark.s3
@pytest.mark.usefixtures("mock_s3_client")
def test_delete_objects_nonexistent_keys_no_error() -> None:
    """Deleting keys that don't exist returns no errors (S3 delete is idempotent)."""
    errors = delete_objects(TEST_BUCKET, ["ghost/a.txt", "ghost/b.txt"])
    assert errors == []
