"""Utilities for s3 interaction."""

from contextlib import suppress
from logging import Logger, getLogger
from pathlib import Path
from types import ModuleType
from typing import IO, Any, Final

from botocore.exceptions import ClientError
from frozendict import frozendict

from cdm_data_loaders.utils.file_transfer.s3 import client
from cdm_data_loaders.utils.file_transfer.s3.transfer_config import build_transfer_config
from cdm_data_loaders.utils.progress import SynchronizedCallback, make_progress_bar

DEFAULT_EXTRA_ARGS: frozendict[str, str] = frozendict({"ChecksumAlgorithm": "CRC64NVME"})

VALID_S3_PREFIXES: list[str] = ["s3://", "s3a://"]

CDM_LAKE_BUCKET: Final[str] = "cdm-lake"
CTS_BUCKET: Final[str] = "cts"
VALID_BUCKETS: list[str] = [CDM_LAKE_BUCKET, CTS_BUCKET]

SUCCESS_RESPONSE: Final[int] = 200

NOT_FOUND_ERROR_CODES: frozenset[str] = frozenset({"404", "NoSuchKey", "NotFound", "NoSuchBucket"})

logger: Logger = getLogger(__name__)


def split_s3_path(s3_path: str, *, allow_bucket_only: bool = False) -> tuple[str, str | None]:
    """Convert a full s3 path (including bucket) into a bucket and key pair.

    Returns a tuple of bucket, key

    :param s3_path: an s3 path, including the bucket name (`s3://bucket/key` or `bucket/key`)
    :type s3_path: str
    :param allow_bucket_only: Allow parsing of a path that only includes a bucket (no key)
    :type allow_bucket_only: bool
    :return: tuple of (bucket, key)
    :rtype: tuple[str, str | None]
    """
    if "://" in s3_path:
        # remove the protocol prefix
        (_, unprefixed_path) = s3_path.split("://", 1)
    else:
        unprefixed_path = s3_path

    if not unprefixed_path:
        # raises a value error
        err_msg = f"Invalid path: '{s3_path}\nNo path found"
        raise ValueError(err_msg)

    if unprefixed_path.startswith("/"):
        err_msg = f"Invalid path: '{s3_path}'\ns3 paths must start with the bucket name"
        raise ValueError(err_msg)

    path_parts = unprefixed_path.split("/", 1)
    # return just the bucket if that is all that was passed
    # allow s3 paths like:
    #   s3://bucket
    #   s3://bucket/
    if allow_bucket_only and (len(path_parts) == 1 or (len(path_parts) == 2 and not path_parts[1])):  # noqa: PLR2004
        return (path_parts[0], "")

    # the first part should be the bucket and the second part the key
    if len(path_parts) != 2 or not path_parts[1]:  # noqa: PLR2004
        err_msg = f"Invalid path: '{s3_path}'\nCould not parse out bucket and key"
        raise ValueError(err_msg)

    return (path_parts[0], path_parts[1])


def list_objects(s3_path: str, *, max_keys: int = 1000) -> list[dict[str, Any]]:
    """List the remote paths that start with ``s3_path``.

    Note: since s3 paths are basically cosmetic, this function returns all paths that start with
    ``s3_path`` minus the bucket name.

    Retrieves all objects under the given prefix; collects all pages of results if there are more than
    1000 files (the max retrievable per `list_objects_v2` query) present.

    :param s3_path: directory to be listed, including the bucket name
    :type s3_path: str
    :param max_keys: maximum number of keys to return. boto3 defaults to 1000 records max.
    :type max_keys: int
    :return: list of object metadata dicts in the directory
    :rtype: list[dict[str, Any]]
    """
    s3 = client.get_s3_client()
    (bucket, key) = split_s3_path(s3_path)
    paginator = s3.get_paginator("list_objects_v2")
    page_iterator = paginator.paginate(Bucket=bucket, Prefix=key, MaxKeys=max_keys)

    contents = []
    for page in page_iterator:
        contents.extend(page.get("Contents", []))

    return contents


def head_object(s3_path: str) -> dict[str, Any]:
    """Check whether an object exists on s3.

    :param s3_path: path to the object on s3, including the bucket name
    :type s3_path: str
    :return: response from the head_object request
    :rtype: dict[str, Any]
    """
    s3 = client.get_s3_client()
    (bucket, key) = split_s3_path(s3_path)
    return s3.head_object(Bucket=bucket, Key=key, ChecksumMode="ENABLED")  # type: ignore[call-arg]


def get_existing_object_info(s3_path: str) -> dict[str, Any] | None:
    """Fetch metadata for an object already in S3, if it exists.

    :param s3_path: path to the object on s3, including the bucket name
    :type s3_path: str
    :raises ClientError: for any S3 error other than "object not found"
    :return: metadata for the existing object, or None if no object exists at `bucket`/`key`
    :rtype: dict[str, Any] | None
    """
    try:
        response = head_object(s3_path)
    except ClientError as exc:
        # AWS service errors, client-side issues
        logger.exception("Error checking for existing object: %s", s3_path)
        error_code = exc.response.get("Error", {}).get("Code")
        status_code = exc.response.get("ResponseMetadata", {}).get("HTTPStatusCode")
        if error_code in NOT_FOUND_ERROR_CODES or status_code == 404:  # noqa: PLR2004
            return None
        raise
    except Exception:
        # any other errors - e.g. invalid/missing args
        logger.exception("Error checking for existing object: %s", s3_path)
        raise

    return response


def object_exists(s3_path: str) -> bool:
    """Check whether an object exists on s3.

    :param s3_path: path to the object on s3, including the bucket name
    :type s3_path: str
    :return: True if the object exists, False otherwise
    :rtype: bool
    """
    return get_existing_object_info(s3_path) is not None


def check_existing_objects(s3_path: str, object_info: dict[str, Any]) -> bool:
    """Check if an object exists in S3 and raise an error if it does.

    :param s3_path: path to the object on s3, including the bucket name
    :type s3_path: str
    :param object_info: information about the object to check
    :type object_info: dict[str, Any]
    :return: True if the object exists, False otherwise
    :rtype: bool
    """
    existing_object = None

    # suppress exceptions so that 404s and other errors are not triggered
    with suppress(Exception):
        existing_object = get_existing_object_info(s3_path)

    if existing_object is None:
        return False

    # compare checksum, metadata, and other attributes if needed to determine if the object is the same or different
    if object_info.get("Metadata") != existing_object["metadata"]:
        logger.debug("Object metadata differs for %s", s3_path)
        return False

    # check the file sizes
    if object_info.get("ContentLength") != existing_object["size"]:
        logger.debug("Object size differs for %s", s3_path)
        return False

    # TODO: all data uploaded to s3 should have a checksum
    checksum_data = {k: v.lower() for k, v in existing_object.items() if k.startswith("Checksum")}
    # if checksum_data == object_info["checksum"]:
    #     logger.debug("Object checksum matches for %s", s3_path)
    #     return True

    return False


def upload_fileobj(  # noqa: PLR0913
    fileobj: IO[bytes],
    s3_path: str,
    *,
    file_path: Path | None = None,
    file_size: int | None = 0,
    user_metadata: dict[str, str] | None = None,
    transfer_config_kwargs: dict[str, Any] | None = None,
    skip_if_exists: bool = False,
    show_progress: bool = True,
) -> bool:
    """Upload data from a binary file-like object to S3.

    Can be used for file uploads via the :func:`upload_file` interface
    or directly for binary data, e.g. data streamed via FTP/HTTP.

    :param fileobj: a binary-mode file-like object supporting ``read()``
    :type fileobj: IO[bytes]
    :param s3_path: s3 destination path, including the bucket name
    :type s3_path: str
    :param file_path: if the input is a file, a Path object for the file (used for calculating checksums, size, etc.)
    :type file_path: Path | None
    :param file_size: size of the file in bytes
    :type file_size: int | None, defaults to 0
    :param user_metadata: user metadata key/value pairs to attach to the object
    :type user_metadata: dict[str, str] | None
    :param transfer_config_kwargs: keyword arguments for configuring the S3 transfer
    :type transfer_config_kwargs: dict[str, Any] | None
    :param skip_if_exists: whether to skip the upload if the object already exists
    :type skip_if_exists: bool
    :param show_progress: whether to display a tqdm progress bar during upload
    :type show_progress: bool
    :return: True if the upload succeeded, else False
    :rtype: bool
    """
    s3 = client.get_s3_client()
    (bucket, key) = split_s3_path(s3_path)

    if skip_if_exists:
        logger.debug("To be implemented")

    upload_args = {
        "Fileobj": fileobj,
        "Bucket": bucket,
        "Key": key,
        "ExtraArgs": {
            **DEFAULT_EXTRA_ARGS,
            **(({"Metadata": user_metadata}) if user_metadata is not None else {}),
        },
        "Config": build_transfer_config(file_size or 0, **(transfer_config_kwargs or {})),
    }

    if file_path:
        logger.debug("uploading %s to %s", str(file_path), s3_path)
    else:
        logger.debug("uploading fileobj to %s", s3_path)
    try:
        with make_progress_bar(total=file_size, desc=s3_path, disable=not show_progress) as pbar:
            s3.upload_fileobj(
                **upload_args,
                Callback=SynchronizedCallback(pbar.update),
            )
    except Exception:
        logger.exception("Error uploading to s3")
        return False
    return True


def upload_file(  # noqa: PLR0913
    local_file_path: Path | str,
    destination_dir: str,
    object_name: str | None = None,
    *,
    user_metadata: dict[str, str] | None = None,
    transfer_config_kwargs: dict[str, Any] | None = None,
    skip_if_exists: bool = False,
    show_progress: bool = True,
) -> bool:
    """Upload a local file to an S3 bucket.

    Internally opens the file and delegates to :func:`upload_fileobj`.

    :param local_file_path: File to upload
    :type local_file_path: Path | str
    :param destination_dir: path to the destination directory on s3, including the bucket name and EXCLUDING the file name
    :type destination_dir: str
    :param object_name: S3 object name. If not specified, the name of the file from local_file_path is used.
    :type object_name: str | None
    :param user_metadata: user metadata key/value pairs to attach to the object; when provided the upload always runs
    :type user_metadata: dict[str, str] | None
    :param transfer_config_kwargs: keyword arguments for configuring the S3 transfer
    :type transfer_config_kwargs: dict[str, Any] | None
    :param skip_if_exists: whether to skip the upload if the object already exists; not yet implemented
    :type skip_if_exists: bool
    :param show_progress: whether to display a tqdm progress bar during upload, defaults to True
    :type show_progress: bool, optional
    :return: True if file was uploaded, else False
    :rtype: bool
    """
    if not local_file_path:
        err_msg = "No local file path specified"
        raise ValueError(err_msg)

    if isinstance(local_file_path, str):
        local_file_path = Path(local_file_path)

    if not object_name:
        object_name = local_file_path.name
    if not object_name:
        msg = "No object_name supplied for the upload"
        raise ValueError(msg)

    if not destination_dir:
        msg = "No destination directory supplied for the file"
        raise ValueError(msg)

    s3_path = f"{destination_dir.removesuffix('/')}/{object_name}"

    with local_file_path.open(mode="rb") as fh:
        return upload_fileobj(
            fh,
            s3_path,
            file_path=local_file_path,
            file_size=local_file_path.stat().st_size,
            user_metadata=user_metadata,
            transfer_config_kwargs=transfer_config_kwargs,
            skip_if_exists=skip_if_exists,
            show_progress=show_progress,
        )


def stream_to_s3(url: str, s3_path: str, requests: ModuleType) -> str:
    """Stream directly from an HTTP download to s3.

    :param url: address of the object to transfer to s3
    :type url: str
    :param s3_path: save path on s3
    :type s3_path: str
    :param requests: module implementing requests.get and returning a response
    :type requests: ModuleType
    :return: path of the file on s3, in the form bucket/key
    :rtype: str
    """
    s3 = client.get_s3_client()
    (bucket, key) = split_s3_path(s3_path)
    with requests.get(url, stream=True) as response:
        response.raise_for_status()
        s3.upload_fileobj(
            # raw stream from urllib3
            response.raw,
            bucket,
            key,
            ExtraArgs={
                **DEFAULT_EXTRA_ARGS,
                "ContentType": response.headers.get("content-type", "application/octet-stream"),
            },
        )
    return f"{bucket}/{key}"


def upload_dir(
    local_dir: Path | str,
    destination_dir: str,
    file_glob: str | None = None,
) -> bool:
    """Upload a directory to an s3 bucket.

        If file_glob is not set, it will default to "**/*", i.e. any path with at least
        one segment (recursive listing of all paths from the current directory).

        Wildcards:
        ** (entire segment)
            Matches any number of file or directory segments, including zero.
            "assets/**" matches any path starting with "assets/", including "assets/"

        * (entire segment)
            Matches one file or directory segment.
            "assets/*" matches any file or directory under "assets/" but none of the children

        **/* (two segments)
            Matches any path with at least one segment
            "assets/**/*" matches any file or directory under "assets/" but not "assets/" itself

        * (part of a segment)
            Matches any number of non-separator characters, including zero.
            "file*.txt" would match "file.txt", "file_type.txt", "file12345.txt", but not "file/b.txt"

        ? (part of a segment)
            Matches one non-separator character.
            "file?.txt" would match "filea.txt" or "file1.txt"

    [seq]
        Matches one character in seq, where seq is a sequence of characters. Range expressions are supported; for example, [a-z] matches any lowercase ASCII letter. Multiple ranges can be combined: [a-zA-Z0-9_] matches any ASCII letter, digit, or underscore.

    [!seq]
        Matches one character not in seq, where seq follows the same rules as above.

    For a literal match, wrap the meta-characters in brackets. For example, "[?]" matches the character "?".

    :param local_dir_path: local directory to upload
    :type local_dir_path: Path | str
    :param destination_dir: remote directory to upload to, including the bucket name
    :type destination_dir: str
    :param file_glob: glob for selecting files to upload
    :type file_glob: str | None
    :return: True or False, depending on the result of the uploads
    :rtype: bool
    """
    if not local_dir:
        msg = "No source directory supplied for the upload"
        raise ValueError(msg)

    if not destination_dir:
        msg = "No destination directory supplied for the upload"
        raise ValueError(msg)

    if not file_glob:
        file_glob = "**/*"

    if isinstance(local_dir, str):
        local_dir = Path(local_dir)

    all_successful = True
    for path in sorted(local_dir.glob(file_glob)):
        if path.is_dir():
            continue
        # get the path of the current file relative to local_dir and use that as the object name
        success = upload_file(path, destination_dir, object_name=str(path.relative_to(local_dir)))
        if not success:
            all_successful = False

    return all_successful


def download_file(
    s3_path: str,
    local_file_path: str | Path,
    transfer_config_kwargs: dict[str, Any] | None = None,
    show_progress: bool = True,  # noqa: FBT001, FBT002
) -> None:
    """Download an object from s3.

    WARNING: will overwrite existing files but will not overwrite a file whilst trying to make a directory

    Will attempt to create the local directory if it does not exist.

    :param s3_path: path to the file on s3, including the bucket name
    :type s3_path: str
    :param local_file_path: local path (including file name) to save the downloaded file to
    :type local_file_path: str | Path
    :param transfer_config_kwargs: keyword arguments for configuring the S3 transfer
    :type transfer_config_kwargs: dict[str, Any] | None
    :param show_progress: whether to display a tqdm progress bar during upload, defaults to True
    :type show_progress: bool, optional
    """
    local_file_path = Path(local_file_path)
    # check whether the parent directory exists
    parent_dir = local_file_path.parent
    if not parent_dir.is_dir():
        try:
            parent_dir.mkdir(parents=True, exist_ok=False)
        except Exception:
            logger.exception("Could not save s3 file to %s", local_file_path)
            raise

    s3 = client.get_s3_client()
    (bucket, key) = split_s3_path(s3_path)

    download_args = {
        "Filename": str(local_file_path),
        "Bucket": bucket,
        "Key": key,
        "ExtraArgs": {"ChecksumMode": "ENABLED"},
    }

    # Get the object size
    existing = get_existing_object_info(s3_path)
    if existing is None:
        msg = f"File not found: {s3_path}"
        raise FileNotFoundError(msg)
    object_size = existing.get("ContentLength", 0)
    download_args["Config"] = build_transfer_config(object_size, **(transfer_config_kwargs or {}))

    with make_progress_bar(total=object_size, desc=str(local_file_path), disable=not show_progress) as pbar:
        s3.download_file(
            **download_args,
            Callback=SynchronizedCallback(pbar.update),
        )


def copy_object(
    current_s3_path: str,
    new_s3_path: str,
) -> dict[str, Any]:
    """Copy an object from one place to another, inheriting the source user metadata.

    Source user metadata (e.g. ``md5``) is preserved on the destination because
    ``MetadataDirective`` is omitted, which defaults to ``COPY``.

    A successful copy operation will return a response where
    resp["ResponseMetadata"]["HTTPStatusCode"] == 200

    Errors (e.g, buckets or keys not existing, wrong credentials, etc.) are passed
    directly to the user without being caught.

    :param current_s3_path: path to the file on s3, including the bucket name
    :type current_s3_path: str
    :param new_s3_path: the desired new file path on s3, including the bucket name
    :type new_s3_path: str
    :return: dictionary containing response from the copy operation
    :rtype: dict[str, Any]
    """
    s3 = client.get_s3_client()
    (current_s3_bucket, current_s3_key) = split_s3_path(current_s3_path)
    (new_s3_bucket, new_s3_key) = split_s3_path(new_s3_path)

    return s3.copy_object(
        CopySource={"Bucket": current_s3_bucket, "Key": current_s3_key},
        Bucket=new_s3_bucket,
        Key=new_s3_key,
        **DEFAULT_EXTRA_ARGS,
    )


def copy_directory(current_s3_path: str, new_s3_path: str) -> tuple[dict[str, str], dict[str, Any]]:
    """Copy all objects under a given S3 prefix to a new prefix.

    Preserves the relative key structure under the source prefix. For example,
    copying s3://my-bucket/foo/ to s3://my-bucket/bar/ will copy
    s3://my-bucket/foo/a/b.txt -> s3://my-bucket/bar/a/b.txt

    If the source bucket does not exist, a NoSuchBucket error will be thrown.

    :param current_s3_path: path to the directory on s3, including the bucket name
    :type current_s3_path: str
    :param new_s3_path: the desired new directory path on s3, including the bucket name
    :type new_s3_path: str
    :return: a tuple of (successes, errors) where:
             - successes maps "bucket/source_key" -> "bucket/dest_key" for each
               successfully copied object
             - errors maps "bucket/source_key" -> the exception or response object
               for each failed copy
    :rtype: tuple[dict[str, str], dict[str, Any]]
    """
    s3 = client.get_s3_client()
    (current_s3_bucket, current_s3_prefix) = split_s3_path(current_s3_path)
    (new_s3_bucket, new_s3_prefix) = split_s3_path(new_s3_path)

    if current_s3_prefix and not current_s3_prefix.endswith("/"):
        current_s3_prefix += "/"
    if new_s3_prefix and not new_s3_prefix.endswith("/"):
        new_s3_prefix += "/"

    paginator = s3.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=current_s3_bucket, Prefix=current_s3_prefix)

    successes: dict[str, str] = {}
    errors: dict[str, Any] = {}

    for page in pages:
        for obj in page.get("Contents", []):
            current_key = obj["Key"]
            relative_key = current_key[len(current_s3_prefix) :]
            new_key = new_s3_prefix + relative_key

            source_path = f"{current_s3_bucket}/{current_key}"
            dest_path = f"{new_s3_bucket}/{new_key}"

            try:
                resp = s3.copy_object(
                    CopySource={"Bucket": current_s3_bucket, "Key": current_key},
                    Bucket=new_s3_bucket,
                    Key=new_key,
                    **DEFAULT_EXTRA_ARGS,
                )
                if resp["ResponseMetadata"]["HTTPStatusCode"] == SUCCESS_RESPONSE:
                    successes[source_path] = dest_path
                else:
                    errors[source_path] = resp
            except Exception as e:
                logger.exception("Failed to copy %s to %s", source_path, dest_path)
                errors[source_path] = e

    return successes, errors


def delete_object(s3_path: str) -> dict[str, Any]:
    """Delete an object from s3.

    A successful deletion will return a response where
    resp["ResponseMetadata"]["HTTPStatusCode"] == 204.

    Errors (e.g, buckets or keys not existing, wrong credentials, etc.) are passed
    directly to the user without being caught.

    :param s3_path: path to the file on s3, including the bucket name
    :type s3_path: str
    :return: dictionary containing response
    :rtype: dict[str, Any]
    """
    s3 = client.get_s3_client()
    (bucket, key) = split_s3_path(s3_path)
    return s3.delete_object(Bucket=bucket, Key=key)


def delete_objects(bucket: str, keys: list[str]) -> list[dict[str, Any]]:
    """Delete multiple objects from an S3 bucket in a single API call.

    Splits into batches of 1000 (the S3 API maximum per request).

    :param bucket: S3 bucket name (no protocol prefix)
    :param keys: list of S3 keys to delete
    :return: list of per-key error dicts returned by S3 (empty if all succeeded)
    :rtype: list[dict[str, Any]]
    """
    if not keys:
        return []

    s3 = client.get_s3_client()
    errors: list[dict[str, Any]] = []
    for i in range(0, len(keys), 1000):
        batch = keys[i : i + 1000]
        resp = s3.delete_objects(
            Bucket=bucket,
            Delete={"Objects": [{"Key": k} for k in batch], "Quiet": False},
        )
        errors.extend(resp.get("Errors", []))
    return errors
