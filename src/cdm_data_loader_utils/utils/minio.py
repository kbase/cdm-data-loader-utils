from pathlib import Path
import botocore
import boto3
import tqdm
from typing import Any
from functools import lru_cache


S3_BUCKET = "cdm-lake"
DEFAULT_EXTRA_ARGS = {"ChecksumAlgorithm": "CRC64NVME"}


@lru_cache
def get_s3_client(args: dict[str, str] | None = None) -> botocore.client.BaseClient:
    if not args:
        try:
            # FIXME: make this less clunky!
            from berdl_notebook_utils.berdl_settings import get_settings

            settings = get_settings()
            args = {
                "endpoint_url": settings.MINIO_ENDPOINT_URL,
                "aws_access_key_id": settings.MINIO_ACCESS_KEY,
                "aws_secret_access_key": settings.MINIO_SECRET_KEY,
            }
        except (ModuleNotFoundError, ImportError, NameError) as e:
            # not in an environment where `get_settings` is defined
            print(e)
            raise
        except Exception:
            raise

    required_args = ["endpoint_url", "aws_access_key_id", "aws_secret_access_key"]
    keyword_args = {kw: args.get(kw) for kw in required_args}
    missing = [kw for kw in required_args if not keyword_args[kw]]
    if missing:
        msg = "Cannot initialise s3 client: missing arguments: " + ", ".join(missing)
        raise ValueError(msg)

    # create the S3 client
    return boto3.client("s3", **keyword_args)


def list_remote_dir_contents(remote_dir: str) -> list[dict[str, Any]]:
    """List the contents of a remote directory.

    :param remote_dir: directory to be listed
    :type remote_dir: str
    :return: list of things in the directory
    :rtype: list[str]
    """
    s3 = get_s3_client()
    response = s3.list_objects_v2(Bucket=S3_BUCKET, Prefix=remote_dir)
    if response["IsTruncated"]:
        print(f"list_remote_dir_contents did not return all files in {remote_dir}")
    return response["Contents"]


def file_exists(path_to_s3_file: str) -> bool:
    s3 = get_s3_client()
    try:
        s3.head_object(Bucket=S3_BUCKET, Key=path_to_s3_file)
    except s3.exceptions.ClientError as e:
        # Error performing head operation on s3 object: An error occurred (404) when calling the HeadObject operation: Not Found
        error_string = str(e)
        if not error_string.startswith("An error occurred (404) when calling the HeadObject operation: Not Found"):
            print(f"Error performing head operation on s3 object: {e!s}")
        return False
    return True


def upload_file(file_path: Path | str, destination_dir: str, object_name: str | None = None) -> bool:
    """Upload a file to an S3 bucket.

    :param file_path: File to upload
    :type file_path: Path | str
    :param destination_dir: location within the cdm-lake bucket to upload to
    :type destination_dir: str
    :param object_name: S3 object name. If not specified, the name of the file from file_path is used.
    :type object_name: str | None
    :return: True if file was uploaded, else False
    :rtype: bool
    """
    if isinstance(file_path, str):
        file_path = Path(file_path)

    if not destination_dir:
        msg = "No destination directory supplied for the file"
        raise ValueError(msg)

    dest_dir = destination_dir.removeprefix(S3_BUCKET)
    dest_dir = dest_dir.removeprefix("/")
    dest_dir = dest_dir.removesuffix("/")

    if not object_name:
        object_name = file_path.name

    s3 = get_s3_client()
    path_to_s3_file = f"{dest_dir}/{object_name}"

    if file_exists(path_to_s3_file):
        print(f"File already present: {S3_BUCKET}/{path_to_s3_file}")
        return True

    # Upload the file
    file_size = file_path.stat().st_size
    with tqdm.tqdm(total=file_size, unit="B", unit_scale=True, desc=str(file_path)) as pbar:
        print(f"uploading {file_path!s} to {S3_BUCKET}/{path_to_s3_file}")
        try:
            # TODO: add in a check whether the obj already exists in the bucket using s3.head_object(Bucket, Key)
            s3.upload_file(
                Filename=str(file_path),
                Bucket=S3_BUCKET,
                Key=path_to_s3_file,
                Callback=pbar.update,
                ExtraArgs=DEFAULT_EXTRA_ARGS,
            )
        except s3.exceptions.ClientError as e:
            print(f"Error uploading to s3: {e!s}")
            return False
        return True


def upload_dir(dir_path: Path | str, destination_dir: str, file_glob: str | None = None) -> bool:
    """Upload a directory to an s3 bucket.

    :param dir_path: local directory to upload
    :type dir_path: Path | str
    :param destination_dir: remote directory to upload to
    :type destination_dir: str
    :param file_glob: glob for selecting files to upload
    :type file_glob: str | None
    :return: True or False, depending on the result of the uploads
    :rtype: bool
    """
    if not dir_path:
        msg = "No source directory supplied for the upload"
        raise ValueError(msg)

    if not destination_dir:
        msg = "No destination directory supplied for the upload"
        raise ValueError(msg)
    destination_dir = destination_dir.removesuffix("/")

    if not file_glob:
        file_glob = "*"

    if isinstance(dir_path, str):
        dir_path = Path(dir_path)

    all_successful = True
    for path in sorted(dir_path.glob("*")):
        if path.is_dir():
            success = upload_dir(path, f"{destination_dir}/{path.name}", file_glob)
        else:
            success = upload_file(path, destination_dir)

        if not success:
            all_successful = False

    return all_successful


def copy_file(current_path: str, new_path: str) -> bool:
    """Copy a file from one place to another, adding in a CRC64NVME checksum."""
    s3 = get_s3_client()
    current_path = current_path.removeprefix(S3_BUCKET)
    new_path = new_path.removeprefix(S3_BUCKET)

    # S3.Client.copy(CopySource, Bucket, Key, ExtraArgs=None, Callback=None, SourceClient=None, Config=None)
    return s3.copy_object(
        CopySource={"Bucket": S3_BUCKET, "Key": current_path},
        Bucket=S3_BUCKET,
        Key=new_path,
        ChecksumAlgorithm="CRC64NVME",
    )


def delete_file(file_path: str) -> bool:
    s3 = get_s3_client()
    file_path = file_path.removeprefix(S3_BUCKET)
    print(f"Deleting file: {file_path}")
    return s3.delete_object(Bucket=S3_BUCKET, Key=file_path)


def download_from_s3(bucket: str, key: str, filename: str, version_id: str | None = None) -> None:
    """
    Download an object from S3 with a progress bar.
    """
    s3 = get_s3_client()
    kwargs = {"Bucket": bucket, "Key": key}
    if version_id is not None:
        kwargs["VersionId"] = version_id

    # Get the object size
    object_size = s3.head_object(**kwargs)["ContentLength"]

    extra_args = {"VersionId": version_id} if version_id is not None else None

    # set ``unit_scale=True`` so tqdm uses SI unit prefixes
    # ``unit="B"`` means it adds the string "B" as a suffix
    # progress is reported as (e.g.) "14.5kB/s".
    with tqdm.tqdm(total=object_size, unit="B", unit_scale=True, desc=filename) as pbar:
        s3.download_file(
            Bucket=bucket,
            Key=key,
            ExtraArgs=extra_args,
            Filename=filename,
            Callback=pbar.update,
        )
