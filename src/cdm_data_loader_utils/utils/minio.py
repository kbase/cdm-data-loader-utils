"""Upload / download utilities for Minio."""

from pathlib import Path
from typing import Any

import boto3
import tqdm
from berdl_notebook_utils.berdl_settings import get_settings

from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger

S3_BUCKET = "cdm-lake"

# Get credentials from environment variables (automatically set in JupyterHub)

settings = get_settings()
logger = get_cdm_logger()


# TODO: cache me?
def get_s3_client():
    """Create an s3 client."""
    return boto3.client(
        "s3",
        endpoint_url=settings.MINIO_ENDPOINT_URL,
        aws_access_key_id=settings.MINIO_ACCESS_KEY,
        aws_secret_access_key=settings.MINIO_SECRET_KEY,
    )


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
        logger.warning("list_remote_dir_contents did not return all files in %s", remote_dir)
    return response["Contents"]


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

    if not object_name:
        object_name = file_path.name

    file_size = file_path.stat().st_size

    s3 = get_s3_client()
    # Upload the file
    with tqdm.tqdm(total=file_size, unit="B", unit_scale=True, desc=str(file_path)) as pbar:
        try:
            # TODO: add in a check whether the obj already exists in the bucket using s3.head_object(Bucket, Key)
            s3.upload_file(
                Filename=str(file_path),
                Bucket=S3_BUCKET,
                Key=f"{destination_dir.removesuffix('/')}/{object_name}",
                Callback=lambda bytes_transferred: pbar.update(bytes_transferred),
            )
        except s3.exceptions.ClientError:
            get_cdm_logger().exception("Error uploading to s3")
            return False
        return True


def download_from_s3(bucket: str, key: str, filename: str, version_id: str | None = None) -> None:
    """
    Download an object from S3 with a progress bar.

    From https://alexwlchan.net/2021/04/s3-progress-bars/
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
            Callback=lambda bytes_transferred: pbar.update(bytes_transferred),
        )
