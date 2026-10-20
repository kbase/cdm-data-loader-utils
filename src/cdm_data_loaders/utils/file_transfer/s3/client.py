"""S3 client creation and reset."""

from typing import Final

import boto3
from botocore.config import Config
from types_boto3_s3.client import S3Client

# "legacy", "standard", "adaptive"
AWS_CLIENT_RETRY_MODE: Final[str] = "adaptive"
# how many times to retry, including the initial attempt
AWS_CLIENT_TOTAL_MAX_ATTEMPTS: Final[int] = 5

_s3_client: S3Client | None = None


def _client_config() -> Config:
    """Configuration for S3 client."""
    return Config(
        retries={"total_max_attempts": AWS_CLIENT_TOTAL_MAX_ATTEMPTS, "mode": AWS_CLIENT_RETRY_MODE},
        request_checksum_calculation="when_supported",
        response_checksum_validation="when_supported",
    )


def get_s3_client(args: dict[str, str | None] | None = None) -> S3Client:
    """Create an S3 client using the provided arguments.

    The client is created once and cached for subsequent calls.
    Call reset_s3_client() to force a new client to be created on the next call.

    To configure the client using arguments, provide a dictionary with the following keys:
        - aws_access_key_id: the access key ID for the S3 client
        - aws_secret_access_key: the secret access key for the S3 client
        - endpoint_url: the endpoint URL for the S3 client (e.g., "https://s3.amazonaws.com" or "https://my-s3-server.com")

    If arguments are not provided, the client will be created using boto3's default
    configuration method, which looks for environment variables (AWS_ACCESS_KEY_ID,
    AWS_SECRET_ACCESS_KEY, and AWS_ENDPOINT_URL_S3 or AWS_ENDPOINT_URL) or an ``./aws`` config directory.
    See the boto3 documentation for more details.

    :param args: arguments for creating the S3 client, defaults to None
    :type args: dict[str, str] | None, optional
    :raises ValueError: if required arguments for creating the S3 client are missing
    :return: initialised s3 client
    :rtype: S3Client
    """
    global _s3_client  # noqa: PLW0603
    if _s3_client is not None:
        return _s3_client

    if not args:
        args = {}

    valid_kwargs = ["aws_access_key_id", "aws_secret_access_key", "endpoint_url"]
    kwargs = {k: v for k, v in args.items() if k in valid_kwargs and v is not None}

    # make sure that if aws_access_key_id or aws_secret_access_key is provided, the other is also provided.
    if bool(kwargs.get("aws_access_key_id")) ^ bool(kwargs.get("aws_secret_access_key")):
        msg = "Cannot initialise s3 client: aws_access_key_id and aws_secret_access_key must be provided together, either via args or environment variables or a config file"
        raise ValueError(msg)

    # initialise using boto3's default config behaviour, plus any overrides from args
    client = boto3.client("s3", config=_client_config(), **kwargs)  # pyright: ignore[reportArgumentType, reportCallIssue]

    missing = []
    # boto3 will not raise an error on client creation if credentials are missing, so throw an error now
    credentials = client._request_signer._credentials  # noqa: SLF001
    if not credentials:
        missing = ["aws_access_key_id", "aws_secret_access_key"]
    else:
        if not credentials.access_key:
            missing.append("aws_access_key_id")
        if not credentials.secret_key:
            missing.append("aws_secret_access_key")

    if missing:
        msg = "Cannot initialise s3 client: missing configuration values: " + ", ".join(missing)
        raise ValueError(msg)

    # nothing missing: we are good to go!
    _s3_client = client
    return _s3_client


def reset_s3_client() -> None:
    """Reset the cached S3 client, forcing a new one to be created on the next call to get_s3_client."""
    global _s3_client  # noqa: PLW0603
    _s3_client = None
