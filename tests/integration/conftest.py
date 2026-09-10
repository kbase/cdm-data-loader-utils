"""Shared fixtures and helpers for S3-backed integration tests.

Most tests use moto to mock S3. A subset of tests (marked `requires_ceph`)
still use a real CEPH instance for end-to-end validation.
"""

import hashlib
import re
from collections.abc import Generator
from functools import lru_cache
from pathlib import Path, PurePosixPath
from typing import Any
from unittest.mock import patch

import boto3
import botocore.config
import pytest
from botocore.exceptions import ClientError
from types_boto3_s3 import S3Client

from cdm_data_loaders.ncbi_ftp.assembly import build_accession_path
from cdm_data_loaders.utils.file_transfer.s3 import client
from cdm_data_loaders.utils.file_transfer.s3.client import _client_config, reset_s3_client

# Maximum length of a bucket name per S3/DNS spec
_MAX_BUCKET_LEN = 63


@lru_cache(maxsize=1)
def ceph_reachable() -> bool:
    """Return True if the CEPH endpoint accepts connections."""
    try:
        client = boto3.client(
            "s3",
            config=botocore.config.Config(
                connect_timeout=1,
                read_timeout=1,
                retries={"max_attempts": 1},
            ),
        )
        client.list_buckets()
    except Exception:  # noqa: BLE001
        return False
    return True


# Fixtures


@pytest.fixture
def ceph_s3_client() -> Generator[S3Client]:
    """Session-scoped real boto3 S3 client pointed at the local CEPH instance.

    Patches ``get_s3_client`` on every module that uses it so internal calls
    are transparently routed to CEPH.
    """
    if not ceph_reachable():
        pytest.skip("CEPH server not available")

    s3_client = boto3.client("s3", config=_client_config())
    reset_s3_client()
    with (
        patch.object(client, "get_s3_client", return_value=s3_client),
        patch.object(client, "_s3_client", s3_client),
    ):
        yield s3_client
    reset_s3_client()


@pytest.fixture
def s3_client(request: pytest.FixtureRequest) -> Generator[S3Client]:
    """Get the appropriate s3 client for the request -- either a moto mock or a ceph client."""
    backend = request.param
    if backend not in ("ceph", "mock"):
        err_msg = f"Invalid s3 client: {backend}"
        raise ValueError(err_msg)
    return request.getfixturevalue(f"{backend}_s3_client")


def _bucket_name_from_node(node_id: str, prefix: str | None = None) -> str:
    """Derive a DNS-compliant S3 bucket name from a pytest node ID.

    :param node_id: e.g. ``tests/integration/test_promote_e2e.py::test_dry_run``
    :param prefix: Optional prefix for the bucket name
    :return: e.g. ``integ-test-dry-run``
    """
    # Extract test function name from the node ID
    parts = node_id.split("::")
    name = parts[-1] if parts else node_id
    name = f"integ-{prefix}-{name}" if prefix else f"integ-{name}"
    # Lowercase, replace non-alphanumeric with hyphens, collapse multiples
    name = re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")
    if len(name) > _MAX_BUCKET_LEN:
        # Truncate but keep it unique via a short hash suffix
        suffix = hashlib.md5(name.encode()).hexdigest()[:6]  # noqa: S324
        name = f"{name[: _MAX_BUCKET_LEN - 7]}-{suffix}"
    return name


def check_existing_bucket(s3: S3Client, bucket: str) -> None:
    """Check whether a bucket exists, and create it if not. If it does exist, empty it.

    :param s3_client: boto3 S3 client
    :param bucket_name: name of the bucket to check
    """
    try:
        s3.head_bucket(Bucket=bucket)
        # Bucket exists — empty it for a clean run
        paginator = s3.get_paginator("list_objects_v2")
        for page in paginator.paginate(Bucket=bucket):
            for obj in page.get("Contents", []):
                s3.delete_object(Bucket=bucket, Key=obj["Key"])
    except s3.exceptions.NoSuchBucket:
        s3.create_bucket(Bucket=bucket)
    except ClientError as e:
        if e.response["Error"]["Code"] in ("404", "NoSuchBucket"):
            s3.create_bucket(Bucket=bucket)
        else:
            raise


@pytest.fixture
def test_bucket(s3_client: S3Client, request: pytest.FixtureRequest) -> PurePosixPath:
    """Create a per-test-method bucket and return its name."""
    bucket = _bucket_name_from_node(request.node.nodeid)
    check_existing_bucket(s3_client, bucket)
    return PurePosixPath(bucket)


@pytest.fixture
def staging_test_bucket(s3_client: S3Client, request: pytest.FixtureRequest) -> PurePosixPath:
    """Create a per-test staging bucket and return its name."""
    bucket = _bucket_name_from_node(request.node.nodeid, prefix="staging")
    check_existing_bucket(s3_client, bucket)
    return PurePosixPath(bucket)


# Helpers


def stage_files_to_s3(
    s3: S3Client,
    bucket: PurePosixPath,
    local_dir: Path,
    staging_prefix: PurePosixPath,
) -> list[PurePosixPath]:
    """Upload a local directory tree to an S3 staging prefix.

    :param s3: boto3 S3 client
    :param bucket: target bucket
    :param local_dir: local root directory to upload
    :param staging_prefix: S3 key prefix (e.g. ``"staging/run1/"``)
    :return: list of S3 keys uploaded
    """
    local_dir = Path(local_dir)
    keys: list[PurePosixPath] = []
    for path in sorted(local_dir.rglob("*")):
        if path.is_dir():
            continue
        rel = path.relative_to(local_dir)
        key = staging_prefix / rel
        s3.upload_file(Filename=str(path), Bucket=str(bucket), Key=str(key))
        keys.append(key)
    return keys


def seed_lakehouse(
    s3: S3Client,
    bucket: PurePosixPath,
    accession: str,
    files: dict[PurePosixPath, str | bytes],
    path_prefix: PurePosixPath,
    assembly_dir: PurePosixPath | None = None,
) -> list[PurePosixPath]:
    """Seed assembly files at the final Lakehouse path.

    :param s3: boto3 S3 client
    :param bucket: target bucket
    :param accession: assembly accession (e.g. ``"GCF_000001215.4"``)
    :param files: mapping of filename → content (str or bytes)
    :param path_prefix: Lakehouse prefix (e.g. ``"tenant-general-warehouse/…/ncbi/"``)
    :param assembly_dir: full assembly dir name; if None, uses ``accession``
    :return: list of S3 keys created
    """
    adir = assembly_dir or PurePosixPath(accession)
    rel = build_accession_path(adir)
    keys: list[PurePosixPath] = []
    prefix = PurePosixPath(path_prefix)
    for fname, content in files.items():
        key = prefix / rel / fname
        body = content.encode() if isinstance(content, str) else content
        md5 = hashlib.md5(body).hexdigest()  # noqa: S324
        s3.put_object(Bucket=str(bucket), Key=str(key), Body=body, Metadata={"md5": md5})
        keys.append(key)
    return keys


def list_all_keys(s3: S3Client, bucket: PurePosixPath, prefix: PurePosixPath | None = None) -> list[PurePosixPath]:
    """List all object keys in a bucket under a prefix.

    :param s3: boto3 S3 client
    :param bucket: bucket name
    :param prefix: optional key prefix filter
    :return: sorted list of keys
    """
    keys: list[PurePosixPath] = []
    paginator = s3.get_paginator("list_objects_v2")
    paginate_kwargs: dict[str, str] = {"Bucket": str(bucket)}
    if prefix is not None:
        paginate_kwargs["Prefix"] = str(prefix)
    for page in paginator.paginate(**paginate_kwargs):
        keys.extend(PurePosixPath(obj["Key"]) for obj in page.get("Contents", []))
    return sorted(keys)


def get_object_metadata(s3: S3Client, bucket: PurePosixPath, key: PurePosixPath) -> dict[str, Any]:
    """Return the S3 user metadata dict for an S3 object (from HeadObject).

    :param s3: boto3 S3 client
    :param bucket: bucket name
    :param key: object key
    :return: user metadata dict
    """
    resp = s3.head_object(Bucket=str(bucket), Key=str(key))
    return resp.get("Metadata", {})
