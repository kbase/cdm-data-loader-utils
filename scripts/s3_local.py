#!/usr/bin/env python3
# ruff: noqa: D103
"""Thin S3 CLI for local CEPH testing (no aws-cli install required).

Usage (all commands assume ``uv run`` from the repo root):

    uv run python scripts/s3_local.py mb  s3://cdm-lake
    uv run python scripts/s3_local.py cp  staging/raw_data/ s3://cdm-lake/staging/run1/raw_data/
    uv run python scripts/s3_local.py cp  my/local/file.txt s3://cdm-lake/staging/staged_file.txt
    uv run python scripts/s3_local.py ls  s3://cdm-lake/staging/run1/
    uv run python scripts/s3_local.py head s3://cdm-lake/some/key.gz

Environment variables (with defaults for the walkthrough):

    AWS_ENDPOINT_URL         http://localhost:9000
    AWS_ACCESS_KEY_ID        test_access_key
    AWS_SECRET_ACCESS_KEY    test_access_secret
"""

import os
import sys

from cdm_data_loaders.utils.file_transfer.s3 import cli, client


def _client() -> None:
    client.reset_s3_client()
    _ = client.get_s3_client(
        {
            "endpoint_url": os.environ.get("AWS_ENDPOINT_URL", "http://localhost:9000"),
            "aws_access_key_id": os.environ.get("AWS_ACCESS_KEY_ID", "test_access_key"),
            "aws_secret_access_key": os.environ.get("AWS_SECRET_ACCESS_KEY", "test_access_secret"),
        }
    )


# dispatch

COMMANDS = {"mb": cli.cmd_mb, "cp": cli.cmd_cp, "ls": cli.cmd_ls, "head": cli.cmd_head}


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] not in COMMANDS:  # noqa: PLR2004
        cmds = ", ".join(COMMANDS)
        err_msg = f"Usage: s3_local.py <{cmds}> [args ...]\n\n{__doc__}"
        raise SystemExit(err_msg)
    _client()
    COMMANDS[sys.argv[1]](sys.argv[2:])


if __name__ == "__main__":
    main()
