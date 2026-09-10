"""End-to-end tests for Phase 2 — FTP download of assemblies.

These tests download real (small) assemblies from the NCBI FTP server.
Marked ``requires_ceph`` (if a running CEPH test store is required) and
``slow_test``.
"""

import json
from functools import lru_cache
from pathlib import Path, PurePosixPath

import pytest
from types_boto3_s3 import S3Client

from cdm_data_loaders.ncbi_ftp.manifest import (
    compute_diff,
    download_assembly_summary,
    filter_by_prefix_range,
    parse_assembly_summary,
    write_transfer_manifest,
)
from cdm_data_loaders.pipelines.ncbi_ftp_download import download_and_stage, download_batch

# Use same stable prefix as manifest tests
STABLE_PREFIX = "900"


@lru_cache
def _manifest_for_one_assembly(tmp_path: Path) -> tuple[Path, str]:
    """Create a transfer manifest containing exactly one FTP path.

    Returns ``(manifest_path, accession)`` for the first latest assembly
    in the stable prefix range.
    """
    raw = download_assembly_summary(database="refseq")
    full = parse_assembly_summary(raw)
    filtered = filter_by_prefix_range(full, prefix_from=STABLE_PREFIX, prefix_to=STABLE_PREFIX)
    diff = compute_diff(filtered, previous_assemblies=None)

    assert len(diff.new) > 0, f"No new assemblies in prefix {STABLE_PREFIX}"

    manifest_path = tmp_path / "transfer_manifest.txt"
    write_transfer_manifest(diff, filtered, manifest_path)

    return manifest_path, diff.new[0]


@pytest.mark.slow_test
@pytest.mark.external_request
def test_download_small_batch(tmp_path: Path) -> None:
    """Download one assembly and verify directory structure and report."""
    manifest_path, _acc = _manifest_for_one_assembly(tmp_path)

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    report = download_batch(
        manifest_path=manifest_path,
        output_dir=output_dir,
        threads=1,
        limit=1,
    )

    assert report["succeeded"] >= 1
    assert report["failed"] == 0

    # Verify directory structure exists
    raw_data = output_dir / "raw_data"
    assert raw_data.exists(), "Expected raw_data/ directory in output"

    # Should have at least one assembly directory with files
    assembly_dirs = list(raw_data.rglob("GCF_*"))
    assert len(assembly_dirs) > 0, "Expected at least one assembly directory"

    # Check for .md5 sidecar files
    md5_files = list(raw_data.rglob("*.md5"))
    assert len(md5_files) > 0, "Expected .md5 sidecar files"

    # Check download report
    report_file = output_dir / "download_report.json"
    assert report_file.exists()
    with report_file.open() as f:
        saved_report = json.load(f)
    assert saved_report["succeeded"] >= 1


@pytest.mark.slow_test
@pytest.mark.external_request
def test_download_resume(tmp_path: Path) -> None:
    """Re-running download on the same manifest succeeds without errors.

    Verify download handles re-runs when some files are already present.
    """
    manifest_path, _acc = _manifest_for_one_assembly(tmp_path)

    output_dir = tmp_path / "output"
    output_dir.mkdir()

    # First download
    report1 = download_batch(
        manifest_path=manifest_path,
        output_dir=output_dir,
        threads=1,
        limit=1,
    )
    assert report1["succeeded"] >= 1

    files_after_first = set(output_dir.rglob("*"))

    # Second download — same manifest
    report2 = download_batch(
        manifest_path=manifest_path,
        output_dir=output_dir,
        threads=1,
        limit=1,
    )

    # Should succeed without errors (files overwritten or skipped)
    assert report2["succeeded"] >= 1
    assert report2["failed"] == 0

    # All original files should still exist
    files_after_second = set(output_dir.rglob("*"))
    assert files_after_first.issubset(files_after_second)


@pytest.mark.requires_ceph
@pytest.mark.slow_test
@pytest.mark.external_request
@pytest.mark.parametrize("s3_client", ["ceph", "mock"], indirect=True)
def test_download_and_stage_e2e(
    tmp_path: Path,
    s3_client: S3Client,
    test_bucket: PurePosixPath,
) -> None:
    """Download one assembly and verify it is staged under the expected S3 prefix."""
    manifest_path, _acc = _manifest_for_one_assembly(tmp_path)

    staging_prefix = PurePosixPath("staging") / "e2e-test"

    # Seed the manifest in CEPH so download_and_stage can read it from S3
    manifest_s3_key = staging_prefix / "input" / "transfer_manifest.txt"
    s3_client.put_object(
        Bucket=str(test_bucket),
        Key=str(manifest_s3_key),
        Body=manifest_path.read_bytes(),
    )

    report = download_and_stage(
        bucket=test_bucket,
        staging_key_prefix=staging_prefix,
        manifest_s3_key=manifest_s3_key,
        threads=1,
        limit=1,
        dry_run=False,
    )

    assert report["succeeded"] >= 1
    assert report["failed"] == 0
    assert report["staged_objects"] > 0
    assert report["staging_key_prefix"] == staging_prefix
    assert report["dry_run"] is False

    # Verify raw_data/ files and .md5 sidecars are staged
    paginator = s3_client.get_paginator("list_objects_v2")
    staged_keys = [
        obj["Key"]
        for page in paginator.paginate(Bucket=str(test_bucket), Prefix=str(staging_prefix / "raw_data"))
        for obj in page.get("Contents", [])
    ]
    assert len(staged_keys) > 0, "Expected staged files under raw_data/"

    data_files = [k for k in staged_keys if not k.endswith(".md5")]
    md5_files = [k for k in staged_keys if k.endswith(".md5")]
    assert len(data_files) > 0, "Expected data files"
    assert len(md5_files) > 0, "Expected .md5 sidecar files"

    # Verify download_report.json was also uploaded
    report_key = staging_prefix / "download_report.json"
    resp = s3_client.get_object(Bucket=str(test_bucket), Key=str(report_key))
    saved_report = json.loads(resp["Body"].read())
    assert saved_report["succeeded"] >= 1
