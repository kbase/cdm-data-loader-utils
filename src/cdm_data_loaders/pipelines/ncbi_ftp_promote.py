"""NCBI FTP staged file promotion pipeline (Phase 3).

Promotes staged files from the containerized NCBI FTP download pipeline to their
final location in the Lakehouse.
"""

import logging
from pathlib import Path, PurePosixPath

from pydantic import AliasChoices, Field
from pydantic_settings import CliImplicitFlag

from cdm_data_loaders.core.settings import LoggerSettings
from cdm_data_loaders.ncbi_ftp.promote import promote_from_s3
from cdm_data_loaders.pipelines.core import PipelineCli

DEFAULT_STAGING_BUCKET: PurePosixPath = PurePosixPath("cts")
DEFAULT_DESTINATION_BUCKET: PurePosixPath = PurePosixPath("cdm-lake")
DEFAULT_DESTINATION_PREFIX: PurePosixPath = PurePosixPath("tenant-general-warehouse/kbase/datasets/ncbi")
DEFAULT_TRANSFER_MANIFEST_FILE: PurePosixPath = PurePosixPath("transfer_manifest.txt")
ESTIMATED_FILES_PER_ASSEMBLY: int = 21

logger = logging.getLogger("dlt")


class PromoteSettings(LoggerSettings):
    """Configuration for the NCBI FTP file promotion pipeline."""

    staging_bucket: PurePosixPath = Field(
        default=DEFAULT_STAGING_BUCKET,
        description="Bucket where staged files are located after download",
        validation_alias=AliasChoices("staging-bucket", "staging_bucket"),
    )
    destination_bucket: PurePosixPath = Field(
        default=DEFAULT_DESTINATION_BUCKET,
        description="Bucket to which staged files will be promoted.",
        validation_alias=AliasChoices("destination-bucket", "destination_bucket"),
    )
    staging_path: PurePosixPath = Field(
        description="Path to folder in the staging bucket where staged files are located; Should contain `raw_data/` folder",
        validation_alias=AliasChoices("s", "staging-path", "staging_path"),
    )
    destination_path: PurePosixPath = Field(
        default=DEFAULT_DESTINATION_PREFIX,
        description="Path to folder in the destination bucket where files will be promoted to; Will contain `raw_data/` folder",
        validation_alias=AliasChoices("destination-path", "destination_path"),
    )
    removed_manifest_path: Path | None = Field(
        default=None,
        description="Local filesystem path to the removed files manifest from Phase 1, or None to skip archiving removed records",
        validation_alias=AliasChoices("r", "removed-manifest", "removed_manifest"),
    )
    updated_manifest_path: Path | None = Field(
        default=None,
        description="Local filesystem path to the updated files manifest from Phase 1, or None to skip archiving updated records",
        validation_alias=AliasChoices("u", "updated-manifest", "updated_manifest"),
    )
    transfer_manifest_path: PurePosixPath | None = Field(
        default_factory=lambda settings: settings["staging_path"] / DEFAULT_TRANSFER_MANIFEST_FILE,
        description="S3 object key of the transfer manifest to trim after promotion, or None to skip pruning staged files",
        validation_alias=AliasChoices("t", "transfer-manifest", "transfer_manifest"),
    )
    dry_run: CliImplicitFlag[bool] = Field(
        default=False,
        description="Log actions without making changes",
        validation_alias=AliasChoices("dry-run", "dry_run"),
    )


def run_promote(config: PromoteSettings) -> None:
    """Main CTS entry point for Phase 3 promotion.

    :param config: validated promote settings
    """
    report = promote_from_s3(
        staging_bucket=config.staging_bucket,
        staging_key_prefix=config.staging_path,
        lakehouse_bucket=config.destination_bucket,
        lakehouse_key_prefix=config.destination_path,
        removed_manifest_path=config.removed_manifest_path,
        updated_manifest_path=config.updated_manifest_path,
        manifest_s3_key=config.transfer_manifest_path,
        dry_run=config.dry_run,
    )
    if report["failed"] > 0:
        msg = f"Promote completed with {report['failed']} failures"
        raise RuntimeError(msg)


cli = PipelineCli(PromoteSettings, run_promote)


if __name__ == "__main__":
    cli()
