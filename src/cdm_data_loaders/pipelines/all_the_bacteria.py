"""All The Bacteria bulk data downloader.

Project documentation: https://allthebacteria.org/

Batch downloading from OSF: https://allthebacteria.org/docs/osf_downloads/

all_atb_files.tsv: https://osf.io/xv7q9/files/r6gcp (or Rg6cp, casing varies)

# TODO: change the output to parquet files.

"""

import csv
import re
from collections.abc import Generator
from functools import cached_property
from logging import Logger, getLogger
from pathlib import Path
from typing import Annotated, Any, Final, Self

import dlt
import dlt.sources.helpers.requests
from dlt.common.pipeline import LoadInfo
from dlt.extract.items import DataItemWithMeta
from dlt.sources.helpers.rest_client.client import RESTClient
from pydantic import Field, StringConstraints, computed_field, model_validator
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.settings import (
    CLI_SHORTCUTS,
    CtsSettings,
    default_settings_with_shortcuts,
)
from cdm_data_loaders.pipelines.core import (
    PipelineCli,
    run_pipeline,
)
from cdm_data_loaders.utils.file_transfer.sync_client import FileDownloader

logger: Logger = getLogger(__name__)

NonEmptyStr = Annotated[str, StringConstraints(min_length=1)]


DATASET_NAME = "all_the_bacteria"
ALL_FILES_TSV_FILE_ID = "R6gcp"
ALL_ATB_FILE_NAME = "all_atb_files.tsv"
REGEX_FILE = "filters.txt"
ATB_VERSION = "2025-05"

VERSION: Final[str] = "version"
PATTERN_FILE: Final[str] = "pattern_file"

# project parts needed:
PROJECT_PARTS = ["Annotation/Bakta", "Assembly", "Metadata"]
PROJECT_PART_REGEX = re.compile(f"^AllTheBacteria/({'|'.join(PROJECT_PARTS)})")

EXPECTED_ATB_FIELDNAMES = ["project", "project_id", "filename", "url", "md5", "size(MB)"]
REQUIRED_ATB_FIELDNAMES = {"project", "filename", "url", "md5"}


class AtbSettings(CtsSettings):
    """Configuration for running the AllTheBacteria import pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name="all_the_bacteria",
        cli_shortcuts={VERSION: "v", **CLI_SHORTCUTS},
    )

    version: NonEmptyStr = Field(
        default=ATB_VERSION,
        description="Name of the current AllTheBacteria version",
    )

    pattern_file: NonEmptyStr | None = Field(
        default=None,
        description="Path, relative to the input dir, of a file containing patterns to match when downloading ATB files",
    )

    @model_validator(mode="after")
    def post_init_checks(self) -> Self:
        """Ensure that the raw data dir can be created."""
        # check that the patterns are OK
        _ = self.pattern_matches

        return self

    @computed_field
    @cached_property
    def raw_data_dir(self) -> Path:
        """Directory to save downloaded raw data files to: ``<output_dir>/raw_data/<version>``."""
        return self.output_dir / "raw_data" / self.version

    @computed_field
    @cached_property
    def pattern_matches(self) -> re.Pattern:
        """The regular expression pattern to be used to select files for download.

        If a pattern_file is supplied, it will read in the file at {input_dir}/{pattern_file}
        and convert the contents into a regular expression. If no file is supplied or the file is empty
        or does not contain any content, the default PROJECT_PART_REGEX will be used instead.
        """
        if self.pattern_file:
            pattern_file = self.input_dir / self.pattern_file
            regex = load_patterns(pattern_file)
            if regex is None:
                err_msg = f"No patterns found in {pattern_file!s}"
                raise ValueError(err_msg)
            # return the compiled regexes from the file
            return regex
        # return the default
        return PROJECT_PART_REGEX


def load_patterns(pattern_file: Path) -> re.Pattern | None:
    """Load the pattern file and convert it into a set of regexes."""
    patterns = []
    try:
        for line in pattern_file.read_text(encoding="utf-8").splitlines():
            trimmed_line = line.strip()
            # skip blank lines
            if not trimmed_line:
                continue
            patterns.append(
                re.escape(trimmed_line[:-1]) + ".*" if trimmed_line.endswith("*") else re.escape(trimmed_line)
            )

        if patterns:
            return re.compile("^(" + "|".join(patterns) + ")$")
    except Exception:
        logger.exception("Could not load patterns from %s", str(pattern_file))

    return None


def download_atb_index_tsv(settings: AtbSettings) -> Path:
    """Download the ATB file index TSV file from the OSF and save it to disk.

    :param settings: pipeline config
    :type settings: AtbSettings
    :raises RuntimeError: if the download URL cannot be found
    :return: path to the downloaded file
    :rtype: Path
    """
    # make sure that the directory structure to save the file in can be written to
    settings.raw_data_dir.mkdir(parents=True, exist_ok=True)

    # get the all_atb_files.tsv file info from the OSF API and retrieve the download link
    osf_client = RESTClient(
        base_url="https://api.osf.io/v2/",
        headers={"accept": "application/json"},
        data_selector="data",
    )
    resp = osf_client.get(f"https://api.osf.io/v2/files/{ALL_FILES_TSV_FILE_ID}/")
    resp.raise_for_status()

    resp_json = resp.json()
    all_files_tsv_download = resp_json.get("data", {}).get("links", {}).get("download")
    if all_files_tsv_download is None:
        logger.error("Could not find download URL in OSF API response:")
        logger.error(resp_json)
        err_msg = f"Could not find download URL in response from 'https://api.osf.io/v2/files/{ALL_FILES_TSV_FILE_ID}/'"
        raise RuntimeError(err_msg)

    atb_files_tsv = settings.raw_data_dir / ALL_ATB_FILE_NAME

    # download the file listing and save it
    # TODO: for non-local destinations, the file should also be copied to the remote
    FileDownloader().download(url=all_files_tsv_download, destination=atb_files_tsv)
    logger.info("Downloaded TSV index file.")
    return atb_files_tsv


def get_file_download_links(settings: AtbSettings, atb_files_tsv: Path) -> Generator[list[dict[str, Any]], Any]:
    """Parse the ATB file index TSV and to yield a list of files to download.

    :param settings: pipeline config
    :type settings: AtbSettings
    :param atb_files_tsv: path to the ATB file index TSV file
    :type atb_files_tsv: Path
    :yield: list of fields to download
    :rtype: Generator[list[dict[str, Any]], Any]
    """
    pattern_to_match = settings.pattern_matches
    with atb_files_tsv.open() as index_file:
        reader = csv.DictReader(index_file, delimiter="\t")
        all_lines = list(reader)
        if not all_lines:
            err_msg = f"No valid TSV data found in {atb_files_tsv!s}"
            logger.error(err_msg)
            raise RuntimeError(err_msg)

        if reader.fieldnames != EXPECTED_ATB_FIELDNAMES:
            err_msg = f"ATB file index TSV headers have changed.\nExpected: {EXPECTED_ATB_FIELDNAMES}\nGot: {reader.fieldnames}"
            logger.warning(err_msg)
            # do we have the essentials? project, filename, url, md5
            missing_required_fields = [f for f in REQUIRED_ATB_FIELDNAMES if f not in (reader.fieldnames or [])]
            if missing_required_fields:
                err_msg = f"Missing required ATB file index TSV headers: {sorted(missing_required_fields)}"
                logger.error(err_msg)
                raise RuntimeError(err_msg)

        yield [row for row in all_lines if pattern_to_match.match(row["project"])]


def osf_file_downloader(settings: AtbSettings, atb_file_list: list[dict[str, Any]]) -> Generator[DataItemWithMeta, Any]:
    """Download files from OSF to the local storage space.

    :param settings: pipeline config
    :type settings: Settings
    :param atb_file_list: info about files to transfer, as a list of dictionaries
    :type atb_file_list: list[dict[str, Any]]
    """
    client = FileDownloader()
    successful_downloads = []
    for f in atb_file_list:
        try:
            project_part = f["project"].removeprefix("AllTheBacteria").removeprefix("/").rstrip("/")
            # TODO: stream to non-local destinations
            save_path = settings.raw_data_dir / project_part / f["filename"]
            client.download(url=f["url"], destination=save_path, expected_checksum=f["md5"], checksum_fn="md5")

            f["path"] = str(save_path)
            successful_downloads.append(f)
        except Exception as e:
            # do something!
            err_msg = f"Could not download file from {f['url']}: {e!s}"
            logger.exception(err_msg)
            continue

    yield dlt.mark.with_table_name(successful_downloads, "downloaded_files")


@dlt.resource(name="atb_file_list")
def atb_file_list(settings: AtbSettings) -> Generator[list[dict[str, Any]], Any, Any]:
    """Generate a list of files to download from the list of all ATB files."""
    atb_files_tsv = download_atb_index_tsv(settings)
    return get_file_download_links(settings, atb_files_tsv)


@dlt.transformer(name="file_downloader", data_from=atb_file_list, parallelized=True)
def file_downloader(
    atb_file_list: list[dict[str, Any]],
    settings: AtbSettings,
) -> Generator[DataItemWithMeta, Any]:
    """Download ATB files to disk.

    :param settings: pipeline config
    :type settings: Settings
    :param atb_file_list: list of files to download
    :type atb_file_list: list[dict[str, Any]]
    :return: output of the osf_file_downloader
    :rtype: Generator[DataItemWithMeta]
    """
    return osf_file_downloader(settings, atb_file_list)


def run_atb_pipeline(settings: AtbSettings) -> LoadInfo | None:
    """Run the AllTheBacteria pipeline.

    :param settings: configuration for the pipeline
    :type settings: AtbSettings
    """
    atb_file_list.bind(settings)
    file_downloader.bind(settings)

    pipeline_kwargs: dict[str, Any] = {
        "pipeline_name": DATASET_NAME,
        "dataset_name": DATASET_NAME,
    }

    return run_pipeline(
        settings=settings,
        resource=file_downloader,
        destination_kwargs={"max_table_nesting": 0},
        pipeline_kwargs=pipeline_kwargs,
        pipeline_run_kwargs=None,
    )


cli = PipelineCli(AtbSettings, run_atb_pipeline)

if __name__ == "__main__":
    cli()
