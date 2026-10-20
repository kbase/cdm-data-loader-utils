"""Pipeline to import data from the NCBI API.

Like other importers in the cdm-data-loaders repo, the NCBI REST API pipeline
expects input files with the naming format `{text_string}_{digits}.{ext}`, e.g.
batch_001.txt, batch_002.txt, batch_003.txt.
"""

import os
import re
from collections.abc import Generator
from functools import cached_property, partial
from itertools import islice
from logging import Logger, getLogger
from pathlib import Path
from typing import Any, Final
from urllib.parse import parse_qs, urlparse

import dlt
from dlt.common.pipeline import LoadInfo
from dlt.extract.items import DataItemWithMeta
from dlt.sources.helpers.requests import Response
from dlt.sources.helpers.rest_client.auth import APIKeyAuth
from dlt.sources.helpers.rest_client.client import RESTClient
from dlt.sources.helpers.rest_client.paginators import (
    JSONResponseCursorPaginator,
)
from pydantic import Field, computed_field, field_validator
from pydantic_settings import SettingsConfigDict
from requests.exceptions import HTTPError

from cdm_data_loaders.core.fields import SaveRawResponses
from cdm_data_loaders.core.settings import CLI_SHORTCUTS, CtsSettings, default_settings_with_shortcuts
from cdm_data_loaders.pipelines.core import (
    PipelineCli,
    run_pipeline,
)
from cdm_data_loaders.utils.batcher import NumericFileSequenceBatcher

logger: Logger = getLogger(__name__)

DATASET_NAME: Final[str] = "ncbi_rest_api"

NCBI_API_KEY = os.environ.get("NCBI_API_KEY") or "DEMO_KEY"

# Max number of items to request per page from the NCBI REST API (max allowed is 1000).
MAX_RESULTS_PER_PAGE: Final[int] = 1000
# Max number of IDs to send in a multi-ID query
# Max URL length seems to be 4611 (as of Jun 2026)
MAX_IDS_PER_QUERY: Final[int] = 250

DATASET: Final[str] = "dataset"
ANNOTATION: Final[str] = "annotation"
ERROR: Final[str] = "error"

QUERY_TYPE_REGEX: Final[re.Pattern[str]] = re.compile(r"^(" + DATASET + r"|" + ANNOTATION + r")$")

REST_CLIENT_HOOKS = {}

BATCH_SIZE: Final[str] = "batch_size"
QUERY_TYPE: Final[str] = "query_type"


class NcbiRestApiSettings(CtsSettings):
    """Configuration for running the NCBI REST API import pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name="ncbi_rest_api",
        cli_shortcuts={**CLI_SHORTCUTS, BATCH_SIZE: "b", QUERY_TYPE: "q"},
    )

    batch_size: int = Field(
        default=MAX_IDS_PER_QUERY,
        description="Number of IDs to send in each request to the NCBI REST API.",
        ge=1,
        le=MAX_IDS_PER_QUERY,
    )

    query_type: str | None = Field(
        default=None,
        description="The type of query to perform, dataset or annotation. By default, both are performed.",
        pattern=QUERY_TYPE_REGEX,
    )

    save_raw_responses: SaveRawResponses

    @computed_field
    @cached_property
    def raw_data_dir(self) -> Path:
        """Directory to save downloaded raw data files to: ``<output_dir>/raw_data/<version>``."""
        return self.output_dir / "raw_data" / "api_responses"

    @field_validator(QUERY_TYPE, mode="before")
    @classmethod
    def trim_lc_query_type(cls, v: str | None) -> str | None:
        """Prepare the query_type parameter for validation.

        :param v: query_type parameter
        :type v: str | None
        :return: cleaned query type param
        :rtype: str | None
        """
        return v.strip().lower() if v and v.strip() else None


def generate_file_path_name_from_url(settings: NcbiRestApiSettings, url_string: str) -> Path:
    """Given an NCBI URL, generate a save directory and file name for the raw HTTP responses.

    :param settings: settings object
    :type  settings: NcbiRestApiSettings
    :param url_string: the URL
    :type url_string: str
    :return: path for the save dir
    :rtype: Path
    """
    url = urlparse(url_string)
    # retrieve the page token
    page_token = "".join(parse_qs(url.query).get("page_token", []))
    # save_parts 0 is the broad type, e.g. gene, protein, genome
    # save_parts 1 is the type of PID
    # save_parts 2 is the value of the PID
    # any subsequent parts are the specific data type being accessed
    save_parts = url.path.replace("/datasets/v2/", "").split("/")
    # make a directory for this ID
    save_dir = settings.raw_data_dir / save_parts[3] / save_parts[2]
    save_dir.mkdir(parents=True, exist_ok=True)
    # file name is the PID plus optional page token if this was a multi-request PID
    file_name = f"{save_parts[2]}__{page_token}.json" if page_token else f"{save_parts[2]}.json"
    return save_dir / file_name


def save_raw_response(settings: NcbiRestApiSettings, response: Response, *_: Any, **__: Any) -> Response:  # noqa: ANN401
    """Save the response content to disk.

    :param response: HTTP Response object
    :type response: Response
    :return: the response
    :rtype: Response
    """
    file_path = generate_file_path_name_from_url(settings, response.url)
    file_path.write_bytes(response.content)
    return response


PIPELINE_SETTINGS: NcbiRestApiSettings | None = None


def set_settings(settings_obj: NcbiRestApiSettings) -> None:
    """Set the global PIPELINE_SETTINGS object to the supplied obj.

    :param settings_obj: settings object
    :type settings_obj: NcbiRestApiSettings
    """
    # ugh, somewhat horrible
    global PIPELINE_SETTINGS
    PIPELINE_SETTINGS = settings_obj


def get_settings() -> NcbiRestApiSettings:
    """Get the pipeline settings.

    return: pipeline settings
    :rtype: NcbiRestApiSettings
    """
    if PIPELINE_SETTINGS is None:
        err_msg = "Pipeline settings have not been initialised"
        raise RuntimeError(err_msg)
    return PIPELINE_SETTINGS


ncbi_genome_client = RESTClient(
    base_url="https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/",
    headers={"accept": "application/json"},
    auth=APIKeyAuth(name="api_key", api_key=NCBI_API_KEY),
    paginator=JSONResponseCursorPaginator(cursor_path="next_page_token", cursor_param="page_token"),
    data_selector="reports",
)


def add_error(
    error_list: list[dict[str, Any]],
    error: Exception,
    error_from: str,
    assembly_id: str | None = None,
    assembly_id_list: list[str] | None = None,
) -> None:
    """Add an error to the list of output errors.

    :param error_list: running list of errors
    :type error_list: list[dict[str, Any]]
    :param error: the error object from the exception handler
    :type error: Exception
    :param error_from: what type of request was being made when the error occurred
    :type error_from: str
    :param assembly_id: ID of the assembly being fetched when the error occurred, defaults to None
    :type assembly_id: str | None, optional
    :param assembly_id_list: list of IDs being fetched when the error occurred, defaults to None
    :type assembly_id_list: list[str] | None, optional
    """
    err_args = {
        "assembly_id": assembly_id or None,
        "assembly_id_list": assembly_id_list or None,
        "error_class": type(error).__name__,
        "error_from": error_from,
        "message": str(error),
        "request_url": None,
        "status": None,
        "reason": None,
    }

    if isinstance(error, HTTPError):
        # save the URL, status code, and error message
        err_args = {
            **err_args,
            "request_url": error.request.url,
            "status": error.response.status_code,
            "reason": error.response.reason,
        }

    error_list.append(err_args)


def get_assembly_reports(assembly_id_list: list[str]) -> dict[str, Any]:
    """Retrieve dataset and annotation reports for a list of IDs from the NCBI datasets API.

    :param assembly_id_list: list of IDs to retrieve data for
    :type assembly_id_list: list[str]
    :return: dictionary with keys DATASET and ANNOTATION containing dataset and annotation data respectively.
    :rtype: dict[str, Any]
    """
    if not assembly_id_list:
        return {}

    settings = get_settings()
    fetch_dataset: bool = settings.query_type in (None, DATASET)
    fetch_annotation: bool = settings.query_type in (None, ANNOTATION)
    errors: list[dict[str, Any]] = []

    # N.b. invalid IDs will not be present in dataset_reports
    dataset_reports = {}
    if fetch_dataset:
        try:
            dataset_reports = get_dataset_reports(assembly_id_list)
        except Exception as e:  # noqa: BLE001
            add_error(errors, e, "dataset_report", assembly_id_list=assembly_id_list)

    annotation_reports: dict[str, Any] = {}
    if fetch_annotation:
        for assembly_id in assembly_id_list:
            try:
                annotation_reports[assembly_id] = get_annotation_report(assembly_id)
            except Exception as e:  # noqa: BLE001
                add_error(errors, e, "annotation_report", assembly_id=assembly_id)

    # ensure every assembly_id in the list has either the downloaded dataset_report or None
    output: dict[str, Any] = {
        ERROR: errors,
    }
    if fetch_annotation:
        output[ANNOTATION] = {assembly_id: annotation_reports.get(assembly_id) for assembly_id in assembly_id_list}

    if fetch_dataset:
        output[DATASET] = {assembly_id: dataset_reports.get(assembly_id) for assembly_id in assembly_id_list}

    return output


def get_dataset_reports(assembly_id_list: list[str]) -> dict[str, dict[str, Any] | None]:
    """Fetch the dataset report for a list of assemblies from the NCBI datasets REST API."""
    if not assembly_id_list:
        return {}

    logger.info("fetching dataset reports for:\n%s", ", ".join(sorted(assembly_id_list)))
    assembly_dataset_reports = []

    for page in ncbi_genome_client.paginate(
        f"{'%2C'.join(assembly_id_list)}/dataset_report",
        params={
            "page_size": MAX_RESULTS_PER_PAGE,
        },
        hooks=REST_CLIENT_HOOKS,  # type: ignore[reportArgumentType]
    ):
        assembly_dataset_reports.extend(page)

    # return dataset reports, indexed by assembly_id
    # invalid IDs are silently dropped by the NCBI REST API
    datasets = {report.get("accession"): report for report in assembly_dataset_reports}
    # fill in the missing gaps in assembly_id_list with None
    return {assembly_id: datasets.get(assembly_id) for assembly_id in assembly_id_list}


def get_annotation_report(assembly_id: str) -> list[dict[str, Any]] | None:
    """Fetch the annotation report for an assembly from the NCBI datasets REST API."""
    logger.info("fetching annotation report for %s", assembly_id)
    page_data = []

    for page in ncbi_genome_client.paginate(
        f"{assembly_id}/annotation_report",
        params={
            "page_size": MAX_RESULTS_PER_PAGE,
        },
        hooks=REST_CLIENT_HOOKS,  # type: ignore[reportArgumentType]
    ):
        page_data.extend(page)

    # page_data is empty if the ID is invalid
    return page_data or None


def assemble_assembly_reports(
    assembly_reports: dict[str, dict[str, Any]],
) -> Generator[DataItemWithMeta, Any]:
    """Parse and export assembly data from the REST API.

    :param assembly_reports: output from the dataset_report and annotation_report API endpoint, indexed by assembly_id
    :type assembly_reports: dict[str, dict[str, None | dict[str, Any]]]
    :yield: table rows for assemblies
    :rtype: Generator[DataItemWithMeta, Any]
    """
    if not assembly_reports:
        return

    dataset_reports: dict[str, dict[str, Any]] = assembly_reports.get(DATASET)  # type: ignore[reportAssignmentType]
    annotation_reports: dict[str, list[dict[str, Any]]] = assembly_reports.get(ANNOTATION)  # type: ignore[reportAssignmentType]
    error_reports: list[dict[str, Any]] = assembly_reports.get(ERROR)  # type: ignore[reportAssignmentType]

    # yield the raw data to save as tables
    # every assembly_id in dataset_reports should have a row in the dataset_report table
    # even if the value is None due to an error or invalid ID
    if dataset_reports:
        yield dlt.mark.with_table_name(
            [{"assembly_id": assembly_id, **(report or {})} for assembly_id, report in dataset_reports.items()],
            f"{DATASET}_report",
        )

    # save the rows from the annotation reports
    # only assemblies with annotation reports have rows in the annotation_report table
    tr_annotation_rows = [
        {"assembly_id": aid, **report}
        for aid, report_list in annotation_reports.items()
        for report in (report_list or [])
    ]
    if tr_annotation_rows:
        yield dlt.mark.with_table_name(tr_annotation_rows, f"{ANNOTATION}_report")

    if error_reports:
        yield dlt.mark.with_table_name(error_reports, "ncbi_import_error")


@dlt.resource(name="assembly_list")
def assembly_list() -> Generator[list[str], Any, Any]:
    """List of assemblies to fetch."""
    settings = get_settings()
    batcher = NumericFileSequenceBatcher(directory=settings.input_dir, batch_size=1)
    while files := batcher.get_batch():
        for file_path in files:
            with file_path.open() as f:
                while next_chunk := list(islice(f, settings.batch_size)):
                    yield [line.strip() for line in next_chunk if line.strip()]


@dlt.transformer(name="assembly_reports", data_from=assembly_list, max_table_nesting=0, parallelized=True)
def assembly_reports(
    assembly_id_list: list[str],
) -> Generator[dict[str, dict[str, dict[str, Any]] | list[dict[str, Any]]], Any]:
    """Fetch reports for a list of assemblies from the NCBI datasets REST API."""
    yield get_assembly_reports(assembly_id_list)


@dlt.transformer(name="assembly_report_parser", data_from=assembly_reports, max_table_nesting=0, parallelized=True)
def assembly_report_parser(
    assembly_reports: dict[str, dict[str, Any]],
) -> Generator[DataItemWithMeta, Any]:
    """Parse and export assembly data from the REST API."""
    return assemble_assembly_reports(assembly_reports)


def run_ncbi_pipeline(settings: NcbiRestApiSettings) -> LoadInfo | None:
    """Run the NCBI datasets API pipeline.

    :param settings: configuration for the pipeline
    :type settings: NcbiRestApiSettings
    """
    set_settings(settings)

    if settings.save_raw_responses:
        REST_CLIENT_HOOKS["response"] = [partial(save_raw_response, settings)]

    pipeline_kwargs = {
        "pipeline_name": DATASET_NAME,
        "dataset_name": DATASET_NAME,
    }

    return run_pipeline(
        settings=settings,
        resource=[assembly_report_parser],
        destination_kwargs={"max_table_nesting": 0},
        pipeline_kwargs=pipeline_kwargs,
        pipeline_run_kwargs=None,
    )


cli = PipelineCli(NcbiRestApiSettings, run_ncbi_pipeline)


if __name__ == "__main__":
    cli()
