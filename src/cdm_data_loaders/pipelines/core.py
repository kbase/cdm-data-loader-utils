"""Common reusable pipeline elements."""

from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime
from logging import Logger, getLogger
from pathlib import Path
from shutil import move
from typing import Any, Final

import dlt
from dlt.common.pipeline import LoadInfo
from dlt.common.runtime.slack import send_slack_message
from dlt.extract import DltResource
from dlt.pipeline.exceptions import PipelineStepFailed
from pydantic import ValidationError
from pydantic_settings import CliApp, SettingsError

from cdm_data_loaders.core.fields import JSONL
from cdm_data_loaders.core.settings import CtsSettings, LoggerSettings
from cdm_data_loaders.utils.cdm_logger import init_logger

WEBHOOK_NOT_CONFIGURED: Final[str] = "Slack webhook not configured"
NO_MESSAGE: Final[str] = "No message supplied"
LOAD_INFO_TABLE_NAME: Final[str] = "_load_info"
DISABLE_COMPRESSION_ENV_VAR: Final[str] = "NORMALIZE__DATA_WRITER__DISABLE_COMPRESSION"
DESTINATION_NOT_CONFIGURED: Final[str] = (
    "Destination {name!r} is not configured: set destination.{name}.destination_type in the dlt config."
)

logger: Logger = getLogger(__name__)


@dataclass(frozen=True, slots=True)
class PipelineCli[S: LoggerSettings]:
    """Console-script entry point for a pipeline; ``run`` returns the LoadInfo."""

    settings_cls: type[S]
    pipeline_fn: Callable[[S], LoadInfo | None]

    def __call__(self) -> None:
        """Parse sys.argv and run the pipeline. Returns None so the process exits with 0."""
        self.run()

    def run(self, cli_args: list[str] | None = None, **settings_kwargs: Any) -> LoadInfo | None:
        """Parse ``cli_args`` (default ``sys.argv[1:]``), run the pipeline and return its LoadInfo."""
        return run_cli(self.settings_cls, self.pipeline_fn, settings_kwargs, cli_args)


@dataclass(frozen=True, slots=True)
class PipelineReport:
    """Summary information for a pipeline run report."""

    pipeline_name: str
    dataset_name: str
    started_at: datetime
    finished_at: datetime
    loads_ids: list[str] | None = None
    failed_messages: list[str] | None = None
    err_msg: str | None = None
    info_load_err_msg: str | None = None


def send_slack_message_carefully(
    slack_hook: str,
    message: str,
    is_markdown: bool = False,  # noqa: FBT001, FBT002
) -> None:
    """Carefully send a slack message by wrapping it in a try/except.

    :param slack_hook: slack webhook URL
    :type slack_hook: str
    :param message: message to be slacked
    :type message: str
    :param is_markdown: whether or not the message is markdown, defaults to False
    :type is_markdown: bool, optional
    """
    if not slack_hook:
        logger.warning("Cannot send slack message: %s", WEBHOOK_NOT_CONFIGURED)
        return
    if not message or not message.strip():
        logger.warning("Cannot send slack message: %s", NO_MESSAGE)
        return

    try:
        send_slack_message(slack_hook, message.strip(), is_markdown)
    except Exception:
        logger.exception("Failed to send slack message")


def run_cli(
    settings_cls: type[LoggerSettings],
    pipeline_fn: Callable[[Any], LoadInfo | None],
    settings_kwargs: dict[str, Any] | None = None,
    cli_args: list[str] | None = None,
) -> LoadInfo | None:
    """Run a DLT pipeline from a generic CLI entry point.

    :param settings_cls: the Settings class to instantiate
    :type  settings_cls: type[LoggerSettings]
    :param pipeline_fn: the run_pipeline function to call with the config
    :type  pipeline_fn: Callable[[Any], LoadInfo | None]
    :param settings_kwargs: any extra non-cli/env var settings to be added
    :type  settings_kwargs: dict[str, Any] | None, optional
    :param cli_args: command line args
    :type  cli_args: list[str] | None, optional
    :return: pipeline load information
    :rtype: LoadInfo | None
    """
    try:
        settings = CliApp.run(settings_cls, cli_args=cli_args, **(settings_kwargs or {}))
        init_logger(settings)
    except (SettingsError, ValidationError, ValueError):
        logger.exception("Error initialising config")
        raise
    except Exception:
        logger.exception("Unexpected error setting up config")
        raise
    return pipeline_fn(settings)


def write_pipeline_report(
    report_path: Path,
    report: PipelineReport,
) -> None:
    """Write a summary report of a pipeline run to ``report_path``.

    :param report_path: path to the report file to write
    :type  report_path: Path
    :param report: pipeline report information to write
    :type  report: PipelineReport
    """
    lines = [
        f"Pipeline: {report.pipeline_name}",
        f"Dataset: {report.dataset_name}",
        f"Started: {report.started_at}",
        f"Finished: {report.finished_at}",
    ]
    if report.loads_ids is None:
        lines.append("Status: error")
        if report.err_msg:
            lines.append(report.err_msg)
    else:
        failed_messages = [msg for msg in (report.failed_messages or []) if msg]
        status = "completed with errors" if failed_messages else "completed"
        lines.append(f"Status: {status}")
        lines.append(f"Completed loads: {', '.join(report.loads_ids)}")
        lines.extend(failed_messages)
        if report.info_load_err_msg is None:
            lines.append("Pipeline load info saved")
        else:
            lines.append(f"Pipeline load info not saved: {report.info_load_err_msg!s}")

    try:
        report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    except Exception as err:
        report_err_msg = f"Could not write pipeline report file: {err!s}"
        logger.exception(report_err_msg)


def save_pipeline_metadata(
    settings: CtsSettings, pipeline: dlt.Pipeline, load_info: LoadInfo | None, err_msg: str | None
) -> None:
    """Save the pipeline metadata to a unique directory in settings.output_dir.

    Writes a brief report of the pipeline state and saves it to report_<unique pipeline load id>.txt.

    :param settings: settings for this pipeline
    :type  settings: CtsSettings
    :param pipeline: dlt pipeline instance
    :type  pipeline: dlt.Pipeline
    :param load_info: pipeline load info, or None if the pipeline did not complete successfully
    :type  load_info: LoadInfo | None
    :param err_msg: error message indicating why pipeline failed, or None if the pipeline was successful
    :type  err_msg: str | None
    """
    started_at = pipeline.last_trace.started_at
    base_id = load_info.loads_ids[0] if load_info and load_info.loads_ids else str(started_at.timestamp())
    finished_at = pipeline.last_trace.finished_at
    pipeline_e2e_id = f"{base_id}__{finished_at}"

    info_load_err_msg: str | None = None
    if load_info is not None:
        try:
            # save the load info as JSONL in the same destination as the ingested data
            load_info_resource: DltResource = dlt.resource(
                [load_info.asdict()],
                name=LOAD_INFO_TABLE_NAME,
                max_table_nesting=0,
                columns={
                    k: {"data_type": "text"}
                    for k in [
                        "environment",
                        "staging_type",
                        "staging_name",
                        "staging_display_credentials",
                        "dataset_name",
                    ]
                },
            )
            load_info_save_info: LoadInfo = pipeline.run(load_info_resource, loader_file_format=str(JSONL))
            if load_info_save_info.has_failed_jobs:
                info_load_err_msg = f"Loading pipeline load info had failed jobs: {load_info_save_info}"
                logger.error(info_load_err_msg)
        except Exception as e:
            info_load_err_msg = f"Loading pipeline load info failed: {e!s}"
            logger.exception(info_load_err_msg)

    if settings.use_output_dir_for_pipeline_metadata and settings.pipeline_dir is not None:
        # copy the pipeline directory so that it is named according to the first load_id
        try:
            move(settings.pipeline_dir, settings.output_dir / pipeline_e2e_id)
        except OSError as e:
            move_err_msg = f"Could not move pipelines dir: {e!s}"
            logger.exception(move_err_msg)

    # write a summary and save it to the output directory
    if load_info is None:
        report = PipelineReport(
            pipeline_name=pipeline.pipeline_name,
            dataset_name=pipeline.dataset_name,
            started_at=started_at,
            finished_at=finished_at,
            err_msg=err_msg,
        )
    else:
        report = PipelineReport(
            pipeline_name=pipeline.pipeline_name,
            dataset_name=pipeline.dataset_name,
            started_at=started_at,
            finished_at=finished_at,
            err_msg=err_msg,
            loads_ids=load_info.loads_ids,
            failed_messages=[
                failed_job.failed_message
                for load_package in load_info.load_packages
                for failed_job in load_package.jobs["failed_jobs"]
                if failed_job.failed_message
            ],
            info_load_err_msg=info_load_err_msg,
        )
    write_pipeline_report(settings.output_dir / f"report_{pipeline_e2e_id}.txt", report)


def run_pipeline(
    settings: CtsSettings,
    resource: DltResource | list[DltResource],
    destination_kwargs: dict[str, Any] | None = None,
    pipeline_kwargs: dict[str, Any] | None = None,
    pipeline_run_kwargs: dict[str, Any] | None = None,
) -> LoadInfo | None:
    """Execute a dlt pipeline.

    :param settings: pipeline settings
    :type settings: CtsSettings
    :param resource: dlt resource to run
    :type resource: DltResource | list[DltResource]
    :param destination_kwargs: keyword arguments for the dlt destination
    :type destination_kwargs: dict[str, Any] | None
    :param pipeline_kwargs: keyword arguments for the dlt pipeline
    :type pipeline_kwargs: dict[str, Any] | None
    :param pipeline_run_kwargs: keyword arguments for the dlt pipeline run
    :type pipeline_run_kwargs: dict[str, Any] | None
    :raises ValueError: if settings.use_destination is not configured in the dlt config
    :return: a load info object, if successful; otherwise None
    :rtype: LoadInfo | None
    """
    pipeline_kwargs = dict(pipeline_kwargs or {})

    # set the output directory for all the pipeline gubbins
    if settings.use_output_dir_for_pipeline_metadata:
        pipeline_kwargs["pipelines_dir"] = settings.pipeline_dir

    # update dev mode
    if settings.dlt_dev_mode:
        pipeline_kwargs["dev_mode"] = settings.dlt_dev_mode

    pipeline = dlt.pipeline(**pipeline_kwargs)

    try:
        pipeline.destination = dlt.destination(settings.use_destination, **(destination_kwargs or {}))
    except Exception as e:
        err_msg = f"Invalid destination configuration for {settings.use_destination}"
        logger.exception(err_msg)
        raise ValueError(err_msg) from e

    slack_hook: str | None = pipeline.runtime_config.slack_incoming_hook

    if not slack_hook:
        logger.info("No Slack alerts will be sent: %s", WEBHOOK_NOT_CONFIGURED)

    load_info: LoadInfo | None = None
    err_msg: str | None = None
    try:
        load_info = pipeline.run(resource, **(pipeline_run_kwargs or {}))
    except PipelineStepFailed as e:
        # errors raise PipelineStepFailed
        err_msg = f"Pipeline failed: {e!s}"
        logger.exception(err_msg)

    if load_info:
        logger.info(load_info)
        if load_info.has_failed_jobs:
            logger.warning("Pipeline contains failed jobs")
        else:
            logger.info("All pipeline jobs completed successfully")

    save_pipeline_metadata(settings, pipeline, load_info, err_msg)

    if slack_hook:
        msg = err_msg or f"Pipeline completed {'with errors' if load_info.has_failed_jobs else 'successfully!'}"
        send_slack_message_carefully(slack_hook, msg)

    return load_info
