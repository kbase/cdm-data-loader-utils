"""Common reusable pipeline elements."""

from collections.abc import Callable
from dataclasses import dataclass
from logging import Logger, getLogger
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


def send_slack_message_carefully(slack_hook: str, message: str, is_markdown: bool = False) -> None:  # noqa: FBT001, FBT002
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
    if settings.pipeline_dir:
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

    try:
        load_info: LoadInfo = pipeline.run(resource, **(pipeline_run_kwargs or {}))
    except PipelineStepFailed as e:
        # errors raise PipelineStepFailed
        err_msg = f"Pipeline failed: {e!s}"
        logger.exception(err_msg)
        if slack_hook:
            send_slack_message_carefully(slack_hook, err_msg)
        return None

    logger.info(load_info)
    if load_info.has_failed_jobs:
        logger.warning("")

    try:
        # save the load info as JSONL in the same destination as the ingested data
        load_info_resource: DltResource = dlt.resource(
            [load_info.asdict()], name=LOAD_INFO_TABLE_NAME, max_table_nesting=0
        )
        load_info_save_info: LoadInfo = pipeline.run(load_info_resource, loader_file_format=JSONL)
        if load_info_save_info.has_failed_jobs:
            logger.error("Loading pipeline load info had failed jobs: %s", load_info_save_info)
    except Exception as e:
        err_msg = f"Loading pipeline load info failed: {e!s}"
        logger.exception(err_msg)

    logger.info("Work complete!")
    if slack_hook:
        send_slack_message_carefully(
            slack_hook, f"Pipeline completed {'with errors' if load_info.has_failed_jobs else 'successfully!'}"
        )

    return load_info
