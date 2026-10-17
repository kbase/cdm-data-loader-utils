"""Common reusable pipeline elements."""

from collections.abc import Callable
from dataclasses import dataclass
from logging import Logger, getLogger
from typing import Any, Final, TypeVar

import dlt
from dlt.common.destination.reference import TDestinationReferenceArg
from dlt.common.pipeline import LoadInfo
from dlt.common.runtime.slack import send_slack_message
from dlt.extract import DltResource
from dlt.pipeline.exceptions import PipelineStepFailed
from dlt.sources.filesystem import filesystem
from pydantic import ValidationError
from pydantic_settings import CliApp, SettingsError

from cdm_data_loaders.core.destination import ConfigLookup, destinations_from_dlt_config, resolve_output_dir
from cdm_data_loaders.core.fields import JSONL, OUTPUT_DIR
from cdm_data_loaders.core.settings import CtsSettings, LoggerSettings
from cdm_data_loaders.utils.cdm_logger import init_logger

WEBHOOK_NOT_CONFIGURED: Final[str] = "Slack webhook not configured"
NO_MESSAGE: Final[str] = "No message supplied"
LOAD_INFO_TABLE_NAME: Final[str] = "_load_info"
DISABLE_COMPRESSION_ENV_VAR: Final[str] = "NORMALIZE__DATA_WRITER__DISABLE_COMPRESSION"
UNRESOLVED_OUTPUT_DIR: Final[str] = (
    "settings.output_dir is not set. Pass settings through resolve_cts_settings() "
    "(run_cli does this automatically), or construct the settings with an explicit output_dir."
)

CtsSettingsT = TypeVar("CtsSettingsT", bound=CtsSettings)

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


def resolve_cts_settings[CtsSettingsT: CtsSettings](
    settings: CtsSettingsT, dlt_config: ConfigLookup | None = None
) -> CtsSettingsT:
    """Return a copy of ``settings`` with output_dir resolved against the dlt destination config.

    An explicit output_dir takes precedence over ``destination.<use_destination>.bucket_url``.
    Neither ``settings`` nor the dlt config is modified.

    :param settings: settings parsed from the CLI and env vars
    :type settings: CtsSettingsT
    :param dlt_config: config to resolve against; defaults to ``dlt.config``. Tests can pass a nested dict.
    :type dlt_config: ConfigLookup | None
    :raises ValueError: if the output location cannot be resolved or is inconsistent with the config
    :return: a copy of the settings with output_dir set
    :rtype: CtsSettingsT
    """
    config = dlt.config if dlt_config is None else dlt_config
    output_dir = resolve_output_dir(
        use_destination=settings.use_destination,
        output_dir=settings.output_dir,
        destinations=destinations_from_dlt_config(config),
        pipeline_metadata_in_output_dir=settings.use_output_dir_for_pipeline_metadata,
    )
    # resolve_output_dir returns a normalised path, so skipping validation in model_copy is safe
    return settings.model_copy(update={OUTPUT_DIR: output_dir})


def filesystem_resource(bucket_url: str, file_glob: str) -> DltResource:
    """Build a dlt filesystem source over an input directory.

    :param bucket_url: directory holding the input files
    :type  bucket_url: str
    :param file_glob: glob pattern selecting the files to read
    :type  file_glob: str
    :return: a dlt filesystem source over the matching files
    :rtype: DltResource
    """
    return filesystem(bucket_url=bucket_url, file_glob=file_glob)


def run_cli(
    settings_cls: type[LoggerSettings],
    pipeline_fn: Callable[[Any], LoadInfo | None],
    settings_kwargs: dict[str, Any] | None = None,
    cli_args: list[str] | None = None,
) -> LoadInfo | None:
    """Run a DLT pipeline from a generic CLI entry point.

    For CtsSettings subclasses, output_dir is resolved against the dlt config before the settings
    are logged and passed to ``pipeline_fn``.

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
        if isinstance(settings, CtsSettings):
            settings = resolve_cts_settings(settings)
    except (SettingsError, ValidationError, ValueError):
        logger.exception("Error initialising config")
        raise
    except Exception:
        logger.exception("Unexpected error setting up config")
        raise
    return pipeline_fn(settings)


def run_pipeline(  # noqa: PLR0913, PLR0917
    settings: CtsSettings,
    resource: DltResource | list[DltResource],
    destination: TDestinationReferenceArg = None,
    destination_kwargs: dict[str, Any] | None = None,
    pipeline_kwargs: dict[str, Any] | None = None,
    pipeline_run_kwargs: dict[str, Any] | None = None,
) -> LoadInfo | None:
    """Execute a dlt pipeline.

    :param settings: resolved pipeline settings (output_dir must be set)
    :type settings: CtsSettings
    :param resource: dlt resource to run
    :type resource: DltResource | list[DltResource]
    :param destination: the dlt destination to use. If given, the caller is responsible for making
        its location consistent with settings.output_dir.
    :type destination: TDestinationReferenceArg
    :param destination_kwargs: keyword arguments for the dlt destination
    :type destination_kwargs: dict[str, Any] | None
    :param pipeline_kwargs: keyword arguments for the dlt pipeline
    :type pipeline_kwargs: dict[str, Any] | None
    :param pipeline_run_kwargs: keyword arguments for the dlt pipeline run
    :type pipeline_run_kwargs: dict[str, Any] | None
    :raises ValueError: if settings.output_dir has not been resolved
    :return: a load info object, if successful; otherwise None
    :rtype: LoadInfo | None
    """
    if settings.output_dir is None:
        raise ValueError(UNRESOLVED_OUTPUT_DIR)

    pipeline_kwargs = dict(pipeline_kwargs or {})

    # set the output directory for all the pipeline gubbins
    if settings.pipeline_dir:
        pipeline_kwargs["pipelines_dir"] = settings.pipeline_dir

    # update dev mode
    if settings.dlt_dev_mode:
        pipeline_kwargs["dev_mode"] = settings.dlt_dev_mode

    if destination is None:
        destination = dlt.destination(settings.use_destination, **(destination_kwargs or {}))

    pipeline = dlt.pipeline(destination=destination, **pipeline_kwargs)

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
