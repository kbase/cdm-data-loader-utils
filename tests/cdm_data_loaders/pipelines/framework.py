"""Unified test execution framework for running CDM data ingestion pipelines in core, cli, or docker mode."""

import gzip
import json
import logging
import os
import shutil
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime
from enum import StrEnum
from importlib import import_module
from pathlib import Path
from types import ModuleType
from typing import Any, Final, Self
from unittest.mock import patch
from uuid import uuid4

import dlt
import pytest
from dlt.common.pipeline import LoadInfo
from dlt.destinations import filesystem

from cdm_data_loaders.core.settings import CtsSettings, LoggerSettings
from cdm_data_loaders.pipelines import core
from cdm_data_loaders.pipelines.core import LOAD_INFO_TABLE_NAME, PipelineCli
from cdm_data_loaders.utils import cdm_logger
from tests.cdm_data_loaders.pipelines.helpers import settings_to_argv
from tests.conftest import BASE_DIR, CONFIG_BUCKET, LOCAL_FS

DEFAULT_DOCKER_IMAGE_TAG: Final[str] = "cdm-data-loaders-test:latest"
CONTAINER_INPUT_DIR: Final[str] = "/input_dir"
CONTAINER_OUTPUT_DIR: Final[str] = "/output_dir"
CONTAINER_AUX_DIR: Final[str] = "/aux_dir"
CONTAINER_PYMODULES_DIR: Final[str] = "/pymodules"

# Flags whose values are local file paths that must be staged into the container in Docker mode
FILE_PATH_FLAGS: Final[frozenset[str]] = frozenset(
    {
        "--xsd-file",
        "--xml2db-config-file",
        "--schema-file",
        "-s",
        "--log-config-file",
    }
)

# Flags whose values are Python module paths that may be dynamically created in sys.modules
MODULE_PATH_FLAGS: Final[frozenset[str]] = frozenset(
    {
        "--entity-models-module",
        "--schema-files-module",
        "-m",
    }
)

INPUT_DIR_FLAGS: Final[frozenset[str]] = frozenset({"--input-dir", "-i"})
OUTPUT_DIR_FLAGS: Final[frozenset[str]] = frozenset({"--output-dir", "-o"})
DESTINATION_FLAGS: Final[frozenset[str]] = frozenset({"--use-destination", "-d"})

# Registry mapping console-script / pipeline names to (cli_command, module_path)
PIPELINE_TARGETS: Final[dict[str, tuple[str, str]]] = {
    "all_the_bacteria": ("all_the_bacteria", "cdm_data_loaders.pipelines.all_the_bacteria"),
    "jsonlines_extract": ("jsonlines_extract", "cdm_data_loaders.pipelines.jsonlines.extract_pipeline"),
    "jsonlines_jsonschema": (
        "jsonlines_jsonschema",
        "cdm_data_loaders.pipelines.jsonlines.extract_jsonschema_validate_pipeline",
    ),
    "jsonlines_pydantic": (
        "jsonlines_pydantic",
        "cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline",
    ),
    "ncbi_ftp_promote": ("ncbi_ftp_promote", "cdm_data_loaders.pipelines.ncbi_ftp_promote"),
    "ncbi_ftp_sync": ("ncbi_ftp_sync", "cdm_data_loaders.pipelines.ncbi_ftp_download"),
    "ncbi_rest_api": ("ncbi_rest_api", "cdm_data_loaders.pipelines.ncbi_rest_api"),
    "pdb_manifest": ("pdb_manifest", "cdm_data_loaders.pipelines.pdb_manifest"),
    "uniprot": ("uniprot", "cdm_data_loaders.pipelines.uniprot_kb"),
    "uniref": ("uniref", "cdm_data_loaders.pipelines.uniref"),
    "xmltodict": ("xmltodict", "cdm_data_loaders.pipelines.xmltodict.pipeline"),
    "xml2db": ("xml2db", "cdm_data_loaders.pipelines.xml2db.pipeline"),
    "xsv": ("xsv", "cdm_data_loaders.pipelines.xsv.pipeline"),
}

MODULE_TO_CLI_COMMAND: Final[dict[str, str]] = {
    module_path: cli_command for cli_command, module_path in PIPELINE_TARGETS.values()
}

_built_docker_images: set[str] = set()


class ExecutionMode(StrEnum):
    """Supported pipeline test execution modes."""

    CORE = "core"
    CLI = "cli"
    DOCKER = "docker"


EXECUTION_MODES: Final[tuple[Any, ...]] = (
    pytest.param(ExecutionMode.CORE, id="core"),
    pytest.param(ExecutionMode.CLI, id="cli"),
    pytest.param(ExecutionMode.DOCKER, id="docker", marks=pytest.mark.docker),
)


@dataclass(frozen=True, slots=True)
class PipelineSpec:
    """Resolved metadata for executing a pipeline across core, cli, and docker modes."""

    cli_command: str
    module_path: str
    settings_cls: type[LoggerSettings] | None = None
    pipeline_fn: Callable[[Any], LoadInfo | None] | None = None
    cli_callable: Callable[[], Any] | None = None


class _PreservedStrDatetime(datetime):
    """Datetime subclass that preserves the exact string representation written to the report file."""

    _raw_str: str

    def __new__(cls, raw_str: str) -> Self:
        """Create a datetime from an ISO-formatted string while preserving its original string spelling."""
        dt = datetime.fromisoformat(raw_str)
        instance = super().__new__(
            cls,
            dt.year,
            dt.month,
            dt.day,
            dt.hour,
            dt.minute,
            dt.second,
            dt.microsecond,
            tzinfo=dt.tzinfo,
        )
        instance._raw_str = raw_str
        return instance

    def __str__(self) -> str:
        """Return the verbatim timestamp string from the pipeline report."""
        return self._raw_str


@dataclass
class ReconstructedLoadInfo:
    """LoadInfo-compatible object reconstructed from pipeline output files on disk."""

    pipeline: dlt.Pipeline
    dataset_name: str
    loads_ids: list[str]
    started_at: datetime
    has_failed_jobs: bool = False
    finished_at: datetime | None = None
    first_run: bool = False
    destination_name: str = LOCAL_FS
    destination_type: str = "filesystem"
    load_packages: list[Any] = field(default_factory=list)
    metrics: dict[str, Any] = field(default_factory=dict)
    _raw_dict: dict[str, Any] = field(default_factory=dict, repr=False)

    def asdict(self) -> dict[str, Any]:
        """Return a dictionary representation matching dlt's LoadInfo.asdict()."""
        if self._raw_dict:
            return dict(self._raw_dict)
        return {
            "dataset_name": self.dataset_name,
            "destination_name": self.destination_name,
            "destination_type": self.destination_type,
            "finished_at": self.finished_at,
            "first_run": self.first_run,
            "has_failed_jobs": self.has_failed_jobs,
            "loads_ids": list(self.loads_ids),
            "pipeline": {"pipeline_name": self.pipeline.pipeline_name},
            "started_at": self.started_at,
        }


type PipelineIdentifier = str | ModuleType | PipelineCli[Any] | PipelineSpec | Callable[..., LoadInfo | None]


def _spec_from_module(cli_command: str, module: ModuleType) -> PipelineSpec:
    """Build a PipelineSpec by inspecting a pipeline module."""
    cli_attr = getattr(module, "cli", None)
    if isinstance(cli_attr, PipelineCli):
        return PipelineSpec(
            cli_command=cli_command,
            module_path=module.__name__,
            settings_cls=cli_attr.settings_cls,
            pipeline_fn=cli_attr.pipeline_fn,
            cli_callable=cli_attr,
        )
    return PipelineSpec(
        cli_command=cli_command,
        module_path=module.__name__,
        settings_cls=getattr(module, "Settings", None),
        pipeline_fn=getattr(module, "run_pipeline", None),
        cli_callable=cli_attr if callable(cli_attr) else None,
    )


def resolve_pipeline(pipeline: PipelineIdentifier) -> PipelineSpec:
    """Resolve a pipeline name, module, PipelineCli, or callable into a PipelineSpec.

    :param pipeline: pipeline identifier (name, module, PipelineCli, PipelineSpec, or run function)
    :type pipeline: PipelineIdentifier
    :raises ValueError: if the pipeline cannot be resolved
    :return: resolved pipeline specification
    :rtype: PipelineSpec
    """
    if isinstance(pipeline, PipelineSpec):
        return pipeline

    if isinstance(pipeline, PipelineCli):
        module_path = pipeline.pipeline_fn.__module__
        cli_command = MODULE_TO_CLI_COMMAND.get(module_path, module_path.rsplit(".", maxsplit=1)[-1])
        return PipelineSpec(
            cli_command=cli_command,
            module_path=module_path,
            settings_cls=pipeline.settings_cls,
            pipeline_fn=pipeline.pipeline_fn,
            cli_callable=pipeline,
        )

    if isinstance(pipeline, ModuleType):
        module_path = pipeline.__name__
        short_key = module_path.removeprefix("cdm_data_loaders.pipelines.")
        cli_command = MODULE_TO_CLI_COMMAND.get(
            module_path,
            PIPELINE_TARGETS.get(short_key, (short_key.rsplit(".", maxsplit=1)[-1], module_path))[0],
        )
        return _spec_from_module(cli_command, pipeline)

    if callable(pipeline):
        module_path = getattr(pipeline, "__module__", "")
        if module_path in MODULE_TO_CLI_COMMAND:
            module = import_module(module_path)
            return _spec_from_module(MODULE_TO_CLI_COMMAND[module_path], module)
        err_msg = f"Cannot resolve pipeline from callable {pipeline!r} (module {module_path!r})"
        raise ValueError(err_msg)

    if isinstance(pipeline, str):
        normalized = pipeline.strip().removeprefix("cdm_data_loaders.pipelines.")
        if normalized in PIPELINE_TARGETS:
            cli_command, module_path = PIPELINE_TARGETS[normalized]
            module = import_module(module_path)
            return _spec_from_module(cli_command, module)
        full_module_path = f"cdm_data_loaders.pipelines.{normalized}"
        if full_module_path in MODULE_TO_CLI_COMMAND:
            cli_command = MODULE_TO_CLI_COMMAND[full_module_path]
            module = import_module(full_module_path)
            return _spec_from_module(cli_command, module)

    err_msg = f"Unknown pipeline {pipeline!r}. Known pipelines: {', '.join(sorted(PIPELINE_TARGETS))}"
    raise ValueError(err_msg)


def _read_saved_load_info_row(dataset_dir: Path) -> dict[str, Any]:
    """Read the serialized load_info row from ``<dataset_dir>/_load_info`` if present."""
    load_info_dir = dataset_dir / LOAD_INFO_TABLE_NAME
    if not load_info_dir.is_dir():
        return {}
    for file_path in sorted(load_info_dir.glob("*.jsonl*")):
        if not file_path.is_file():
            continue
        with file_path.open("rb") as fh:
            is_gzipped = fh.read(2) == b"\x1f\x8b"
        opener = gzip.open if is_gzipped else open
        with opener(file_path, "rt", encoding="utf-8") as fh:
            for line in fh:
                stripped = line.strip()
                if stripped:
                    return json.loads(stripped)
    return {}


def reconstruct_load_info(output_dir: Path, bucket_dir: Path | None = None) -> LoadInfo | None:
    """Reconstruct a LoadInfo-compatible object and attached dlt.Pipeline from disk artifacts.

    Reads the ``report_<pipeline_e2e_id>.txt`` file written by ``save_pipeline_metadata`` in
    ``output_dir`` and attaches a ``dlt.Pipeline`` pointing at the filesystem destination
    containing ``<dataset_name>``.

    :param output_dir: directory containing ``report_<pipeline_e2e_id>.txt``
    :type output_dir: Path
    :param bucket_dir: optional separate filesystem destination bucket directory
    :type bucket_dir: Path | None
    :return: reconstructed LoadInfo, or None if the report indicates an error or no report exists
    :rtype: LoadInfo | None
    """
    output_dir = Path(output_dir)
    report_files = sorted(output_dir.glob("report_*.txt"), key=lambda path: path.stat().st_mtime)
    if not report_files:
        return None

    report_file = report_files[-1]
    lines = report_file.read_text(encoding="utf-8").splitlines()
    if len(lines) < 3:  # noqa: PLR2004
        return None

    # Line 0: "Pipeline <pipeline_name> dataset <dataset_name>"
    header_prefix = "Pipeline "
    dataset_sep = " dataset "
    if not lines[0].startswith(header_prefix) or dataset_sep not in lines[0]:
        return None
    pipeline_name, dataset_name = lines[0].removeprefix(header_prefix).rsplit(dataset_sep, maxsplit=1)
    pipeline_name = pipeline_name.strip()
    dataset_name = dataset_name.strip()

    # Line 1: "Started: <started_at>"
    started_str = lines[1].removeprefix("Started: ").strip()
    started_at = _PreservedStrDatetime(started_str)

    # Line 2: "Status: completed" | "Status: completed with errors" | "Status: error"
    status_str = lines[2].removeprefix("Status: ").strip()
    if status_str == "error":
        return None
    has_failed_jobs = status_str != "completed"

    # Line 3+: "Completed loads: <id1>, <id2>, ..."
    loads_ids: list[str] = []
    for line in lines[3:]:
        if line.startswith("Completed loads:"):
            loads_raw = line.removeprefix("Completed loads:").strip()
            loads_ids = [item.strip() for item in loads_raw.split(",") if item.strip()]
            break

    # Determine where the dataset directory actually lives on disk
    if (output_dir / dataset_name).is_dir():
        actual_bucket_dir = output_dir
    elif bucket_dir is not None and (Path(bucket_dir) / dataset_name).is_dir():
        actual_bucket_dir = Path(bucket_dir)
    else:
        actual_bucket_dir = Path(bucket_dir) if bucket_dir is not None else output_dir

    saved_row = _read_saved_load_info_row(actual_bucket_dir / dataset_name)
    finished_at_raw = saved_row.get("finished_at")
    finished_at = _PreservedStrDatetime(str(finished_at_raw)) if finished_at_raw else None

    # Use a fresh working directory outside output_dir so dlt syncs the schema from the destination
    # without colliding with stale schemas in ~/.dlt/pipelines or adding files to output_dir.
    pipeline_e2e_id = report_file.stem.removeprefix("report_")
    reconstructed_pipelines_dir = output_dir.parent / f".reconstructed_dlt_{pipeline_e2e_id}"
    pipeline = dlt.pipeline(
        pipeline_name=pipeline_name,
        destination=filesystem(bucket_url=str(actual_bucket_dir)),
        dataset_name=dataset_name,
        pipelines_dir=str(reconstructed_pipelines_dir),
    )

    return ReconstructedLoadInfo(  # pyright: ignore[reportReturnType]
        pipeline=pipeline,
        dataset_name=dataset_name,
        loads_ids=loads_ids,
        started_at=started_at,
        has_failed_jobs=has_failed_jobs,
        finished_at=finished_at,
        first_run=bool(saved_row.get("first_run", False)),
        destination_name=str(saved_row.get("destination_name", LOCAL_FS)),
        destination_type=str(saved_row.get("destination_type", "filesystem")),
        _raw_dict=saved_row,
    )


def ensure_docker_image(image_tag: str | None = None, *, force_rebuild: bool = False) -> str:
    """Build the Docker image from the repository Dockerfile if not already built in this session.

    :param image_tag: optional image tag; defaults to ``cdm-data-loaders-test:latest``
    :type image_tag: str | None
    :param force_rebuild: force running ``docker build`` even if already built in this session
    :type force_rebuild: bool
    :raises RuntimeError: if ``docker build`` fails
    :return: Docker image tag
    :rtype: str
    """
    tag = image_tag or os.environ.get("CDM_TEST_DOCKER_IMAGE") or DEFAULT_DOCKER_IMAGE_TAG
    if not force_rebuild and tag in _built_docker_images:
        return tag

    cmd = ["docker", "build", "-t", tag, str(BASE_DIR.resolve())]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)  # noqa: S603
    if result.returncode != 0:
        err_msg = (
            f"Docker image build failed (exit {result.returncode}):\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
        raise RuntimeError(err_msg)

    _built_docker_images.add(tag)
    return tag


def _extract_flag_value(cli_args: list[str], flags: frozenset[str]) -> str | None:
    """Extract the value immediately following any of ``flags`` in ``cli_args``."""
    for idx, arg in enumerate(cli_args):
        if arg in flags and idx + 1 < len(cli_args):
            return cli_args[idx + 1]
        for flag in flags:
            if arg.startswith(f"{flag}="):
                return arg.split("=", maxsplit=1)[1]
    return None


def _normalize_cli_args(cli_args: list[str]) -> list[str]:
    """Strip 0-byte ``.conf`` placeholder log files created by test fixtures that bypass ``init_logger``."""
    normalized: list[str] = []
    idx = 0
    while idx < len(cli_args):
        arg = cli_args[idx]
        if arg == "--log-config-file" and idx + 1 < len(cli_args):
            log_path = Path(cli_args[idx + 1])
            if log_path.suffix.lower() == ".conf" and log_path.exists() and log_path.stat().st_size == 0:
                idx += 2
                continue
        normalized.append(arg)
        idx += 1
    return normalized


def _materialize_dynamic_module(module_name: str, pymodules_dir: Path) -> bool:
    """If ``module_name`` is an in-memory module in ``sys.modules``, write a ``.py`` file for Docker."""
    module = sys.modules.get(module_name)
    if module is None or getattr(module, "__file__", None) is not None:
        return False

    pymodules_dir.mkdir(parents=True, exist_ok=True)
    target_file = pymodules_dir / f"{module_name}.py"

    if hasattr(module, "ENTITY_SCHEMAS"):
        schemas = module.ENTITY_SCHEMAS
        target_file.write_text(f"ENTITY_SCHEMAS = {schemas!r}\n", encoding="utf-8")
        return True

    if hasattr(module, "ENTITY_MODELS"):
        models: dict[str, Any] = module.ENTITY_MODELS
        imports: list[str] = []
        entries: list[str] = []
        for idx, (table_name, model_cls) in enumerate(models.items()):
            alias = f"_Model_{idx}"
            imports.append(f"from {model_cls.__module__} import {model_cls.__name__} as {alias}")
            entries.append(f"    {table_name!r}: {alias},")
        code = "\n".join([*imports, "", "ENTITY_MODELS = {", *entries, "}", ""])
        target_file.write_text(code, encoding="utf-8")
        return True

    return False


def _write_docker_dlt_config(config_path: Path, use_destination: str) -> None:
    """Write ``config.toml`` for the container with ``/input_dir`` and ``/output_dir`` buckets."""
    config_path.parent.mkdir(parents=True, exist_ok=True)
    destinations = {"cts", LOCAL_FS, use_destination}
    blocks: list[str] = []
    for dest in sorted(destinations):
        blocks.append(f'[destination.{dest}]\nbucket_url = "{CONTAINER_OUTPUT_DIR}"\ndestination_type = "filesystem"\n')
        blocks.append(f'[sources.{dest}]\nbucket_url = "{CONTAINER_INPUT_DIR}"\n')
    config_path.write_text("\n".join(blocks), encoding="utf-8")


class PipelineTestRunner:
    """Execute any CDM ingestion pipeline in ``core``, ``cli``, or ``docker`` mode."""

    def __init__(self, mode: ExecutionMode | str = ExecutionMode.CORE, docker_image: str | None = None) -> None:
        """Initialize the runner with the target execution mode and optional Docker image tag."""
        self.mode = ExecutionMode(mode)
        self.docker_image = docker_image

    def __call__(
        self,
        pipeline: PipelineIdentifier,
        settings: CtsSettings | LoggerSettings | None = None,
        *,
        cli_args: list[str] | None = None,
        **overrides: Any,  # noqa: ANN401
    ) -> LoadInfo | None:
        """Alias for :meth:`run` so the runner can be passed as a ``run_fn`` callback."""
        return self.run(pipeline, settings, cli_args=cli_args, **overrides)

    def run(
        self,
        pipeline: PipelineIdentifier,
        settings: CtsSettings | LoggerSettings | None = None,
        *,
        cli_args: list[str] | None = None,
        **overrides: Any,  # noqa: ANN401
    ) -> LoadInfo | None:
        """Execute ``pipeline`` in the configured mode and return its ``LoadInfo``.

        :param pipeline: pipeline name, module, PipelineCli, PipelineSpec, or pipeline run function
        :type pipeline: PipelineIdentifier
        :param settings: optional pre-built settings instance
        :type settings: CtsSettings | LoggerSettings | None
        :param cli_args: optional explicit CLI argument list (excluding the program name)
        :type cli_args: list[str] | None
        :param overrides: optional field overrides applied if building or updating settings
        :type overrides: Any
        :raises ValueError: if neither settings nor cli_args are provided and settings cannot be built
        :return: pipeline LoadInfo (or ReconstructedLoadInfo in cli/docker modes)
        :rtype: LoadInfo | None
        """
        spec = resolve_pipeline(pipeline)
        effective_settings = settings
        if effective_settings is not None and overrides:
            effective_settings = effective_settings.model_copy(update=overrides)
        elif effective_settings is None and cli_args is None:
            if spec.settings_cls is None:
                err_msg = f"Pipeline {spec.cli_command!r} has no settings_cls; supply cli_args explicitly."
                raise ValueError(err_msg)
            effective_settings = spec.settings_cls(**overrides)

        raw_args = list(cli_args) if cli_args is not None else settings_to_argv(effective_settings)
        # Strip leading program name if the caller passed a full sys.argv list
        if raw_args and not raw_args[0].startswith("-"):
            raw_args = raw_args[1:]
        args = _normalize_cli_args(raw_args)

        output_dir_str = _extract_flag_value(args, OUTPUT_DIR_FLAGS) or (
            str(getattr(effective_settings, "output_dir", "")) if effective_settings else None
        )
        if not output_dir_str:
            err_msg = "Could not determine output_dir from settings or cli_args"
            raise ValueError(err_msg)
        output_dir = Path(output_dir_str)
        output_dir.mkdir(parents=True, exist_ok=True)

        use_destination = (
            _extract_flag_value(args, DESTINATION_FLAGS)
            or (str(getattr(effective_settings, "use_destination", LOCAL_FS)) if effective_settings else None)
            or LOCAL_FS
        )

        if self.mode == ExecutionMode.CORE:
            return self._run_core(spec, args, output_dir, use_destination)
        if self.mode == ExecutionMode.CLI:
            return self._run_cli(spec, args, output_dir, use_destination)
        return self._run_docker(spec, args, output_dir, use_destination)

    def _resolve_bucket_dir(self, output_dir: Path, use_destination: str) -> Path:
        """Determine the active host bucket directory for ``use_destination``."""
        configured = dlt.config.get(f"destination.{use_destination}.bucket_url")
        default_buckets = {CONFIG_BUCKET.get(LOCAL_FS), "/configured/out/", None}
        if configured and configured not in default_buckets and Path(str(configured)).exists():
            return Path(str(configured))
        return output_dir

    def _run_core(
        self,
        spec: PipelineSpec,
        args: list[str],
        output_dir: Path,
        use_destination: str,
    ) -> LoadInfo | None:
        """Execute the pipeline in-process via ``cdm_data_loaders.pipelines.core.run_cli``."""
        if spec.settings_cls is None or spec.pipeline_fn is None:
            err_msg = f"Pipeline {spec.cli_command!r} does not expose settings_cls and pipeline_fn for core.run_cli"
            raise ValueError(err_msg)

        bucket_dir = self._resolve_bucket_dir(output_dir, use_destination)
        root = logging.getLogger()
        handlers, level = root.handlers[:], root.level
        applied_key = cdm_logger.applied_key
        try:
            with dlt.config.values(
                {
                    f"destination.{use_destination}.destination_type": "filesystem",
                    f"destination.{use_destination}.bucket_url": str(bucket_dir),
                }
            ):
                return core.run_cli(spec.settings_cls, spec.pipeline_fn, cli_args=args)
        finally:
            for handler in root.handlers[:]:
                root.removeHandler(handler)
            for handler in handlers:
                root.addHandler(handler)
            root.setLevel(level)
            cdm_logger.applied_key = applied_key

    def _run_cli(
        self,
        spec: PipelineSpec,
        args: list[str],
        output_dir: Path,
        use_destination: str,
    ) -> LoadInfo | None:
        """Execute the pipeline via ``cdm_data_loaders.pipelines.<pipeline>.cli()`` and reconstruct LoadInfo."""
        if spec.cli_callable is None:
            err_msg = f"Pipeline {spec.cli_command!r} does not expose a cli() entry point"
            raise ValueError(err_msg)

        bucket_dir = self._resolve_bucket_dir(output_dir, use_destination)
        root = logging.getLogger()
        handlers, level = root.handlers[:], root.level
        applied_key = cdm_logger.applied_key
        try:
            with (
                patch.object(sys, "argv", [spec.cli_command, *args]),
                dlt.config.values(
                    {
                        f"destination.{use_destination}.destination_type": "filesystem",
                        f"destination.{use_destination}.bucket_url": str(bucket_dir),
                    }
                ),
            ):
                spec.cli_callable()
        finally:
            for handler in root.handlers[:]:
                root.removeHandler(handler)
            for handler in handlers:
                root.addHandler(handler)
            root.setLevel(level)
            cdm_logger.applied_key = applied_key

        return reconstruct_load_info(output_dir, bucket_dir=bucket_dir)

    def _run_docker(
        self,
        spec: PipelineSpec,
        args: list[str],
        output_dir: Path,
        use_destination: str,
    ) -> LoadInfo | None:
        """Execute the pipeline inside the Docker container and reconstruct LoadInfo from mounted output."""
        image_tag = ensure_docker_image(self.docker_image)
        input_dir_str = _extract_flag_value(args, INPUT_DIR_FLAGS)
        if not input_dir_str:
            err_msg = "Docker execution mode requires --input-dir / -i in settings or cli_args"
            raise ValueError(err_msg)
        host_input_dir = Path(input_dir_str).resolve()
        host_output_dir = output_dir.resolve()

        stage_dir = host_output_dir.parent / f".docker_stage_{uuid4().hex[:8]}"
        staged_input_dir = stage_dir / "input_dir"
        aux_dir = stage_dir / "aux"
        pymodules_dir = stage_dir / "pymodules"
        staged_input_dir.mkdir(parents=True, exist_ok=True)

        if host_input_dir.is_dir():
            shutil.copytree(host_input_dir, staged_input_dir, dirs_exist_ok=True)

        _write_docker_dlt_config(staged_input_dir / "dlt_conf" / "config.toml", use_destination)

        container_args = self._translate_args_for_docker(args, aux_dir, pymodules_dir)

        docker_cmd: list[str] = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{staged_input_dir}:{CONTAINER_INPUT_DIR}",
            "-v",
            f"{host_output_dir}:{CONTAINER_OUTPUT_DIR}",
            "-v",
            f"{(BASE_DIR / 'tests').resolve()}:/app/tests:ro",
        ]
        if aux_dir.is_dir():
            docker_cmd.extend(["-v", f"{aux_dir.resolve()}:{CONTAINER_AUX_DIR}:ro"])
        if pymodules_dir.is_dir():
            docker_cmd.extend(
                [
                    "-v",
                    f"{pymodules_dir.resolve()}:{CONTAINER_PYMODULES_DIR}:ro",
                    "-e",
                    f"PYTHONPATH={CONTAINER_PYMODULES_DIR}:/app/src:/app",
                ]
            )
        else:
            docker_cmd.extend(["-e", "PYTHONPATH=/app/src:/app"])

        docker_cmd.extend([image_tag, "pipeline", spec.cli_command, *container_args])

        result = subprocess.run(docker_cmd, capture_output=True, text=True, check=False)  # noqa: S603
        if result.returncode != 0:
            stderr_and_stdout = f"{result.stdout}\n{result.stderr}".strip()
            if "Unknown table name" in stderr_and_stdout or "ValueError" in stderr_and_stdout:
                raise ValueError(stderr_and_stdout)
            err_msg = (
                f"Docker pipeline '{spec.cli_command}' failed with exit code {result.returncode}:\n"
                f"Command: {' '.join(docker_cmd)}\n"
                f"STDOUT:\n{result.stdout}\n"
                f"STDERR:\n{result.stderr}"
            )
            raise RuntimeError(err_msg)

        return reconstruct_load_info(host_output_dir, bucket_dir=host_output_dir)

    def _translate_args_for_docker(
        self,
        args: list[str],
        aux_dir: Path,
        pymodules_dir: Path,
    ) -> list[str]:
        """Rewrite host file paths and materialize dynamic modules for container execution."""
        translated: list[str] = []
        idx = 0
        aux_counter = 0
        while idx < len(args):
            token = args[idx]
            if token in INPUT_DIR_FLAGS and idx + 1 < len(args):
                translated.extend([token, CONTAINER_INPUT_DIR])
                idx += 2
                continue
            if token in OUTPUT_DIR_FLAGS and idx + 1 < len(args):
                translated.extend([token, CONTAINER_OUTPUT_DIR])
                idx += 2
                continue
            if token in FILE_PATH_FLAGS and idx + 1 < len(args):
                host_file = Path(args[idx + 1]).resolve()
                aux_dir.mkdir(parents=True, exist_ok=True)
                staged_name = f"{aux_counter}_{host_file.name}"
                aux_counter += 1
                shutil.copy2(host_file, aux_dir / staged_name)
                translated.extend([token, f"{CONTAINER_AUX_DIR}/{staged_name}"])
                idx += 2
                continue
            if token in MODULE_PATH_FLAGS and idx + 1 < len(args):
                module_name = args[idx + 1]
                _materialize_dynamic_module(module_name, pymodules_dir)
                translated.extend([token, module_name])
                idx += 2
                continue
            translated.append(token)
            idx += 1
        return translated
