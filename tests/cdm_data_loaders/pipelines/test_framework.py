"""Unit and integration tests for the unified pipeline test execution framework."""

import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest
from pydantic import BaseModel

from cdm_data_loaders.pipelines import uniprot_kb_cdm
from cdm_data_loaders.pipelines.jsonlines.extract_pipeline import (
    cli as jsonl_extract_cli,
)
from cdm_data_loaders.pipelines.jsonlines.extract_pipeline import (
    run_jsonlines_pipeline,
)
from cdm_data_loaders.pipelines.jsonlines.settings import JsonlExtractSettings, JsonlJsonschemaIngestSettings
from cdm_data_loaders.pipelines.xmltodict.settings import XmlToDictSettings
from tests.cdm_data_loaders.pipelines import framework
from tests.cdm_data_loaders.pipelines.conftest import DEFAULT_DLT_TABLES
from tests.cdm_data_loaders.pipelines.framework import (
    CONTAINER_AUX_DIR,
    CONTAINER_INPUT_DIR,
    CONTAINER_OUTPUT_DIR,
    CONTAINER_PYMODULES_DIR,
    PIPELINE_TARGETS,
    ExecutionMode,
    PipelineSpec,
    PipelineTestRunner,
    ensure_docker_image,
    reconstruct_load_info,
    resolve_pipeline,
)
from tests.cdm_data_loaders.pipelines.helpers import (
    assert_pipeline_report_written,
    assert_run_ok,
    read_data_rows,
    settings_to_argv,
)
from tests.conftest import LOCAL_FS


class _SampleEntityModel(BaseModel):
    """Sample Pydantic model for dynamic module materialization tests."""

    record_id: str


@pytest.mark.parametrize("target_key", sorted(PIPELINE_TARGETS))
def test_resolve_pipeline_pass_all_registered_targets(target_key: str) -> None:
    """Every registered pipeline resolves to a valid PipelineSpec with a synchronized cli_prog_name."""
    spec = resolve_pipeline(target_key)
    expected_cli, expected_module = PIPELINE_TARGETS[target_key]
    assert spec.cli_command == expected_cli
    assert spec.module_path == expected_module
    assert callable(spec.cli_callable)
    assert spec.settings_cls is not None
    assert spec.settings_cls.model_config.get("cli_prog_name") == expected_cli


def test_resolve_pipeline_pass_from_module_cli_spec_and_callable() -> None:
    """resolve_pipeline accepts ModuleType, PipelineCli, PipelineSpec, and pipeline run functions."""
    from_module = resolve_pipeline(uniprot_kb_cdm)
    assert from_module.cli_command == "uniprot_cdm"
    assert from_module.settings_cls is uniprot_kb_cdm.UniProtSettings
    assert from_module.pipeline_fn is uniprot_kb_cdm.run_uniprot_pipeline

    from_cli = resolve_pipeline(jsonl_extract_cli)
    assert from_cli.cli_command == "jsonlines_extract"
    assert from_cli.settings_cls is JsonlExtractSettings
    assert from_cli.pipeline_fn is run_jsonlines_pipeline

    from_callable = resolve_pipeline(run_jsonlines_pipeline)
    assert from_callable.cli_command == "jsonlines_extract"
    assert from_callable.pipeline_fn is run_jsonlines_pipeline

    from_full_module_str = resolve_pipeline("cdm_data_loaders.pipelines.uniprot_kb_cdm")
    assert from_full_module_str.cli_command == "uniprot_cdm"

    from_spec = resolve_pipeline(from_cli)
    assert from_spec is from_cli


def test_resolve_pipeline_fail_unknown_identifier() -> None:
    """resolve_pipeline raises ValueError for unknown pipeline names or unmapped callables."""
    with pytest.raises(ValueError, match="Unknown pipeline 'nonexistent_pipeline'"):
        resolve_pipeline("nonexistent_pipeline")

    with pytest.raises(ValueError, match="Cannot resolve pipeline from callable"):
        resolve_pipeline(lambda _settings: None)


def test_settings_to_argv_pass_excludes_computed_fields_and_expands_lists(tmp_path: Path) -> None:
    """settings_to_argv omits computed fields and None values, formats bools, and expands list fields."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    settings = JsonlJsonschemaIngestSettings(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        use_destination=LOCAL_FS,
        use_output_dir_for_pipeline_metadata=True,
        dlt_dev_mode=False,
        dataset_name="custom_ds",
        loader_file_format="jsonl",
        schema_files_module="some.module",
        table_names=["alpha", "beta"],
    )
    argv = settings_to_argv(settings)
    assert "--pipeline-dir" not in argv
    assert "--log-config-file" not in argv
    assert "--use-output-dir-for-pipeline-metadata" in argv
    idx_meta = argv.index("--use-output-dir-for-pipeline-metadata")
    assert argv[idx_meta + 1] == "true"
    idx_tables = argv.index("--table-names")
    assert argv[idx_tables + 1 : idx_tables + 3] == ["alpha", "beta"]


def test_reconstruct_load_info_pass_and_error_cases(tmp_path: Path) -> None:
    """reconstruct_load_info returns None on missing or error reports and reconstructs completed reports."""
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    assert reconstruct_load_info(empty_dir) is None

    error_dir = tmp_path / "error_out"
    error_dir.mkdir()
    (error_dir / "report_123_456.txt").write_text(
        "Pipeline: test_pipe\n"
        "Dataset: test_ds\n"
        "Started: 2026-01-01 00:00:00+00:00\n"
        "Finished: 2026-01-01 00:00:05+00:00\n"
        "Status: error\n"
        "Pipeline failed\n",
        encoding="utf-8",
    )
    assert reconstruct_load_info(error_dir) is None

    completed_dir = tmp_path / "completed_out"
    completed_dir.mkdir()
    started_str = "2026-01-01 12:34:56.789000+00:00"
    finished_str = "2026-01-01 12:35:00.000000+00:00"
    (completed_dir / "report_111_222.txt").write_text(
        f"Pipeline: test_pipe\n"
        f"Dataset: test_ds\n"
        f"Started: {started_str}\n"
        f"Finished: {finished_str}\n"
        f"Status: completed with errors\n"
        f"Completed loads: 111, 222\n"
        f"Pipeline load info saved\n",
        encoding="utf-8",
    )
    reconstructed = reconstruct_load_info(completed_dir)
    assert reconstructed is not None
    assert reconstructed.pipeline.pipeline_name == "test_pipe"
    assert reconstructed.dataset_name == "test_ds"
    assert str(reconstructed.started_at) == started_str
    assert str(reconstructed.finished_at) == finished_str
    assert reconstructed.has_failed_jobs is True
    assert reconstructed.loads_ids == ["111", "222"]
    as_dict = reconstructed.asdict()
    assert as_dict["dataset_name"] == "test_ds"
    assert as_dict["loads_ids"] == ["111", "222"]


def test_ensure_docker_image_pass_caches_and_fail_raises(monkeypatch: pytest.MonkeyPatch) -> None:
    """ensure_docker_image builds once per tag and raises RuntimeError when docker build fails."""
    monkeypatch.setattr(framework, "_built_docker_images", set())
    calls: list[list[str]] = []

    def fake_run(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(cmd)
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="built", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    tag = ensure_docker_image("custom-image:test")
    assert tag == "custom-image:test"
    assert len(calls) == 1
    assert calls[0][:4] == ["docker", "build", "-t", "custom-image:test"]

    # Second call uses cache without invoking subprocess.run
    assert ensure_docker_image("custom-image:test") == "custom-image:test"
    assert len(calls) == 1

    # Failing build raises RuntimeError
    def fake_fail(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(cmd, returncode=1, stdout="", stderr="build boom")

    monkeypatch.setattr(subprocess, "run", fake_fail)
    with pytest.raises(RuntimeError, match="Docker image build failed"):
        ensure_docker_image("failing-image:test")


@pytest.mark.parametrize("mode", [ExecutionMode.CORE, ExecutionMode.CLI], ids=str)
def test_pipeline_test_runner_pass_core_and_cli_modes_end_to_end(
    tmp_path: Path,
    mode: ExecutionMode,
) -> None:
    """PipelineTestRunner executes a pipeline in core and cli modes and hooks into all result assertions."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    expected_records = [{"id": "rec_1", "value": 10}, {"id": "rec_2", "value": 20}]
    (input_dir / "data.jsonl").write_text(
        "\n".join(f'{{"id": "{r["id"]}", "value": {r["value"]}}}' for r in expected_records) + "\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    settings = JsonlExtractSettings(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        use_destination=LOCAL_FS,
        use_output_dir_for_pipeline_metadata=True,
        dlt_dev_mode=False,
        dataset_name="framework_test_ds",
        table_name="items",
        loader_file_format="jsonl",
    )

    runner = PipelineTestRunner(mode=mode)
    load_info = assert_run_ok(
        runner(
            "jsonlines_extract",
            settings=settings,
        )
    )

    assert_pipeline_report_written(output_dir, load_info, metadata_moved=True)
    assert set(load_info.pipeline.dataset().tables) >= set(DEFAULT_DLT_TABLES)
    assert read_data_rows(output_dir / settings.dataset_name) == {"items": expected_records}


def test_pipeline_test_runner_pass_docker_mode_command_and_staging(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Docker mode stages input, dlt config, auxiliary files, dynamic modules, and reconstructs LoadInfo."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    (input_dir / "sample.xml").write_text("<root><entry><id>1</id></entry></root>", encoding="utf-8")
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    xsd_file = tmp_path / "schema.xsd"
    xsd_file.write_text("<xs:schema xmlns:xs='http://www.w3.org/2001/XMLSchema'/>", encoding="utf-8")

    # Create dynamic in-memory modules for both ENTITY_SCHEMAS and ENTITY_MODELS
    dyn_schema_mod_name = "_test_dynamic_schema_mod"
    dyn_schema_mod = ModuleType(dyn_schema_mod_name)
    dyn_schema_mod.ENTITY_SCHEMAS = {"widget": {"type": "object"}}  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, dyn_schema_mod_name, dyn_schema_mod)

    dyn_model_mod_name = "_test_dynamic_model_mod"
    dyn_model_mod = ModuleType(dyn_model_mod_name)
    dyn_model_mod.ENTITY_MODELS = {"widget": _SampleEntityModel}  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, dyn_model_mod_name, dyn_model_mod)

    monkeypatch.setattr(framework, "_built_docker_images", {"test-image:latest"})
    captured_cmds: list[list[str]] = []

    def fake_docker_run(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        captured_cmds.append(cmd)
        # Simulate the container writing a completed report to the mounted output_dir
        (output_dir / "report_999_1000.txt").write_text(
            "Pipeline: xmltodict\n"
            "Dataset: xml_docker_ds\n"
            "Started: 2026-02-01 00:00:00+00:00\n"
            "Finished: 2026-02-01 00:00:05+00:00\n"
            "Status: completed\n"
            "Completed loads: 999\n"
            "Pipeline load info saved\n",
            encoding="utf-8",
        )
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_docker_run)

    settings = XmlToDictSettings(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        use_destination=LOCAL_FS,
        use_output_dir_for_pipeline_metadata=True,
        dlt_dev_mode=False,
        dataset_name="xml_docker_ds",
        table_name="entries",
        xml_tag="entry",
        loader_file_format="parquet",
    )

    runner = PipelineTestRunner(mode=ExecutionMode.DOCKER, docker_image="test-image:latest")
    cli_args = [
        *settings_to_argv(settings),
        "--xsd-file",
        str(xsd_file),
        "--schema-files-module",
        dyn_schema_mod_name,
        "--entity-models-module",
        dyn_model_mod_name,
    ]
    load_info = runner.run("xmltodict", cli_args=cli_args)

    assert load_info is not None
    assert load_info.pipeline.pipeline_name == "xmltodict"
    assert load_info.dataset_name == "xml_docker_ds"
    assert len(captured_cmds) == 1

    docker_cmd = captured_cmds[0]
    assert docker_cmd[:3] == ["docker", "run", "--rm"]
    idx_image = docker_cmd.index("test-image:latest")
    assert docker_cmd[idx_image : idx_image + 3] == ["test-image:latest", "pipeline", "xmltodict"]

    container_args = docker_cmd[idx_image + 3 :]
    idx_in = container_args.index("--input-dir")
    assert container_args[idx_in + 1] == CONTAINER_INPUT_DIR
    idx_out = container_args.index("--output-dir")
    assert container_args[idx_out + 1] == CONTAINER_OUTPUT_DIR
    idx_xsd = container_args.index("--xsd-file")
    assert container_args[idx_xsd + 1].startswith(f"{CONTAINER_AUX_DIR}/")
    assert f"PYTHONPATH={CONTAINER_PYMODULES_DIR}:/app/src:/app" in docker_cmd


def test_pipeline_test_runner_fail_docker_mode_nonzero_exit(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Docker mode raises ValueError on validation errors and RuntimeError on general container failures."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    monkeypatch.setattr(framework, "_built_docker_images", {"test-image:latest"})

    def fake_fail_runtime(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(cmd, returncode=2, stdout="", stderr="container crashed")

    monkeypatch.setattr(subprocess, "run", fake_fail_runtime)
    runner = PipelineTestRunner(mode=ExecutionMode.DOCKER, docker_image="test-image:latest")
    with pytest.raises(RuntimeError, match="failed with exit code 2"):
        runner.run(
            "jsonlines_extract",
            cli_args=["--input-dir", str(input_dir), "--output-dir", str(output_dir)],
        )

    def fake_fail_value_error(cmd: list[str], **_kwargs: object) -> subprocess.CompletedProcess[str]:
        return subprocess.CompletedProcess(cmd, returncode=1, stdout="", stderr="ValueError: Unknown table name: foo")

    monkeypatch.setattr(subprocess, "run", fake_fail_value_error)
    with pytest.raises(ValueError, match="Unknown table name"):
        runner.run(
            "jsonlines_extract",
            cli_args=["--input-dir", str(input_dir), "--output-dir", str(output_dir)],
        )


def test_pipeline_test_runner_fail_missing_required_arguments(tmp_path: Path) -> None:
    """PipelineTestRunner raises ValueError when output_dir or input_dir (in docker mode) is missing."""
    runner = PipelineTestRunner(mode=ExecutionMode.CORE)
    with pytest.raises(ValueError, match="Could not determine output_dir"):
        runner.run("jsonlines_extract", cli_args=["--input-dir", str(tmp_path)])

    spec_no_settings = PipelineSpec(cli_command="custom", module_path="custom.mod")
    with pytest.raises(ValueError, match="has no settings_cls"):
        runner.run(spec_no_settings)

    with pytest.raises(ValueError, match="does not expose settings_cls and pipeline_fn"):
        runner.run(spec_no_settings, cli_args=["--output-dir", str(tmp_path)])

    cli_runner = PipelineTestRunner(mode=ExecutionMode.CLI)
    with pytest.raises(ValueError, match="does not expose a cli\\(\\) entry point"):
        cli_runner.run(spec_no_settings, cli_args=["--output-dir", str(tmp_path)])

    docker_runner = PipelineTestRunner(mode=ExecutionMode.DOCKER, docker_image="test-image:latest")
    monkeypatch = pytest.MonkeyPatch()
    monkeypatch.setattr(framework, "_built_docker_images", {"test-image:latest"})
    try:
        with pytest.raises(ValueError, match="Docker execution mode requires --input-dir"):
            docker_runner.run("jsonlines_extract", cli_args=["--output-dir", str(tmp_path)])
    finally:
        monkeypatch.undo()


def test_resolve_pipeline_pass_legacy_module_without_pipeline_cli() -> None:
    """resolve_pipeline handles modules whose cli attribute is a plain function rather than PipelineCli."""
    legacy_mod = ModuleType("cdm_data_loaders.pipelines.custom_legacy")
    legacy_mod.Settings = JsonlExtractSettings  # type: ignore[attr-defined]
    legacy_mod.run_pipeline = run_jsonlines_pipeline  # type: ignore[attr-defined]
    legacy_mod.cli = lambda: None  # type: ignore[attr-defined]
    spec = resolve_pipeline(legacy_mod)
    assert spec.cli_command == "custom_legacy"
    assert spec.settings_cls is JsonlExtractSettings
    assert spec.pipeline_fn is run_jsonlines_pipeline
    assert spec.cli_callable is legacy_mod.cli


@pytest.mark.docker
def test_pipeline_test_runner_pass_real_docker_container(tmp_path: Path, docker_image: str) -> None:
    """Execute a real Docker container run via PipelineTestRunner when --run-docker is enabled."""
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    expected_records = [{"id": "docker_1", "value": 42}]
    (input_dir / "data.jsonl").write_text('{"id": "docker_1", "value": 42}\n', encoding="utf-8")
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    settings = JsonlExtractSettings(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        use_destination=LOCAL_FS,
        use_output_dir_for_pipeline_metadata=True,
        dlt_dev_mode=False,
        dataset_name="docker_framework_ds",
        table_name="items",
        loader_file_format="jsonl",
    )

    runner = PipelineTestRunner(mode=ExecutionMode.DOCKER, docker_image=docker_image)
    load_info = assert_run_ok(runner.run("jsonlines_extract", settings=settings))

    assert_pipeline_report_written(output_dir, load_info, metadata_moved=True)
    assert set(load_info.pipeline.dataset().tables) >= set(DEFAULT_DLT_TABLES)
    assert read_data_rows(output_dir / settings.dataset_name) == {"items": expected_records}
