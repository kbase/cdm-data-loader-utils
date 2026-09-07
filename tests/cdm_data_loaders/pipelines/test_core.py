"""Tests for the shared core DLT pipeline functions."""

import logging
import os
from copy import deepcopy
from pathlib import Path
from typing import Any, Final
from unittest.mock import MagicMock, call, patch

import dlt
import pytest
from pydantic import ValidationError
from pydantic_settings import SettingsError

from cdm_data_loaders.core.fields import (
    DEV_MODE,
    OUTPUT_DIR,
    USE_DESTINATION,
    VALID_DESTINATIONS,
)
from cdm_data_loaders.core.settings import (
    BatchedFileInputSettings,
    CtsSettings,
    InputOutputSettings,
    LoggerSettings,
)
from cdm_data_loaders.pipelines import core
from cdm_data_loaders.pipelines.core import (
    NO_MESSAGE,
    WEBHOOK_NOT_CONFIGURED,
    construct_env_var,
    run_cli,
    run_pipeline,
    send_slack_message_carefully,
    sync_configs,
)
from tests.cdm_data_loaders.core.test_settings import SETTINGS_CLASSES, TEST_CTS_SETTINGS


def make_batched_settings(**kwargs: str | int) -> BatchedFileInputSettings:
    """Generate a validated BatchedFileInputSettings object with a valid dlt config."""
    return BatchedFileInputSettings.model_validate(kwargs)


@pytest.fixture
def empty_dlt_config() -> dict[str, Any]:
    """A completely empty dlt config dict.

    Represents dlt.config before any environment variables or config files have contributed
    values, so that sync_configs tests can confirm behaviour is driven purely by the settings
    object and not by any pre-existing config state.
    """
    return {}


@pytest.fixture
def test_bfi_settings(tmp_path: Path) -> BatchedFileInputSettings:
    """Minimal valid BatchedFileInputSettings (no start_at, no output_dir)."""
    return make_batched_settings(input_dir="/fake/input", output_dir=str(tmp_path))


@pytest.fixture
def test_cts_settings() -> CtsSettings:
    """A fully validated CtsSettings instance using the test dlt config."""
    return CtsSettings.model_validate(TEST_CTS_SETTINGS)


@pytest.fixture(
    params=[
        pytest.param({"input_dir": "/fake/input"}, id="default"),
        pytest.param(
            {
                "input_dir": "/path/to/dir",
                "use_destination": "local_fs",
                "start_at": 15,
                "output_dir": "/some/dir",
            },
            id="alt",
        ),
    ]
)
def config(request: pytest.FixtureRequest) -> BatchedFileInputSettings:
    """Parametrized fixture providing default and non-default settings."""
    return make_batched_settings(**request.param)


@pytest.fixture(
    params=[
        pytest.param(lambda: LoggerSettings.model_validate({}), id="no-relevant-attrs"),
        pytest.param(
            lambda: InputOutputSettings.model_validate({"input_dir": "/fake/input", "output_dir": "/fake/output"}),
            id="output_dir-only-no-destination",
        ),
    ]
)
def settings_missing_sync_attrs(request: pytest.FixtureRequest) -> LoggerSettings | InputOutputSettings:
    """Settings instances missing one or more attribute that sync_configs checks for via hasattr."""
    return request.param()


def assert_pipeline_run_correctly(
    mock_dlt: MagicMock,
    fake_resource: MagicMock,
    destination: str,
    destination_kwargs: dict[str, Any] | None,
    pipeline_kwargs: dict[str, Any] | None,
    pipeline_run_kwargs: dict[str, Any] | None,
) -> None:
    """Shared assertion block for run_pipeline tests."""
    assert mock_dlt.destination.call_args_list == [call(destination, **destination_kwargs or {})]
    assert mock_dlt.pipeline.call_args_list == [
        call(destination=mock_dlt.destination.return_value, **pipeline_kwargs or {})
    ]
    mock_pipeline = mock_dlt.pipeline.return_value
    assert mock_pipeline.run.call_args_list == [call(fake_resource, **pipeline_run_kwargs or {})]


# construct_env_var
@pytest.mark.parametrize(
    ("b_var", "t_var", "char_str", "expected"),
    [
        ("B123", "T456", "C789", "https://hooks.slack.com/services/B123/T456/C789/"),
        ("one", "two", "three", "https://hooks.slack.com/services/one/two/three/"),
    ],
)
def test_construct_env_var_sets_runtime_slack_incoming_hook(
    b_var: str, t_var: str, char_str: str, expected: str
) -> None:
    """Test successful setting of the RUNTIME__SLACK_INCOMING_HOOK env var when all three variables are present."""
    with patch.dict(
        os.environ,
        {
            "VARIABLE_B": b_var,
            "VARIABLE_T": t_var,
            "CHAR_STR": char_str,
        },
        clear=True,
    ):
        # function returns None
        assert construct_env_var() is None
        assert os.environ["RUNTIME__SLACK_INCOMING_HOOK"] == expected


@pytest.mark.parametrize(
    ("b_var", "t_var", "char_str"),
    [
        (None, "T456", "C789"),
        ("B123", None, "C789"),
        ("B123", "T456", None),
        (None, None, None),
        ("", "T456", "C789"),
    ],
)
def test_construct_env_var_does_not_set_when_any_variable_missing(
    b_var: str | None, t_var: str | None, char_str: str | None
) -> None:
    """Test unsuccessful interpolation of env vars if any of the three variables are missing or empty."""
    env = {}
    if b_var is not None:
        env["VARIABLE_B"] = b_var
    if t_var is not None:
        env["VARIABLE_T"] = t_var
    if char_str is not None:
        env["CHAR_STR"] = char_str

    with patch.dict(os.environ, env, clear=True):
        assert construct_env_var() is None
        assert "RUNTIME__SLACK_INCOMING_HOOK" not in os.environ


def test_construct_env_var_overwrites_existing_hook() -> None:
    """An existing RUNTIME__SLACK_INCOMING_HOOK is overwritten when all vars are present."""
    with patch.dict(
        os.environ,
        {
            "VARIABLE_B": "B",
            "VARIABLE_T": "T",
            "CHAR_STR": "C",
            "RUNTIME__SLACK_INCOMING_HOOK": "https://old.hook/",
        },
        clear=True,
    ):
        assert construct_env_var() is None
        assert os.environ["RUNTIME__SLACK_INCOMING_HOOK"] == "https://hooks.slack.com/services/B/T/C/"


# send slack message carefully
SLACK_HOOK: Final[str] = "https://slack.hook/what/ever"
TEST_MESSAGE: Final[str] = "3... 2... 1... testing?"


@pytest.mark.parametrize(
    ("slack_hook", "message", "err_msg"),
    [
        ("", TEST_MESSAGE, WEBHOOK_NOT_CONFIGURED),
        (SLACK_HOOK, "", NO_MESSAGE),
        ("", "", WEBHOOK_NOT_CONFIGURED),
    ],
)
def test_send_slack_message_carefully_params_fail(
    slack_hook: str,
    message: str,
    err_msg: str,
    caplog: pytest.LogCaptureFixture,
    mock_send_slack_message: MagicMock,
) -> None:
    """Test sending a slack message ever so carefully, but with incorrect parameters."""
    send_slack_message_carefully(slack_hook, message)
    mock_send_slack_message.assert_not_called()
    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.WARNING
    assert caplog.records[-1].message == f"Cannot send slack message: {err_msg}"


@pytest.mark.parametrize("markdown", [True, False, None])
def test_send_slack_message_carefully_markdown_params(
    markdown: None | bool, mock_send_slack_message: MagicMock
) -> None:
    """Test that the markdown param is correctly passed on to send_slack_message.

    :param markdown: markdown parameter
    :type markdown: None | bool
    """
    if markdown is None:
        send_slack_message_carefully(SLACK_HOOK, TEST_MESSAGE)
    else:
        send_slack_message_carefully(SLACK_HOOK, TEST_MESSAGE, markdown)

    if markdown:
        mock_send_slack_message.assert_called_once_with(SLACK_HOOK, TEST_MESSAGE, True)  # noqa: FBT003
    else:
        mock_send_slack_message.assert_called_once_with(SLACK_HOOK, TEST_MESSAGE, False)  # noqa: FBT003


def test_send_slack_message_fail_error_oh_no(monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture) -> None:
    """Ensure that errors are caught and don't crash the whole goddamn ship of fools."""
    slack_mock = MagicMock(side_effect=ValueError("Oh no! An error!"))
    monkeypatch.setattr(core, "send_slack_message", slack_mock)
    send_slack_message_carefully(SLACK_HOOK, TEST_MESSAGE)
    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.ERROR
    assert caplog.records[-1].message == "Failed to send slack message"


# sync_configs
def test_sync_configs_mutates_dlt_config_in_place(dlt_config: dict[str, Any], test_cts_settings: CtsSettings) -> None:
    """sync_configs mutates the supplied dlt_config dict in-place."""
    original_id = id(dlt_config)
    sync_configs(test_cts_settings, dlt_config)
    assert id(dlt_config) == original_id


def test_sync_configs_with_mock_dlt_config_object(test_cts_settings: CtsSettings) -> None:
    """sync_configs works with any mapping that supports __setitem__."""
    mock_cfg = MagicMock()
    sync_configs(test_cts_settings, mock_cfg)
    mock_cfg.__setitem__.assert_any_call("normalize.data_writer.disable_compression", test_cts_settings.dev_mode)
    mock_cfg.__setitem__.assert_any_call(
        f"destination.{test_cts_settings.use_destination}.bucket_url",
        test_cts_settings.output_dir,
    )


@pytest.mark.parametrize("use_destination", VALID_DESTINATIONS)
@pytest.mark.parametrize("dev_mode", [True, False])
@pytest.mark.parametrize("output_dir", ["/some/path", "s3://bucket/whatever"])
def test_sync_configs_both_keys_set_in_single_call(
    dlt_config: dict[str, Any],
    test_cts_settings: CtsSettings,
    dev_mode: bool,
    use_destination: str,
    output_dir: str,
) -> None:
    """Test that sync_configs changes the disable_compression and bucket_url values."""
    original_dlt_config = deepcopy(dlt_config)
    test_cts_settings.dev_mode = dev_mode
    test_cts_settings.use_destination = use_destination
    test_cts_settings.output_dir = output_dir
    sync_configs(test_cts_settings, dlt_config)
    assert dlt_config == {
        **original_dlt_config,
        "normalize.data_writer.disable_compression": dev_mode,
        f"destination.{use_destination}.bucket_url": output_dir,
    }


def test_sync_configs_attrs_missing_as_expected(
    settings_missing_sync_attrs: LoggerSettings | InputOutputSettings,
) -> None:
    """Sanity-check the fixture: confirm which attributes are genuinely absent."""
    assert hasattr(settings_missing_sync_attrs, OUTPUT_DIR) == isinstance(
        settings_missing_sync_attrs, InputOutputSettings
    )
    assert not hasattr(settings_missing_sync_attrs, DEV_MODE)
    assert not hasattr(settings_missing_sync_attrs, USE_DESTINATION)


def test_sync_configs_no_op_when_relevant_attrs_missing(
    settings_missing_sync_attrs: LoggerSettings | InputOutputSettings, empty_dlt_config: dict[str, Any]
) -> None:
    """sync_configs must not write any keys when the settings object lacks the attrs it needs.

    This covers both the case where dev_mode/output/use_destination are all absent, and the case
    where output_dir is present but use_destination is not (so the bucket_url key must still be skipped).
    """
    sync_configs(settings_missing_sync_attrs, empty_dlt_config)  # pyright: ignore[reportArgumentType]
    assert empty_dlt_config == {}


def test_sync_configs_no_op_with_mock_config_when_relevant_attrs_missing(
    settings_missing_sync_attrs: LoggerSettings | InputOutputSettings,
) -> None:
    """As above, but verified via a MagicMock config to confirm __setitem__ is never even attempted."""
    mock_cfg = MagicMock()
    sync_configs(settings_missing_sync_attrs, mock_cfg)  # pyright: ignore[reportArgumentType]
    mock_cfg.__setitem__.assert_not_called()


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_sync_configs_sets_both_keys_from_a_de_novo_dlt_config(
    settings_cls: type[CtsSettings], empty_dlt_config: dict[str, Any]
) -> None:
    """When dev_mode/output_dir/use_destination are all present, sync_configs sets exactly those two keys.

    ``empty_dlt_config`` starts out completely empty, distinct from the pre-populated ``dlt_config``
    fixture used (via the autouse patch) to validate ``settings_cls`` on construction. This confirms
    sync_configs's own behaviour depends only on its arguments, not on any leftover config state.
    """
    settings = settings_cls(dev_mode=True, output_dir="/some/output", use_destination="local_fs")  # pyright: ignore[reportCallIssue]

    sync_configs(settings, empty_dlt_config)

    assert empty_dlt_config == {
        "normalize.data_writer.disable_compression": True,
        "destination.local_fs.bucket_url": "/some/output",
    }


# tests for run_cli()
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_run_cli_calls_settings_cls_with_dlt_config(
    settings_cls: type[CtsSettings], dlt_config: dict[str, Any]
) -> None:
    """Ensure run_cli instantiates the supplied settings class with dlt.config."""
    captured: list[CtsSettings] = []

    class _CaptureCls(settings_cls):  # type: ignore[valid-type]
        def __init__(self, **data: Any) -> None:  # noqa: ANN401
            super().__init__(**data)
            assert data == {}
            captured.append(self)

    run_cli(_CaptureCls, MagicMock())

    assert len(captured) == 1
    captured_config = captured[0]
    assert isinstance(captured_config, settings_cls)
    # The object passed to sync_configs must be a fully initialised instance
    for attr in settings_cls.model_fields:
        assert hasattr(captured_config, attr)
    assert captured_config._dlt_config == dlt_config  # noqa: SLF001


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_run_cli_function_calls_args(settings_cls: type[CtsSettings]) -> None:
    """Ensure run_cli calls sync_configs with the instantiated config and dlt.config."""
    instantiated_cls = settings_cls()  # pyright: ignore[reportCallIssue]
    pipeline_fn_mock = MagicMock()
    settings_cls_mock = MagicMock(return_value=instantiated_cls)

    with (
        patch("cdm_data_loaders.pipelines.core.construct_env_var") as mock_env_var,
        patch("cdm_data_loaders.pipelines.core.sync_configs") as mock_sync,
    ):
        run_cli(settings_cls_mock, pipeline_fn_mock)  # type: ignore[reportArgumentType]

    mock_env_var.assert_called_once_with()
    settings_cls_mock.assert_called_once_with()
    mock_sync.assert_called_once_with(instantiated_cls, dlt.config)
    pipeline_fn_mock.assert_called_once_with(instantiated_cls)


# error handling: SettingsError/ValidationError/ValueError
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_run_cli_reraises_settings_error(
    settings_cls: type[CtsSettings],
    caplog: pytest.LogCaptureFixture,
) -> None:
    """SettingsError is printed and re-raised."""
    err = SettingsError("bad CLI arg")

    with (
        patch.object(settings_cls, "__init__", side_effect=err),
        pytest.raises(SettingsError, match="bad CLI arg"),
    ):
        run_cli(settings_cls, MagicMock())

    log_records = caplog.records
    assert log_records[-1].levelno == logging.ERROR
    assert log_records[-1].message == "Error initialising config"


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize(
    ("bad_dlt_config", "error", "err_msg"),
    [
        (None, ValidationError, "dlt_config must be defined"),
        ({}, ValueError, "No valid destinations found in dlt configuration"),
        ({"destination": {}}, ValueError, "No valid destinations found in dlt configuration"),
    ],
)
def test_run_cli_reraises_validation_errors(
    settings_cls: type[CtsSettings],
    bad_dlt_config: None | dict[str, Any],
    error: type[Exception],
    err_msg: str,
    caplog: pytest.LogCaptureFixture,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure that errors in instantiating the configuration are re-raised.

    See also the cts_defaults test ``test_cli_app_run_dlt_config_errors``.
    """
    monkeypatch.setattr(core.dlt, "config", bad_dlt_config)
    with pytest.raises(error, match=err_msg):
        run_cli(settings_cls, MagicMock())

    log_records = caplog.records
    assert log_records[-1].levelno == logging.ERROR
    assert log_records[-1].message == "Error initialising config"


# error handling: unexpected Exception
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_run_cli_reraises_unexpected_exception(
    settings_cls: type[CtsSettings], caplog: pytest.LogCaptureFixture
) -> None:
    """Ensure that other exceptions are caught and re-raised."""
    boom = RuntimeError("disk on fire")
    mock_pipeline_fn = MagicMock()
    with (
        patch.object(settings_cls, "__init__", side_effect=boom),
        patch("cdm_data_loaders.pipelines.core.sync_configs") as mock_sync,
        pytest.raises(RuntimeError, match="disk on fire"),
    ):
        run_cli(settings_cls, mock_pipeline_fn)

    log_records = caplog.records
    assert log_records[-1].levelno == logging.ERROR
    assert log_records[-1].message == "Unexpected error setting up config"

    mock_sync.assert_not_called()
    mock_pipeline_fn.assert_not_called()


# sync_configs not called on error
@pytest.mark.parametrize(
    "exc",
    [
        SettingsError("bad"),
        ValueError("bad value"),
    ],
)
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
def test_run_cli_sync_configs_not_called_on_settings_instantiation_error(
    exc: Exception,
    settings_cls: type[CtsSettings],
) -> None:
    """Ensure that further execution is stopped if settings instantiation fails."""
    mock_pipeline_fn = MagicMock()
    with (
        patch.object(settings_cls, "__init__", side_effect=exc),
        patch("cdm_data_loaders.pipelines.core.sync_configs") as mock_sync,
        pytest.raises(type(exc)),
    ):
        run_cli(settings_cls, mock_pipeline_fn)

    mock_sync.assert_not_called()
    mock_sync.assert_not_called()
    mock_pipeline_fn.assert_not_called()


def test_run_cli_uses_slack_env_var_if_set(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """RUNTIME__SLACK_INCOMING_HOOK built by construct_env_var is visible during run_cli."""
    monkeypatch.setenv("VARIABLE_B", "BBB")
    monkeypatch.setenv("VARIABLE_T", "TTT")
    monkeypatch.setenv("CHAR_STR", "CCC")

    run_cli(CtsSettings, MagicMock())

    expected = "https://hooks.slack.com/services/BBB/TTT/CCC/"
    assert os.environ.get("RUNTIME__SLACK_INCOMING_HOOK") == expected


def test_run_cli_no_slack_env_var_when_vars_missing(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """RUNTIME__SLACK_INCOMING_HOOK is not set when source vars are absent."""
    for var in ("VARIABLE_B", "VARIABLE_T", "CHAR_STR", "RUNTIME__SLACK_INCOMING_HOOK"):
        monkeypatch.delenv(var, raising=False)

    run_cli(CtsSettings, MagicMock())

    assert "RUNTIME__SLACK_INCOMING_HOOK" not in os.environ


# run_cli + sync_configs interaction: hasattr guards
SETTINGS_CLASSES_MISSING_SYNC_ATTRS: Final[list[type[LoggerSettings]]] = [LoggerSettings, InputOutputSettings]


def test_settings_classes_missing_sync_attrs_sanity_check() -> None:
    """Sanity-check that LoggerSettings/InputOutputSettings genuinely lack the attrs sync_configs checks.

    LoggerSettings has none of dev_mode/output/use_destination; InputOutputSettings has output but not
    dev_mode or use_destination, so the bucket_url branch (which requires both) should still be skipped.
    """
    logger_settings = LoggerSettings()  # pyright: ignore[reportCallIssue]
    io_settings = InputOutputSettings()  # pyright: ignore[reportCallIssue]

    assert not hasattr(logger_settings, DEV_MODE)
    assert not hasattr(logger_settings, OUTPUT_DIR)
    assert not hasattr(logger_settings, USE_DESTINATION)

    assert hasattr(io_settings, OUTPUT_DIR)
    assert not hasattr(io_settings, DEV_MODE)
    assert not hasattr(io_settings, USE_DESTINATION)


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES_MISSING_SYNC_ATTRS)
def test_run_cli_calls_sync_configs_even_when_attrs_missing(settings_cls: type[LoggerSettings]) -> None:
    """run_cli still calls sync_configs for settings classes lacking dev_mode/use_destination.

    Mirrors test_run_cli_function_calls_args: confirms sync_configs is invoked with the instantiated
    settings object and dlt.config regardless of which attributes that object actually has - the
    hasattr guards live inside sync_configs itself, not in run_cli, so the call should always happen.
    """
    instantiated_cls = settings_cls()  # pyright: ignore[reportCallIssue]
    pipeline_fn_mock = MagicMock()
    settings_cls_mock = MagicMock(return_value=instantiated_cls)

    with (
        patch("cdm_data_loaders.pipelines.core.construct_env_var") as mock_env_var,
        patch("cdm_data_loaders.pipelines.core.sync_configs") as mock_sync,
    ):
        run_cli(settings_cls_mock, pipeline_fn_mock)  # type: ignore[reportArgumentType]

    mock_env_var.assert_called_once_with()
    mock_sync.assert_called_once_with(instantiated_cls, dlt.config)
    pipeline_fn_mock.assert_called_once_with(instantiated_cls)


@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES_MISSING_SYNC_ATTRS)
def test_run_cli_sync_configs_no_op_on_de_novo_config_when_attrs_missing(
    settings_cls: type[LoggerSettings], monkeypatch: pytest.MonkeyPatch
) -> None:
    """End-to-end: running the full run_cli path against a de novo dlt.config leaves it untouched.

    dlt.config is replaced here with a brand new, empty dict - as if no environment variables or
    config files had ever contributed to it - confirming that sync_configs's hasattr guards correctly
    prevent any spurious key creation (or AttributeError) when exercised through the real CLI entry
    point, for settings classes that don't carry all of dev_mode/output_dir/use_destination.
    """
    empty_dlt_config: dict[str, Any] = {}
    monkeypatch.setattr(core.dlt, "config", empty_dlt_config)
    pipeline_fn_mock = MagicMock()

    run_cli(settings_cls, pipeline_fn_mock)  # type: ignore[reportArgumentType]

    assert empty_dlt_config == {}
    pipeline_fn_mock.assert_called_once()
    called_settings = pipeline_fn_mock.call_args.args[0]
    assert isinstance(called_settings, settings_cls)


# dlt.config state after successful run
@pytest.mark.parametrize("settings_cls", SETTINGS_CLASSES)
@pytest.mark.parametrize("dev_mode", [True, False])
@pytest.mark.parametrize(
    ("use_destination", "output_dir"), [("local_fs", "/some/path"), ("s3", "s3://bucket/whatever")]
)
def test_run_cli_dlt_config_updated_after_success(
    dlt_config: dict[str, Any],
    settings_cls: type[CtsSettings],
    dev_mode: bool,
    use_destination: str,
    output_dir: str,
) -> None:
    """Test that sync_configs changes the disable_compression and bucket_url values."""
    original_dlt_config = deepcopy(dlt_config)
    settings = settings_cls(dev_mode=dev_mode, output_dir=output_dir, use_destination=use_destination)  # pyright: ignore[reportCallIssue]
    assert settings.dev_mode == dev_mode
    assert settings.output_dir == output_dir
    assert settings.use_destination == use_destination
    settings_cls_mock = MagicMock(return_value=settings)

    run_cli(settings_cls_mock, MagicMock())  # type: ignore[reportArgumentType]

    assert dlt_config == {
        **original_dlt_config,
        "normalize.data_writer.disable_compression": dev_mode,
        f"destination.{use_destination}.bucket_url": output_dir,
    }


# run_pipeline tests
def test_run_pipeline_minimal(test_bfi_settings: BatchedFileInputSettings, mock_dlt: MagicMock) -> None:
    """Ensure pipeline.run is called with correct args in the simplest case."""
    fake_resource = MagicMock()
    run_pipeline(test_bfi_settings, fake_resource)
    assert_pipeline_run_correctly(mock_dlt, fake_resource, test_bfi_settings.use_destination, {}, {}, {})


@pytest.mark.parametrize("destination_kwargs", [None, {}, {"max_table_nesting": 0}])
@pytest.mark.parametrize("pipeline_kwargs", [None, {}, {"pipeline_name": "p", "dataset_name": "d"}])
@pytest.mark.parametrize("pipeline_run_kwargs", [None, {}, {"table_format": "delta"}])
def test_run_pipeline_destination_pipeline_pipeline_run_kwargs_set(
    mock_dlt: MagicMock,
    destination_kwargs: dict[str, Any] | None,
    pipeline_kwargs: dict[str, Any] | None,
    pipeline_run_kwargs: dict[str, Any] | None,
) -> None:
    """Ensure a non-empty output_dir sets the correct dlt.config bucket_url key."""
    settings = make_batched_settings(input_dir="/i", output_dir="/custom/output", use_destination="local_fs")
    fake_resource = MagicMock()
    run_pipeline(
        settings,
        fake_resource,
        destination_kwargs=destination_kwargs,
        pipeline_kwargs=pipeline_kwargs,
        pipeline_run_kwargs=pipeline_run_kwargs,
    )
    assert_pipeline_run_correctly(
        mock_dlt,
        fake_resource,
        "local_fs",
        destination_kwargs,
        pipeline_kwargs,
        pipeline_run_kwargs,
    )


def test_run_pipeline_graceful_fail(
    test_bfi_settings: BatchedFileInputSettings,
    mock_dlt: MagicMock,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Ensure that errors during pipeline runs are caught and do not cause the entire pipeline to go ka-boom."""
    error = RuntimeError("Oh crap!!")
    fake_resource = MagicMock()

    mock_dlt.pipeline.return_value.run.side_effect = error

    output = run_pipeline(test_bfi_settings, fake_resource)
    assert output is None

    assert caplog.records[-1].levelno == logging.ERROR
    assert caplog.records[-1].message.startswith("Pipeline failed: ")

    for m in caplog.records:
        assert not m.message.startswith("Work complete")


@pytest.mark.parametrize("slack_configured", [True, False])
@pytest.mark.parametrize("success", [True, False])
def test_run_pipeline_slack_configured(
    test_bfi_settings: BatchedFileInputSettings,
    mock_dlt: MagicMock,
    mock_send_slack_message: MagicMock,
    slack_configured: bool,
    success: bool,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test that a slack message is sent if the slack_incoming_hook runtime config value is available."""
    error = None
    # set up slack config
    if slack_configured:
        slack_hook = "http://some.url.slack.com"
        mock_dlt.pipeline.return_value.runtime_config.slack_incoming_hook = slack_hook

    # set an error to be triggered if success is false
    if not success:
        error = RuntimeError("Oh crap!!")
        mock_dlt.pipeline.return_value.run.side_effect = error

    run_pipeline(test_bfi_settings, MagicMock())

    if slack_configured:
        if success:
            mock_send_slack_message.assert_called_once_with(
                "http://some.url.slack.com",
                "Pipeline completed successfully!",
                False,  # noqa: FBT003
            )
        else:
            mock_send_slack_message.assert_called_once_with(
                "http://some.url.slack.com",
                f"Pipeline failed: {error!s}",
                False,  # noqa: FBT003
            )
        assert f"No Slack alerts will be sent: {WEBHOOK_NOT_CONFIGURED}" not in caplog.messages
    else:
        mock_send_slack_message.assert_not_called()
        assert f"No Slack alerts will be sent: {WEBHOOK_NOT_CONFIGURED}" in caplog.messages

    if success:
        # log messages on success
        assert caplog.records[-1].levelno == logging.INFO
        assert caplog.records[-1].message.startswith("Work complete!")
    else:
        assert caplog.records[-1].levelno == logging.ERROR
        assert caplog.records[-1].message.startswith("Pipeline failed: ")


def test_run_pipeline_sets_pipelines_dir_when_pipeline_dir_set(mock_dlt: MagicMock) -> None:
    """pipelines_dir is injected into pipeline_kwargs when config.pipeline_dir is set."""
    settings = make_batched_settings(input_dir="/i", output_dir="/out", use_output_dir_for_pipeline_metadata=True)
    assert settings.pipeline_dir is not None

    run_pipeline(settings, MagicMock())

    pipeline_call_kwargs = mock_dlt.pipeline.call_args.kwargs
    assert pipeline_call_kwargs["pipelines_dir"] == settings.pipeline_dir


def test_run_pipeline_no_pipelines_dir_when_pipeline_dir_none(
    test_bfi_settings: BatchedFileInputSettings, mock_dlt: MagicMock
) -> None:
    """pipelines_dir is absent from pipeline_kwargs when config.pipeline_dir is None."""
    assert test_bfi_settings.pipeline_dir is None
    run_pipeline(test_bfi_settings, MagicMock())
    pipeline_call_kwargs = mock_dlt.pipeline.call_args.kwargs
    assert "pipelines_dir" not in pipeline_call_kwargs


def test_run_pipeline_sets_dev_mode_in_pipeline_kwargs_when_true(mock_dlt: MagicMock) -> None:
    """dev_mode=True is forwarded to dlt.pipeline()."""
    settings = make_batched_settings(input_dir="/i", output_dir="/out", dev_mode=True)
    run_pipeline(settings, MagicMock())
    pipeline_call_kwargs = mock_dlt.pipeline.call_args.kwargs
    assert pipeline_call_kwargs.get("dev_mode") is True


def test_run_pipeline_dev_mode_absent_from_pipeline_kwargs_when_false(
    test_bfi_settings: BatchedFileInputSettings, mock_dlt: MagicMock
) -> None:
    """dev_mode=False is NOT forwarded to dlt.pipeline()."""
    assert test_bfi_settings.dev_mode is False
    run_pipeline(test_bfi_settings, MagicMock())
    pipeline_call_kwargs = mock_dlt.pipeline.call_args.kwargs
    assert "dev_mode" not in pipeline_call_kwargs
