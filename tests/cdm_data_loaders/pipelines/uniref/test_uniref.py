"""Tests for the UniRef DLT pipeline."""

from collections.abc import Callable
from unittest.mock import MagicMock, patch

import pytest

from cdm_data_loaders.parsers.uniprot.uniref import ENTRY_XML_TAG
from cdm_data_loaders.pipelines import uniref as uniref_module
from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.pipelines.uniref import (
    UnirefSettings,
    cli,
    parse_uniref,
    run_uniref_pipeline,
)


def test_cli_passes_settings_class_to_run_cli() -> None:
    """Ensure that cli is a PipelineCli wired to UnirefSettings and run_uniref_pipeline."""
    assert isinstance(cli, PipelineCli)
    assert cli.settings_cls is UnirefSettings
    assert cli.pipeline_fn is run_uniref_pipeline


# Tests for running the pipeline itself
def test_run_uniref_pipeline_args_set_correctly(test_settings: UnirefSettings) -> None:
    """Ensure that the pipeline arguments are set correctly, and each pipeline has a different name."""
    with patch.object(uniref_module, "run_pipeline") as mock_run_pipeline:
        run_uniref_pipeline(test_settings)

    assert mock_run_pipeline.call_count == 1
    _, kwargs = mock_run_pipeline.call_args
    assert kwargs.keys() == {"settings", "resource", "pipeline_kwargs"}
    assert kwargs["pipeline_kwargs"] == {
        "pipeline_name": "uniref",
        "dataset_name": f"uniref_{test_settings.variant}",
    }
    assert kwargs["settings"] == test_settings
    assert isinstance(kwargs["resource"], Callable)


def test_run_uniref_pipeline_sets_core_run_pipeline_args_correctly(
    test_settings: UnirefSettings, mock_dlt: MagicMock, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure that run_uniref_pipeline calls core.run_pipeline with the correct args."""
    mock_parse_uniref = MagicMock()
    monkeypatch.setattr(uniref_module, "parse_uniref", mock_parse_uniref)

    run_uniref_pipeline(test_settings)

    # parse_uniref was called once with the settings to produce the resource
    mock_parse_uniref.assert_called_once_with(test_settings)

    # the return value of parse_uniref(settings) is what gets passed to pipeline.run
    expected_resource = mock_parse_uniref.return_value

    mock_dlt.destination.assert_called_once_with(test_settings.use_destination, bucket_url=test_settings.output_dir)
    mock_dlt.pipeline.assert_called_once_with(
        destination=mock_dlt.destination.return_value,
        pipeline_name="uniref",
        dataset_name=f"uniref_{test_settings.variant}",
    )
    first_run, load_info_save_run = mock_dlt.pipeline.return_value.run.call_args_list
    assert first_run.args == (expected_resource,)
    assert first_run.kwargs == {}
    assert load_info_save_run.args[0] is not None
    assert load_info_save_run.kwargs == {"loader_file_format": "jsonl"}


def test_parse_uniref_resource(test_settings: UnirefSettings) -> None:
    """Ensure that parse_uniref calls build_xml_file_resource with the namespaced UniRef XML tag."""
    with patch.object(uniref_module, "build_xml_file_resource") as mock_build:
        resource = parse_uniref(test_settings)

    assert resource is mock_build.return_value
    assert mock_build.call_count == 1
    kwargs = mock_build.call_args.kwargs
    assert kwargs.keys() == {"settings", "xml_tag", "parse_fn", "resource_name"}
    assert kwargs["xml_tag"] == ENTRY_XML_TAG
    assert kwargs["settings"] == test_settings
    assert kwargs["resource_name"] == "parse_uniref"
    assert isinstance(kwargs["parse_fn"], Callable)
