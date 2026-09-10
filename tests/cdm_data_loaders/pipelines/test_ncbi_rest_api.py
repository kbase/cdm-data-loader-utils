"""Tests for the NCBI datasets API pipeline functions."""

import re
from pathlib import Path
from typing import Any
from unittest import mock
from unittest.mock import MagicMock, call, patch

import pytest
from dlt.extract.items import DataItemWithMeta
from frozendict import frozendict
from pydantic import ValidationError
from pydantic_settings import CliApp
from requests import HTTPError

from cdm_data_loaders.pipelines import ncbi_rest_api as ncbi_module
from cdm_data_loaders.pipelines.ncbi_rest_api import (
    ANNOTATION,
    DATASET,
    DATASET_NAME,
    ERROR,
    MAX_IDS_PER_QUERY,
    NcbiRestApiSettings,
    assemble_assembly_reports,
    assembly_list,
    cli,
    get_annotation_report,
    get_assembly_reports,
    get_dataset_reports,
    get_settings,
    run_ncbi_pipeline,
    set_settings,
)
from tests.cdm_data_loaders.core.conftest import (
    TEST_CTS_SETTINGS,
    TEST_CTS_SETTINGS_RECONCILED,
    check_settings,
    make_settings_autofill_config,
)
from tests.conftest import DEFAULT_VCR_CONFIG
from tests.helpers import assert_cli_field_roundtrips, assert_no_cli_clashes


@pytest.fixture(autouse=True)
def patch_rest_client_hooks(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that the REST_CLIENT_HOOKS dict is empty for tests."""
    monkeypatch.setattr("cdm_data_loaders.pipelines.ncbi_rest_api.REST_CLIENT_HOOKS", {})


ID_WITH_2K_ANNOTS = "GCF_000003135.1"
ID_WITH_500_ANNOTS = "GCF_000007725.1"
ID_TRIGGERS_500_ERR = "GCF_500_ERROR"
VALID_IDS = [ID_WITH_500_ANNOTS, ID_WITH_2K_ANNOTS]
INVALID_ID = "invalid_id"
ALL_IDS = [*VALID_IDS, INVALID_ID]
BATCH_SIZE = 200
BATCH_SIZE_STRING = str(BATCH_SIZE)

TEST_NCBI_SETTINGS = frozendict(**TEST_CTS_SETTINGS, batch_size=BATCH_SIZE_STRING, query_type="")

TEST_NCBI_SETTINGS_RECONCILED = frozendict(**TEST_CTS_SETTINGS_RECONCILED, batch_size=BATCH_SIZE, query_type=None)

TEST_NCBI_SETTINGS_V1 = frozendict(**TEST_CTS_SETTINGS, batch_size=BATCH_SIZE_STRING, query_type=" ANNOTATION ")

TEST_NCBI_SETTINGS_RECONCILED_V1 = frozendict(
    **TEST_CTS_SETTINGS_RECONCILED, batch_size=BATCH_SIZE, query_type=ANNOTATION
)


@pytest.fixture(scope="module")
def vcr_config() -> dict[str, Any]:
    """VCR config for tests that make HTTP requests."""
    return {**DEFAULT_VCR_CONFIG}


@pytest.fixture(scope="module")
def valid_assembly_ids() -> list[str]:
    """A list of assembly IDs."""
    return VALID_IDS


@pytest.fixture(scope="module")
def assembly_id(valid_assembly_ids: list[str]) -> str:
    """Single valid assembly ID."""
    return valid_assembly_ids[0]


@pytest.fixture(scope="module")
def invalid_assembly_id() -> str:
    """Invalid assembly ID."""
    return INVALID_ID


@pytest.fixture(scope="module")
def assembly_ids(valid_assembly_ids: list[str], invalid_assembly_id: str) -> list[str]:
    """List of assembly IDs including both valid and invalid IDs."""
    return [*valid_assembly_ids, invalid_assembly_id]


@pytest.fixture(scope="module")
def test_settings() -> NcbiRestApiSettings:
    """Generate a test settings class."""
    return make_settings_autofill_config(NcbiRestApiSettings)  # type: ignore[reportReturnType]


@pytest.fixture(autouse=True)
def reset_settings(monkeypatch: pytest.MonkeyPatch) -> None:
    """Reset the pipeline settings prior to running tests.

    :param monkeypatch: monkeypatch
    :type monkeypatch: pytest.MonkeyPatch
    """
    monkeypatch.setattr(ncbi_module, "PIPELINE_SETTINGS", None)


def test_get_set_settings(test_settings: NcbiRestApiSettings) -> None:
    """Test setting and retrieval of settings."""
    assert ncbi_module.PIPELINE_SETTINGS is None
    with pytest.raises(RuntimeError, match="Pipeline settings have not been initialised"):
        get_settings()

    set_settings(test_settings)
    assert test_settings == ncbi_module.PIPELINE_SETTINGS
    assert get_settings() == test_settings

    set_settings(None)  # type: ignore[reportArgumentType]
    assert ncbi_module.PIPELINE_SETTINGS is None
    with pytest.raises(RuntimeError, match="Pipeline settings have not been initialised"):
        get_settings()


@pytest.mark.parametrize(
    ("batch_size", "parsed_batch_size"),
    [
        ("10", 10),
        (10, 10),
        (None, MAX_IDS_PER_QUERY),
    ],
)
def test_settings_valid_batch_size(
    batch_size: int | str | None, parsed_batch_size: int, test_settings: NcbiRestApiSettings
) -> None:
    """Ensure that a valid batch size is correctly parsed."""
    if batch_size is None:
        settings = test_settings
    else:
        settings: NcbiRestApiSettings = make_settings_autofill_config(NcbiRestApiSettings, {"batch_size": batch_size})  # type: ignore[reportReturnType]
    assert settings.batch_size == parsed_batch_size


@pytest.mark.parametrize(
    ("bad_batch_size", "message"),
    [
        ("0", "Input should be greater than or equal to 1"),
        ("-1", "Input should be greater than or equal to 1"),
        ("1001", f"Input should be less than or equal to {MAX_IDS_PER_QUERY}"),
        ("notanint", "Input should be a valid integer"),
        ("", "Input should be a valid integer"),
        ("1.2345", "Input should be a valid integer"),
    ],
)
@pytest.mark.parametrize("use_cliapp", [True, False])
def test_cli_invalid_batch_size_via_cli_raises(bad_batch_size: str, message: str, use_cliapp: bool) -> None:
    """Ensure that an invalid batch size passed via CLI raises an error."""
    if use_cliapp:
        with pytest.raises(ValidationError, match=message):
            CliApp.run(NcbiRestApiSettings, cli_args=["--batch-size", bad_batch_size])
    else:
        with pytest.raises(ValidationError, match=message):
            make_settings_autofill_config(NcbiRestApiSettings, {"batch_size": bad_batch_size})


@pytest.mark.parametrize(
    ("query_type", "parsed_query_type"),
    [
        (ANNOTATION, ANNOTATION),
        (DATASET, DATASET),
        ("  ANNOTATION  ", ANNOTATION),
        ("\n\nDataset\t\n", DATASET),
        ("   ", None),
        ("", None),
        (None, None),
    ],
)
def test_cli_valid_query_type(
    query_type: str | None, parsed_query_type: str | None, test_settings: NcbiRestApiSettings
) -> None:
    """Ensure that a valid batch size is correctly parsed."""
    settings: NcbiRestApiSettings = make_settings_autofill_config(NcbiRestApiSettings, {"query_type": query_type})  # type: ignore[reportReturnType]
    assert settings.query_type == parsed_query_type

    if query_type is None:
        alt_settings: NcbiRestApiSettings = test_settings
        assert alt_settings.query_type == parsed_query_type


STRING_MATCH_MESSAGE = re.escape("String should match pattern '^(dataset|annotation)$'")


@pytest.mark.parametrize(
    "query_type",
    [
        "annot",
        "anotation",
        "data set",
    ],
)
@pytest.mark.parametrize("use_cliapp", [True, False])
def test_cli_invalid_query_type(query_type: str, use_cliapp: bool) -> None:
    """Ensure that an invalid batch size passed via CLI raises an error."""
    if use_cliapp:
        with pytest.raises(ValidationError, match=STRING_MATCH_MESSAGE):
            CliApp.run(NcbiRestApiSettings, cli_args=["--query-type", query_type])
    else:
        with pytest.raises(ValidationError, match=STRING_MATCH_MESSAGE):
            make_settings_autofill_config(NcbiRestApiSettings, {"query_type": query_type})


@pytest.mark.parametrize(
    ("settings", "reconciled"),
    [
        (TEST_NCBI_SETTINGS, TEST_NCBI_SETTINGS_RECONCILED),
        (TEST_NCBI_SETTINGS_V1, TEST_NCBI_SETTINGS_RECONCILED_V1),
    ],
)
def test_settings_all_params_set(settings: frozendict, reconciled: frozendict) -> None:
    """Ensure that settings are set correctly when all args are specified."""
    s = make_settings_autofill_config(NcbiRestApiSettings, settings)
    check_settings(s, reconciled)


# model_fields
def test_no_cli_param_clashes() -> None:
    """Ensure that the settings object does not have any internal alias collisions."""
    assert_no_cli_clashes(NcbiRestApiSettings)


@pytest.mark.parametrize(
    ("settings", "reconciled"),
    [
        (TEST_NCBI_SETTINGS, TEST_NCBI_SETTINGS_RECONCILED),
        (TEST_NCBI_SETTINGS_V1, TEST_NCBI_SETTINGS_RECONCILED_V1),
    ],
)
@pytest.mark.parametrize("field_name", list(TEST_NCBI_SETTINGS))
def test_cli_fields_parse_correctly(field_name: str, settings: frozendict, reconciled: frozendict) -> None:
    """Ensure that all fields and variants are parsed correctly."""
    assert_cli_field_roundtrips(
        NcbiRestApiSettings,
        field_name,
        settings[field_name],
        reconciled[field_name],
    )


def test_cli_passes_settings_class_to_run_cli() -> None:
    """Ensure that cli() calls run_cli with NcbiRestApiSettings as the settings class."""
    with patch.object(ncbi_module, "run_cli") as mock_run_cli:
        cli()

    mock_run_cli.assert_called_once()
    assert mock_run_cli.call_args[0] == (NcbiRestApiSettings, run_ncbi_pipeline)


def test_cli_calls_run_ncbi_pipeline(monkeypatch: pytest.MonkeyPatch, dlt_config: dict[str, Any]) -> None:
    """Ensure that cli() calls run_ncbi_pipeline with the settings."""
    mock_settings_instance = MagicMock()
    mock_settings_cls = MagicMock(return_value=mock_settings_instance)
    mock_run_ncbi_pipeline = MagicMock()

    monkeypatch.setattr(ncbi_module, "NcbiRestApiSettings", mock_settings_cls)
    monkeypatch.setattr(ncbi_module, "run_ncbi_pipeline", mock_run_ncbi_pipeline)

    cli()

    mock_settings_cls.assert_called_once_with()
    mock_run_ncbi_pipeline.assert_called_once_with(mock_settings_instance)


def check_dataset_report(dataset_report: dict[str, Any] | None, assembly_id: str) -> None:
    """Check the basic structure of a dataset report."""
    assert dataset_report is not None
    assert dataset_report["accession"] == assembly_id
    for key in [
        "current_accession",
        "source_database",
        "organism",
        "assembly_info",
        "assembly_stats",
    ]:
        assert key in dataset_report


def check_annotation_report(annotation_report: list[dict[str, Any]] | None, assembly_id: str) -> None:
    """Check the basic structure of an annotation report."""
    assert annotation_report is not None
    for item in annotation_report:
        assert isinstance(item, dict)
        assert "row_id" in item
        assert "annotation" in item
        assert item.get("annotation", {}).get("annotations", [{}])[0].get("assembly_accession") == assembly_id
    all_row_ids = [int(item["row_id"]) for item in annotation_report]
    assert all_row_ids == list(range(1, len(all_row_ids) + 1))


def test_assembly_list_resource() -> None:
    """Test that the assembly list resource yields the expected assembly IDs."""
    settings: NcbiRestApiSettings = make_settings_autofill_config(
        NcbiRestApiSettings, {"input_dir": "tests/data/ncbi_rest_api/input"}
    )  # type: ignore[reportAssignmentType]
    set_settings(settings)

    ass_list = list(assembly_list())
    assert ass_list == [
        "GCF_029958545.3",
        "GCF_029958565.3",
        "GCF_029958585.3",
        "invalid_id",
        "GCF_029958605.3",
        "GCF_029958625.3",
        "GCF_029958645.3",
        "GCF_029958665.3",
    ]


@pytest.mark.parametrize("dev_mode", [False, True, None])
@pytest.mark.parametrize("use_pipeline_dir", [False, True, None])
def test_run_ncbi_pipeline_sets_core_run_pipeline_args_correctly(
    dev_mode: bool | None,
    use_pipeline_dir: bool | None,
    mock_dlt: MagicMock,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Ensure that run_ncbi_pipeline calls core.run_pipeline with the correct args."""
    mock_assembly_report_parser = MagicMock()
    monkeypatch.setattr(ncbi_module, "assembly_report_parser", mock_assembly_report_parser)
    mock_assembly_list = MagicMock()
    monkeypatch.setattr(ncbi_module, "assembly_list", mock_assembly_list)

    base_settings: dict[str, str | bool] = {
        "input_dir": "tests/data/ncbi_rest_api/input",
        "output_dir": "/some/dir",
    }
    if dev_mode is not None:
        base_settings["dev_mode"] = dev_mode
    if use_pipeline_dir is not None:
        base_settings["use_output_dir_for_pipeline_metadata"] = use_pipeline_dir

    settings: NcbiRestApiSettings = make_settings_autofill_config(NcbiRestApiSettings, base_settings)  # type: ignore[reportAssignmentType]

    check_settings(
        settings,
        {
            "dev_mode": bool(dev_mode),
            "input_dir": "tests/data/ncbi_rest_api/input",
            "log_config_file": None,
            "output_dir": "/some/dir",
            "pipeline_dir": "/some/dir/.dlt_conf" if use_pipeline_dir else None,
            "raw_data_dir": "/some/dir/raw_data",
            "use_destination": "local_fs",
            "use_output_dir_for_pipeline_metadata": bool(use_pipeline_dir),
            "batch_size": MAX_IDS_PER_QUERY,
            "query_type": None,
        },
    )

    run_ncbi_pipeline(settings)

    mock_dlt.destination.assert_called_once_with(settings.use_destination, max_table_nesting=0)
    mock_dlt.destination.assert_called_once()
    assert mock_dlt.destination.call_args_list[0].kwargs == {"max_table_nesting": 0}
    assert mock_dlt.destination.call_args_list[0].args == ("local_fs",)

    mock_dlt.pipeline.assert_called_once()
    assert mock_dlt.pipeline.call_args.kwargs["destination"] == mock_dlt.destination.return_value
    assert mock_dlt.pipeline.call_args.kwargs["pipeline_name"] == DATASET_NAME
    assert mock_dlt.pipeline.call_args.kwargs["dataset_name"] == DATASET_NAME
    if dev_mode:  # truthy
        assert mock_dlt.pipeline.call_args.kwargs["dev_mode"] is True
    else:
        assert "dev_mode" not in mock_dlt.pipeline.call_args.kwargs
    if use_pipeline_dir:  # truthy
        assert mock_dlt.pipeline.call_args.kwargs["pipelines_dir"] == f"{settings.output_dir}/.dlt_conf"  # type: ignore[reportArgumentType]
    else:
        assert "pipelines_dir" not in mock_dlt.pipeline.call_args.kwargs

    mock_dlt.pipeline.return_value.run.assert_called_once_with([mock_assembly_report_parser])


@pytest.mark.default_cassette("test_get_assembly_reports.yaml")
@pytest.mark.vcr
def test_get_dataset_reports() -> None:
    """Ensure that every assembly ID appears as a key in the output dict with the appropriate dict output."""
    dataset_report = get_dataset_reports(ALL_IDS)
    assert set(dataset_report.keys()) == set(ALL_IDS)
    assert dataset_report[INVALID_ID] is None
    for assembly_id in VALID_IDS:
        check_dataset_report(dataset_report.get(assembly_id), assembly_id)


def test_get_dataset_reports_empty_id_list_yields_empty_dict() -> None:
    """An empty input list produces an empty output dict."""
    assert get_dataset_reports([]) == {}


@pytest.mark.default_cassette("test_get_assembly_reports.yaml")
@pytest.mark.vcr
def test_get_annotation_report_single_page() -> None:
    """Test the retrieval of an annotation report with a single page."""
    annotation_report = get_annotation_report(ID_WITH_500_ANNOTS)
    check_annotation_report(annotation_report, ID_WITH_500_ANNOTS)


@pytest.mark.default_cassette("test_get_assembly_reports.yaml")
@pytest.mark.vcr
def test_get_annotation_report_multi_page() -> None:
    """Test the retrieval of an annotation report with multiple pages."""
    annotation_report = get_annotation_report(ID_WITH_2K_ANNOTS)
    assert isinstance(annotation_report, list)
    check_annotation_report(annotation_report, ID_WITH_2K_ANNOTS)


@mock.patch("tenacity.nap.time.sleep", MagicMock())
@pytest.mark.vcr
def test_get_annotation_report_multi_page_err() -> None:
    """An error in the middle of a multi-page retrieval should stop the whole retrieval process."""
    with pytest.raises(HTTPError, match="500 Server Error: Internal Server Error for url"):
        get_annotation_report(ID_TRIGGERS_500_ERR)


@pytest.mark.default_cassette("test_get_assembly_reports.yaml")
@pytest.mark.vcr
def test_get_annotation_report_invalid_id() -> None:
    """Test the retrieval of an annotation report for an invalid ID."""
    assert get_annotation_report(INVALID_ID) is None


def test_get_assembly_reports_empty_id_list(test_settings: NcbiRestApiSettings) -> None:
    """Ensure that getting reports for an empty list returns nothing."""
    set_settings(test_settings)
    assert get_assembly_reports([]) == {}


@pytest.mark.vcr
def test_get_assembly_reports(test_settings: NcbiRestApiSettings) -> None:
    """Test the retrieval of annotation and dataset reports."""
    set_settings(test_settings)
    assembly_reports = get_assembly_reports(ALL_IDS)
    assert set(assembly_reports) == {DATASET, ANNOTATION, ERROR}
    for datatype in [DATASET, ANNOTATION]:
        assert set(assembly_reports[datatype]) == set(ALL_IDS)
        assert assembly_reports[datatype][INVALID_ID] is None
    for assembly_id in [ID_WITH_2K_ANNOTS, ID_WITH_500_ANNOTS]:
        check_annotation_report(assembly_reports[ANNOTATION][assembly_id], assembly_id)
        check_dataset_report(assembly_reports[DATASET][assembly_id], assembly_id)
    assert assembly_reports[ERROR] == []


@pytest.mark.parametrize("query_type", [None, DATASET, ANNOTATION])
def test_get_assembly_reports_mock_subs(query_type: str | None, monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that the correct subs are called and the correct subs are not called with the query_type parameter."""
    settings: NcbiRestApiSettings = make_settings_autofill_config(NcbiRestApiSettings, {"query_type": query_type})  # type: ignore[reportAssignmentType]
    set_settings(settings)
    mock_get_annotation_report = MagicMock(return_value={"this": "that"})
    mock_get_dataset_reports = MagicMock(return_value=dict.fromkeys(ALL_IDS, "blob"))

    monkeypatch.setattr(ncbi_module, "get_annotation_report", mock_get_annotation_report)
    monkeypatch.setattr(ncbi_module, "get_dataset_reports", mock_get_dataset_reports)

    output = get_assembly_reports(ALL_IDS)
    assert output[ERROR] == []

    if query_type == ANNOTATION:
        mock_get_dataset_reports.assert_not_called()
        assert DATASET not in output

    else:
        mock_get_dataset_reports.assert_called_once_with(ALL_IDS)
        assert output[DATASET] == dict.fromkeys(ALL_IDS, "blob")

    if query_type == DATASET:
        mock_get_annotation_report.assert_not_called()
        assert ANNOTATION not in output
    else:
        assert mock_get_annotation_report.call_args_list == [
            call(
                n,
            )
            for n in ALL_IDS
        ]
        assert output[ANNOTATION] == {this_id: {"this": "that"} for this_id in ALL_IDS}


RECORDED_ERRORS = {
    "dataset_404": {
        "assembly_id": None,
        "assembly_id_list": ALL_IDS,
        "error_class": "HTTPError",
        "error_from": "dataset_report",
        "message": '404 Client Error: Not Found for url: https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000007725.1%2CGCF_000003135.1%2Cinvalid_id/dataset_report?page_size=1000\nResponse: {"error":"Not Found","code":404,"message":"Your request is invalid. (For more help, see the NCBI Datasets Documentation at https://www.ncbi.nlm.nih.gov/datasets/docs/)"}\n',
        "request_url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000007725.1%2CGCF_000003135.1%2Cinvalid_id/dataset_report?page_size=1000",
        "status": 404,
        "reason": "Not Found",
    },
    "annotation_report_500": {
        "assembly_id": ID_WITH_2K_ANNOTS,
        "assembly_id_list": None,
        "error_class": "HTTPError",
        "error_from": "annotation_report",
        "message": '500 Server Error: Internal Server Error for url: https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000003135.1/annotation_report?page_size=1000&page_token=eNrjYos2NDAwjAUABagBiw\nResponse: {"error":"Internal Server Error","code":500,"message":"Internal Server Error (For more help, see the NCBI Datasets Documentation at https://www.ncbi.nlm.nih.gov/datasets/docs/)"}\n',
        "request_url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000003135.1/annotation_report?page_size=1000&page_token=eNrjYos2NDAwjAUABagBiw",
        "status": 500,
        "reason": "Internal Server Error",
    },
    "annotation_report_404": {
        "assembly_id": ID_WITH_500_ANNOTS,
        "assembly_id_list": None,
        "error_class": "HTTPError",
        "error_from": "annotation_report",
        "message": '404 Client Error: Not Found for url: https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000007725.1/annotation_report?page_size=1000\nResponse: {"error":"Not Found","code":404,"message":"Your request is invalid. (For more help, see the NCBI Datasets Documentation at https://www.ncbi.nlm.nih.gov/datasets/docs/)"}\n',
        "request_url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_000007725.1/annotation_report?page_size=1000",
        "status": 404,
        "reason": "Not Found",
    },
    "value_error": {
        "assembly_id": INVALID_ID,
        "assembly_id_list": None,
        "error_class": "ValueError",
        "error_from": "annotation_report",
        "message": f"Some error message involving {INVALID_ID}.",
        "request_url": None,
        "status": None,
        "reason": None,
    },
}


@mock.patch("tenacity.nap.time.sleep", MagicMock())
@pytest.mark.default_cassette("test_get_assembly_reports_annotation_report_errors.yaml")
@pytest.mark.vcr
def test_get_assembly_reports_annotation_report_errors(test_settings: NcbiRestApiSettings) -> None:
    """Test the retrieval of assembly data when errors occur fetching annotation reports."""
    set_settings(test_settings)
    original_get_annotation_report = get_annotation_report

    def patched_get_annotation_report(assembly_id: str) -> list[dict[str, Any]] | None:
        """Patched version of get_annotation_report that throws a value error with a certain input.

        :param assembly_id: assembly ID
        :type assembly_id: str
        :raises ValueError: if the ID is INVALID_ID
        :return: output from the real get_annotation_report
        :rtype: list[dict[str, Any]] | None
        """
        if assembly_id == INVALID_ID:
            err_msg = f"Some error message involving {INVALID_ID}."
            raise ValueError(err_msg)
        return original_get_annotation_report(assembly_id)

    with mock.patch(
        "cdm_data_loaders.pipelines.ncbi_rest_api.get_annotation_report",
        side_effect=patched_get_annotation_report,
    ):
        assembly_reports = get_assembly_reports(ALL_IDS)

    assert set(assembly_reports) == {DATASET, ANNOTATION, ERROR}
    for datatype in [DATASET, ANNOTATION]:
        assert set(assembly_reports[datatype]) == set(ALL_IDS)
        assert assembly_reports[datatype][INVALID_ID] is None
    for assembly_id in [ID_WITH_2K_ANNOTS, ID_WITH_500_ANNOTS]:
        check_dataset_report(assembly_reports[DATASET][assembly_id], assembly_id)
    # ID_WITH_500 succeeds, ID_WITH_2K does not
    check_annotation_report(assembly_reports[ANNOTATION][ID_WITH_500_ANNOTS], ID_WITH_500_ANNOTS)
    assert assembly_reports[ANNOTATION][ID_WITH_2K_ANNOTS] is None

    assert assembly_reports[ERROR] == [
        RECORDED_ERRORS["annotation_report_500"],
        RECORDED_ERRORS["value_error"],
    ]


@mock.patch("tenacity.nap.time.sleep", MagicMock())
@pytest.mark.vcr
def test_get_assembly_reports_dataset_report_errors(test_settings: NcbiRestApiSettings) -> None:
    """Test the retrieval of assembly data when an error occurs fetching dataset reports."""
    set_settings(test_settings)
    assembly_reports = get_assembly_reports(ALL_IDS)
    assert set(assembly_reports) == {DATASET, ANNOTATION, ERROR}
    for datatype in [DATASET, ANNOTATION]:
        assert set(assembly_reports[datatype]) == set(ALL_IDS)
        assert assembly_reports[datatype][INVALID_ID] is None
    for assembly_id in [ID_WITH_2K_ANNOTS, ID_WITH_500_ANNOTS]:
        check_annotation_report(assembly_reports[ANNOTATION][assembly_id], assembly_id)
        assert assembly_reports[DATASET][assembly_id] is None

    assert assembly_reports[ERROR] == [RECORDED_ERRORS["dataset_404"]]


@mock.patch("tenacity.nap.time.sleep", MagicMock())
@pytest.mark.vcr
def test_get_assembly_reports_total_wipeout(test_settings: NcbiRestApiSettings) -> None:
    """Test the retrieval of assembly data when all queries fail."""
    set_settings(test_settings)
    original_get_annotation_report = get_annotation_report

    def patched_get_annotation_report(assembly_id: str) -> list[dict[str, Any]] | None:
        """Patched version of get_annotation_report that throws a value error with a certain input.

        :param assembly_id: assembly ID
        :type assembly_id: str
        :raises ValueError: if the ID is INVALID_ID
        :return: output from the real get_annotation_report
        :rtype: list[dict[str, Any]] | None
        """
        if assembly_id == INVALID_ID:
            err_msg = f"Some error message involving {INVALID_ID}."
            raise ValueError(err_msg)
        return original_get_annotation_report(assembly_id)

    with mock.patch(
        "cdm_data_loaders.pipelines.ncbi_rest_api.get_annotation_report",
        side_effect=patched_get_annotation_report,
    ):
        output = get_assembly_reports(ALL_IDS)

    assert output == {
        DATASET: dict.fromkeys(ALL_IDS),
        ANNOTATION: dict.fromkeys(ALL_IDS),
        ERROR: [
            RECORDED_ERRORS["dataset_404"],
            RECORDED_ERRORS["annotation_report_404"],
            RECORDED_ERRORS["annotation_report_500"],
            RECORDED_ERRORS["value_error"],
        ],
    }


@pytest.mark.skip("FIXME: not working, possibly due to parallelization?")
@pytest.mark.vcr
def test_get_assembly_report_parser_with_cassette(tmp_path: Path) -> None:
    with patch("dlt.mark"):
        settings: NcbiRestApiSettings = make_settings_autofill_config(
            NcbiRestApiSettings,
            {"input_dir": "tests/data/ncbi_rest_api/input", "output_dir": str(tmp_path)},
        )  # type: ignore[reportAssignmentType]
        run_ncbi_pipeline(settings)


def collect_results(reports: dict) -> dict[str, list]:
    """Drain the generator returned by assemble_assembly_reports into a dict keyed by table name."""
    results: dict[str, list] = {}
    for item in assemble_assembly_reports(reports):
        assert isinstance(item, DataItemWithMeta), f"Expected DataItemWithMeta, got {type(item)}"
        table_name = item.meta.table_name
        results.setdefault(table_name, [])
        results[table_name].extend(item.data)
    return results


# assemble_assembly_reports tests

DATASET_REPORT_1 = {
    "accession": ID_WITH_2K_ANNOTS,
    "organism": {"tax_id": 9606, "organism_name": "Homo sapiens"},
    "assembly_stats": {"total_sequence_length": 3099734149},
}

DATASET_REPORT_2 = {
    "accession": ID_WITH_500_ANNOTS,
    "organism": {"tax_id": 10090, "organism_name": "Mus musculus"},
    "assembly_stats": {"total_sequence_length": 2728222451},
}

ANNOTATION_REPORT_1 = [
    {"release_date": "2022-01-01", "annotation_name": "Annotation A"},
    {"release_date": "2022-06-01", "annotation_name": "Annotation B"},
]

ANNOTATION_REPORT_2 = [
    {"release_date": "2023-01-01", "annotation_name": "Annotation C"},
]

ERROR_REPORT = {
    "assembly_id": INVALID_ID,
    "assembly_id_list": None,
    "error_class": "HTTPError",
    "error_from": "dataset_report",
    "message": "404 Not Found",
    "request_url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCA_000001635.9/dataset_report",
    "status": 404,
    "reason": "Not Found",
}


@pytest.fixture
def full_assembly_reports() -> dict[str, Any]:
    """Assembly reports fixture with valid data for two assemblies and no errors."""
    return {
        DATASET: {
            ID_WITH_2K_ANNOTS: DATASET_REPORT_1,
            ID_WITH_500_ANNOTS: DATASET_REPORT_2,
        },
        ANNOTATION: {
            ID_WITH_2K_ANNOTS: ANNOTATION_REPORT_1,
            ID_WITH_500_ANNOTS: ANNOTATION_REPORT_2,
        },
        ERROR: [],
    }


@pytest.fixture
def reports_with_errors() -> dict[str, Any]:
    """Assembly reports fixture with one None dataset report, one None annotation report, and one error entry."""
    return {
        DATASET: {
            ID_WITH_2K_ANNOTS: DATASET_REPORT_1,
            ID_WITH_500_ANNOTS: DATASET_REPORT_2,
            INVALID_ID: None,
        },
        ANNOTATION: {
            ID_WITH_2K_ANNOTS: ANNOTATION_REPORT_1,
            ID_WITH_500_ANNOTS: ANNOTATION_REPORT_2,
            INVALID_ID: None,
        },
        ERROR: [ERROR_REPORT],
    }


EXPECTED_DB_TABLES = {
    f"{DATASET}_report": [
        {"assembly_id": ID_WITH_2K_ANNOTS, **DATASET_REPORT_1},
        {
            "assembly_id": ID_WITH_500_ANNOTS,
            **DATASET_REPORT_2,
        },
    ],
    f"{ANNOTATION}_report": [
        {"assembly_id": ID_WITH_2K_ANNOTS, **ANNOTATION_REPORT_1[0]},
        {"assembly_id": ID_WITH_2K_ANNOTS, **ANNOTATION_REPORT_1[1]},
        {"assembly_id": ID_WITH_500_ANNOTS, **ANNOTATION_REPORT_2[0]},
    ],
}

EXPECTED_DB_TABLES_WITH_ERROR = {
    f"{DATASET}_report": [*EXPECTED_DB_TABLES[f"{DATASET}_report"], {"assembly_id": INVALID_ID}],
    f"{ANNOTATION}_report": [*EXPECTED_DB_TABLES[f"{ANNOTATION}_report"]],
    "ncbi_import_error": [ERROR_REPORT],
}


@pytest.mark.parametrize("reports", [{}, None])
def test_assemble_assembly_reports_empty_dict_yields_nothing(
    reports: None | dict[str, Any],
) -> None:
    """Ensure that empty or None as input produces no output items."""
    assert list(assemble_assembly_reports(reports)) == []  # type: ignore[reportArgumentType]


def test_assemble_assembly_reports_yields_two_items_when_no_errors(
    full_assembly_reports: dict[str, Any],
) -> None:
    """When the error list is empty, the generator should yield only two DataItemWithMeta objects."""
    assert len(list(assemble_assembly_reports(full_assembly_reports))) == 2


def test_assemble_assembly_reports_all_items_are_data_item_with_meta(
    full_assembly_reports: dict[str, Any],
) -> None:
    """Every item yielded by the generator should be a DataItemWithMeta instance."""
    for item in assemble_assembly_reports(full_assembly_reports):
        assert isinstance(item, DataItemWithMeta)


def test_assemble_assembly_reports_table_names_when_errors_present(
    reports_with_errors: dict[str, Any],
) -> None:
    """The output should contain dataset_report, annotation_report, and ncbi_import_error table names."""
    results = collect_results(reports_with_errors)
    assert set(results) == {"dataset_report", "annotation_report", "ncbi_import_error"}
    assert results == EXPECTED_DB_TABLES_WITH_ERROR


def test_assemble_assembly_reports_no_error_table_when_error_list_empty(
    full_assembly_reports: dict[str, Any],
) -> None:
    """When the error list is empty, no ncbi_import_error table should be present in the output."""
    results = collect_results(full_assembly_reports)
    assert set(results) == {"dataset_report", "annotation_report"}
    assert results == EXPECTED_DB_TABLES


def test_assemble_assembly_reports_dataset_report_none_report_emits_row_with_only_assembly_id() -> None:
    """A None dataset report should still produce a row containing only the assembly_id key."""
    reports = {
        DATASET: {INVALID_ID: None},
        ANNOTATION: {INVALID_ID: None},
        ERROR: [],
    }
    results = collect_results(reports)
    # no ncbi_import_error or annotation_report
    assert results["dataset_report"] == [{"assembly_id": INVALID_ID}]


def test_assemble_assembly_reports_multiple_errors_all_yielded() -> None:
    """All entries in the error list should appear as individual rows in the ncbi_import_error table."""
    errors = [
        {"error_class": "HTTPError", "error_from": "dataset_report", "message": "err1"},
        {"error_class": "ValueError", "error_from": "annotation_report", "message": "err2"},
    ]
    reports = {
        DATASET: {},
        ANNOTATION: {},
        ERROR: errors,
    }
    results = collect_results(reports)
    assert results["ncbi_import_error"] == errors
    assert set(results) == {"ncbi_import_error"}
