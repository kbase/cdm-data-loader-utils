"""Tests for the all_the_bacteria pipeline."""

import csv
import logging
import re
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, call, patch

import pytest
from frozendict import frozendict
from requests.exceptions import HTTPError

from cdm_data_loaders.core.fields import VALID_DESTINATIONS
from cdm_data_loaders.pipelines import all_the_bacteria
from cdm_data_loaders.pipelines.all_the_bacteria import (
    ALL_ATB_FILE_NAME,
    DATASET_NAME,
    AtbSettings,
    cli,
    download_atb_index_tsv,
    get_file_download_links,
    load_patterns,
    osf_file_downloader,
    run_atb_pipeline,
)
from cdm_data_loaders.utils.file_transfer.core import NonRetryableDownloadError
from tests.cdm_data_loaders.core.conftest import (
    TEST_CTS_SETTINGS,
    TEST_CTS_SETTINGS_RECONCILED,
)
from tests.conftest import DEFAULT_VCR_CONFIG
from tests.helpers import assert_cli_field_roundtrips, assert_no_cli_clashes

TEST_SETTINGS = frozendict(**TEST_CTS_SETTINGS, version="1.2.3")

TEST_SETTINGS_RECONCILED = frozendict(**TEST_CTS_SETTINGS_RECONCILED, version="1.2.3", pattern_file=None)

TEST_SETTINGS_V1 = frozendict(**TEST_CTS_SETTINGS, version=1, pattern_file="some/path/or/other.txt")

TEST_SETTINGS_RECONCILED_V1 = frozendict(
    **TEST_CTS_SETTINGS_RECONCILED, version="1", pattern_file="some/path/or/other.txt"
)


@pytest.fixture(scope="module")
def vcr_config() -> dict[str, Any]:
    """VCR config for tests that make HTTP requests."""
    return {**DEFAULT_VCR_CONFIG}


@pytest.fixture
def test_settings(tmp_path: Path) -> AtbSettings:
    """Generate fake settings for testing."""
    return AtbSettings(output_dir=str(tmp_path))  # pyright: ignore[reportCallIssue]


@pytest.fixture
def test_s3_settings() -> AtbSettings:
    """Generate fake settings that use s3."""
    return AtbSettings(use_destination="s3")  # pyright: ignore[reportCallIssue]


@pytest.fixture
def pattern_file(tmp_path: Path) -> Path:
    """Pattern file for testing load_patterns."""
    p = tmp_path / "patterns.txt"
    p.write_text(
        "hello world\nfoo.bar\nstarts with*\n2+2=4\n",
        encoding="utf-8",
    )
    return p


# model_fields
def test_no_cli_param_clashes() -> None:
    """Ensure that the settings object does not have any internal alias collisions."""
    assert_no_cli_clashes(AtbSettings)


@pytest.mark.parametrize(
    ("settings", "reconciled"),
    [
        (TEST_SETTINGS, TEST_SETTINGS_RECONCILED),
        (TEST_SETTINGS_V1, TEST_SETTINGS_RECONCILED_V1),
    ],
)
@pytest.mark.parametrize("field_name", list(TEST_SETTINGS))
def test_cli_fields_parse_correctly(
    field_name: str, settings: frozendict[str, Any], reconciled: frozendict[str, Any]
) -> None:
    """Ensure that all fields and variants are parsed correctly."""
    assert_cli_field_roundtrips(
        AtbSettings,
        field_name,
        settings[field_name],
        reconciled[field_name],
    )


def test_cli_calls_run_atb_pipeline(monkeypatch: pytest.MonkeyPatch) -> None:
    """Ensure that cli() calls run_atb_pipeline with the settings."""
    mock_settings_instance = MagicMock()
    mock_settings_cls = MagicMock(return_value=mock_settings_instance)
    mock_run_atb_pipeline = MagicMock()

    monkeypatch.setattr(all_the_bacteria, "AtbSettings", mock_settings_cls)
    monkeypatch.setattr(all_the_bacteria, "run_atb_pipeline", mock_run_atb_pipeline)

    cli()

    mock_settings_cls.assert_called_once_with()
    mock_run_atb_pipeline.assert_called_once_with(mock_settings_instance)


def test_load_patterns_returns_compiled_pattern(pattern_file: Path) -> None:
    """load_patterns() should return a compiled re.Pattern object."""
    assert isinstance(load_patterns(pattern_file), re.Pattern)


def test_load_patterns_exact_match(pattern_file: Path) -> None:
    """Plain text lines should match exactly against their literal content."""
    pattern = load_patterns(pattern_file)
    assert isinstance(pattern, re.Pattern)
    assert pattern.match("hello world")
    assert pattern.match("foo.bar")
    assert pattern.match("2+2=4")

    # patterns are anchored with ^ and $, so partial matches should not succeed
    assert not pattern.match("hello world!!!")
    assert not pattern.match("well, hello world")

    # dot is literal, not a wildcard
    assert not pattern.match("fooXbar")
    # + is literal, not a quantifier
    assert not pattern.match("22=4")

    # line ending with * should match the prefix followed by any suffix
    assert pattern.match("starts with")
    assert pattern.match("starts with anything")
    assert pattern.match("starts with 123!@#")

    # line ending with * should not match strings with a different prefix
    assert not pattern.match("ends with")


def test_load_patterns_blank_lines_are_ignored(tmp_path: Path) -> None:
    """Blank lines in the file should be silently skipped."""
    p = tmp_path / "patterns.txt"
    p.write_text("\nhello\n\nworld\n", encoding="utf-8")
    pattern = load_patterns(p)
    assert isinstance(pattern, re.Pattern)
    assert pattern.match("hello")
    assert pattern.match("world")
    assert not pattern.match("")


def test_load_patterns_only_wildcard_matches_anything(tmp_path: Path) -> None:
    """A file containing only '*' should produce a pattern that matches any string."""
    p = tmp_path / "patterns.txt"
    p.write_text("*\n", encoding="utf-8")
    pattern = load_patterns(p)
    assert isinstance(pattern, re.Pattern)
    assert pattern.match("")
    assert pattern.match("anything at all")


def test_load_patterns_alternation(tmp_path: Path) -> None:
    """Each line in the file should become an alternative in the combined pattern."""
    p = tmp_path / "patterns.txt"
    p.write_text("cat\ndog\nbird\n", encoding="utf-8")
    pattern = load_patterns(p)
    assert isinstance(pattern, re.Pattern)
    assert pattern.match("cat")
    assert pattern.match("dog")
    assert pattern.match("bird")
    assert not pattern.match("fish")


def test_load_patterns_no_file_returns_none(tmp_path: Path) -> None:
    """Ensure that loading a non-existent file returns None."""
    pattern = load_patterns(tmp_path / "some" / "path")
    assert pattern is None


def test_load_patterns_touched_file_returns_none(tmp_path: Path) -> None:
    """Ensure that loading an empty file returns None."""
    p = tmp_path / "patterns.txt"
    p.touch()
    pattern = load_patterns(p)
    assert pattern is None


def test_load_patterns_empty_file_returns_none(tmp_path: Path) -> None:
    """Ensure that loading an empty file returns None."""
    p = tmp_path / "patterns.txt"
    p.write_text("\n\n   \t\n  \n   \n\t\t\n", encoding="utf-8")
    pattern = load_patterns(p)
    assert pattern is None


@pytest.mark.vcr
def test_download_atb_index_tsv_vcr(test_settings: AtbSettings) -> None:
    """Ensure that the download_atb_index function fetches the correct file."""
    output_file = download_atb_index_tsv(test_settings)
    raw_data_dir = Path(test_settings.raw_data_dir)
    assert raw_data_dir.is_dir()
    assert output_file.exists()
    assert output_file.parent == raw_data_dir


@pytest.mark.default_cassette("test_download_atb_index_tsv_vcr.yaml")
def test_download_atb_index_tsv_vcr_destination_s3(test_s3_settings: AtbSettings) -> None:
    """Ensure that the download_atb_index function fetches the correct file."""
    mock_download_client = MagicMock()
    mock_stream_to_s3 = MagicMock()
    mock_requests = MagicMock()
    with (
        patch(
            "cdm_data_loaders.pipelines.all_the_bacteria.FileDownloader",
            return_value=mock_download_client,
        ),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.stream_to_s3", mock_stream_to_s3),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.requests", mock_requests),
    ):
        output_file = download_atb_index_tsv(test_s3_settings)

    download_url = "https://osf.io/download/r6gcp/"
    mock_stream_to_s3.assert_called_once_with(
        url=download_url,
        s3_path=f"{test_s3_settings.raw_data_dir}/{ALL_ATB_FILE_NAME}",
        requests=mock_requests,
    )
    mock_download_client.download.assert_called_once_with(url=download_url, destination=Path(ALL_ATB_FILE_NAME))
    assert output_file == Path(ALL_ATB_FILE_NAME)


@pytest.mark.vcr
def test_download_atb_index_tsv_error_404(test_settings: AtbSettings) -> None:
    """Ensure that a 404 response causes an error and the function to die."""
    with pytest.raises(HTTPError, match="404 Client Error"):
        download_atb_index_tsv(test_settings)


@pytest.mark.default_cassette("test_download_atb_index_tsv_vcr.yaml")
def test_download_atb_index_s3_error_boom(test_s3_settings: AtbSettings, caplog: pytest.LogCaptureFixture) -> None:
    """Ensure that an error in the s3 upload causes things to die unpleasantly."""
    mock_download_client = MagicMock()
    mock_stream_to_s3 = MagicMock(side_effect=ValueError("ZOMG!"))
    mock_requests = MagicMock()
    with (
        patch(
            "cdm_data_loaders.pipelines.all_the_bacteria.FileDownloader",
            return_value=mock_download_client,
        ),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.stream_to_s3", mock_stream_to_s3),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.requests", mock_requests),
        pytest.raises(ValueError, match="ZOMG!"),
    ):
        download_atb_index_tsv(test_s3_settings)

    download_url = "https://osf.io/download/r6gcp/"
    mock_stream_to_s3.assert_called_once_with(
        url=download_url,
        s3_path=f"{test_s3_settings.raw_data_dir}/{ALL_ATB_FILE_NAME}",
        requests=mock_requests,
    )
    mock_download_client.assert_not_called()
    mock_download_client.download.assert_not_called()
    last_log_record = caplog.records.pop()
    assert last_log_record.levelno == logging.ERROR
    assert last_log_record.message == f"Could not transfer {ALL_ATB_FILE_NAME} to s3"


@pytest.mark.vcr
def test_download_atv_index_tsv_error_missing_key(test_settings: AtbSettings) -> None:
    """Ensure that the lack of a download link in the response throws an error."""
    with pytest.raises(RuntimeError, match="Could not find download URL in response from "):
        download_atb_index_tsv(test_settings)


def test_download_atv_index_cannot_create_dir() -> None:
    """Ensure that the output raw_data_dir directory can be saved to."""
    with pytest.raises(OSError, match=r"(Read-only file system|Permission denied)"):
        download_atb_index_tsv(AtbSettings(output_dir="/path/to/file"))


@pytest.mark.vcr
def test_download_atv_index_tsv_error_cannot_download_tsv(test_settings: AtbSettings) -> None:
    """Ensure that the lack of a download link in the response throws an error."""
    with pytest.raises(NonRetryableDownloadError, match="Client error: 404 Not Found"):
        download_atb_index_tsv(test_settings)


def test_get_file_download_links(test_settings: AtbSettings) -> None:
    """Ensure that the appropriate files are picked out of the ATB file index TSV file using the default matcher."""
    file_path = Path("tests") / "data" / "atb" / "all_atb_files.tsv"
    filtered_files = list(get_file_download_links(test_settings, file_path))
    # load the expected results
    expected = Path("tests") / "data" / "atb" / "filtered_files.tsv"
    with expected.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        expected_files = list(reader)
    assert len(filtered_files[0]) > 1
    assert filtered_files[0] == expected_files


EXPECTED_LINES = {
    "*": "all_atb_files.tsv",
    "AllTheBacteria/Annotation/Bakta\nAllTheBacteria/Assembly\nAllTheBacteria/Metadata\n": "assembly_bakta_metadata_exact.tsv",
    "AllTheBacteria/Annotation/Bakta*\nAllTheBacteria/Assembly\nAllTheBacteria/Metadata\n": "assembly_bakta_star_metadata.tsv",
    "AllTheBacteria/Annotation/Bakta": "bakta_exact.tsv",
    "AllTheBacteria/Annotation/Bakta*": "bakta_star.tsv",
}


@pytest.mark.parametrize("pattern_lines", EXPECTED_LINES)
def test_get_file_download_links_use_pattern_file(tmp_path: Path, pattern_lines: str) -> None:
    """Generate a pattern file from EXPECTED_LINES and check that the output from get_file_download_links is correct."""
    # create the pattern file
    p = tmp_path / "patterns.txt"
    p.write_text(f"{pattern_lines}\n", encoding="utf-8")

    settings = AtbSettings(input_dir=str(tmp_path), pattern_file="patterns.txt")
    file_path = Path("tests") / "data" / "atb" / "all_atb_files.tsv"
    filtered_files = list(get_file_download_links(settings, file_path))
    # load the expected results
    expected = Path("tests") / "data" / "atb" / EXPECTED_LINES[pattern_lines]
    with expected.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        expected_files = list(reader)
    assert len(filtered_files[0]) > 1
    assert filtered_files[0] == expected_files


def test_get_file_download_links_invalid_file(test_settings: AtbSettings, caplog: pytest.LogCaptureFixture) -> None:
    """Ensure that the correct fields are present in the ATB TSV file and throw an error if not."""
    file_path = Path("tests") / "data" / "atb" / "invalid_atb_files.tsv"
    with pytest.raises(RuntimeError, match="Missing required ATB file index TSV headers"):
        list(get_file_download_links(test_settings, file_path))
    records = caplog.records
    assert records[-1].levelno == logging.ERROR
    assert records[-1].message.startswith(
        "Missing required ATB file index TSV headers: ['filename', 'md5', 'project', 'url']"
    )
    assert records[-2].levelno == logging.WARNING
    assert records[-2].message.startswith("ATB file index TSV headers have changed.")


def test_get_file_download_links_empty_file(
    test_settings: AtbSettings, caplog: pytest.LogCaptureFixture, tmp_path: Path
) -> None:
    """Ensure that an empty file causes a runtime error."""
    file_path = tmp_path / "fake_file.tsv"
    file_path.touch()
    with pytest.raises(RuntimeError, match=f"No valid TSV data found in {file_path!s}"):
        list(get_file_download_links(test_settings, file_path))
    records = caplog.records
    assert records[-1].levelno == logging.ERROR
    assert records[-1].message == f"No valid TSV data found in {file_path!s}"


def test_get_file_download_links_no_file(test_settings: AtbSettings) -> None:
    """Ensure an error is thrown if the file cannot be found."""
    file_path = Path("/path") / "to" / "file"
    with pytest.raises(FileNotFoundError, match="No such file or directory"):
        list(get_file_download_links(test_settings, file_path))


FILE_DOWNLOADER_OUTPUT = [
    frozendict(
        {
            "filename": "file1.txt",
            "url": "https://osf.io/file1",
            "md5": "md5sum1",
            "project": "AllTheBacteria/Annotation/Project",
            "path": "Annotation/Project/file1.txt",
        }
    ),
    frozendict(
        {
            "filename": "file2.txt",
            "url": "https://osf.io/file2",
            "md5": "md5sum2",
            "project": "AllTheBacteria/Side/Project/",
            "path": "Side/Project/file2.txt",
        }
    ),
    frozendict(
        {
            "filename": "some/path/to/file3.txt",
            "url": "https://osf.io/file3",
            "md5": "md5sum3",
            "project": "AllTheBacteria",
            "path": "some/path/to/file3.txt",
        }
    ),
    frozendict(
        {
            "filename": "not/least/file4.txt",
            "url": "https://osf.io/file4",
            "md5": "md5sum4",
            "project": "AllTheBacteria/last/but",
            "path": "last/but/not/least/file4.txt",
        }
    ),
]


# osf_file_downloader tests
@pytest.mark.parametrize(
    "atb_input",
    [
        # each of the files singly and then the whole lot as a batch
        *[[f] for f in FILE_DOWNLOADER_OUTPUT],
        FILE_DOWNLOADER_OUTPUT,
    ],
)
@pytest.mark.parametrize("use_destination", VALID_DESTINATIONS)
def test_osf_file_downloader_success(
    request: pytest.FixtureRequest,
    atb_input: list[dict[str, Any]],
    use_destination: str,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Ensure that the osf_file_downloader function correctly calls the download client for each file."""
    settings = request.getfixturevalue("test_s3_settings" if use_destination == "s3" else "test_settings")

    atb_file_list = [{k: v for k, v in f.items() if k != "path"} for f in atb_input]

    # expected_output path needs the raw data dir adding to it
    expected_output = [dict(f.items()) for f in atb_input]
    for f in expected_output:
        f["path"] = f"{settings.raw_data_dir}/{f['path']}"

    mock_download_client = MagicMock()
    mock_stream_to_s3 = MagicMock()
    mock_requests = MagicMock()
    with (
        patch(
            "cdm_data_loaders.pipelines.all_the_bacteria.FileDownloader",
            return_value=mock_download_client,
        ),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.stream_to_s3", mock_stream_to_s3),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.requests", mock_requests),
    ):
        # get the output from the generator
        output = list(osf_file_downloader(settings, atb_file_list))

    # should have a separate call for each file downloaded
    if use_destination == "s3":
        mock_download_client.download.assert_not_called()
        call_args = [c.kwargs for c in mock_stream_to_s3.call_args_list]
        assert call_args == [
            {
                "url": f["url"],
                "s3_path": f["path"],
                "requests": mock_requests,
            }
            for f in atb_file_list
        ]
    else:
        assert mock_download_client.download.call_count == len(atb_file_list)

        call_list = [c.kwargs for c in mock_download_client.download.call_args_list]
        expected_calls = [
            {
                "url": f["url"],
                "destination": Path(settings.raw_data_dir) / f["path"],
                "expected_checksum": f["md5"],
                "checksum_fn": "md5",
            }
            for f in atb_input
        ]
        assert call_list == expected_calls

    # output from dlt.mark.with_table_name
    assert output[0].data == expected_output
    # the input args are mutated in place
    assert atb_file_list == expected_output

    # no logs should be emitted for successful downloads
    assert caplog.records == []


@pytest.mark.parametrize(
    ("atb_file_list", "expected_exceptions", "expected_paths"),
    [
        (
            [
                {
                    "project": "AllTheBacteria/One",
                    "filename": "Two/good_file.txt",
                    "url": "https://osf.io/good",
                    "md5": "md5sum1",
                },
                {
                    "project": "AllTheBacteria/One/Two",
                    "filename": "great_file.txt",
                    "url": "https://osf.io/great",
                    "md5": "md5sum2",
                },
                {
                    "project": "AllTheBacteria/One",
                    "filename": "fail.txt",
                    "url": "https://osf.io/fail",
                    "md5": "badmd5",
                },
            ],
            ["Could not download file from https://osf.io/fail: Loser!"],
            {"Two/good_file.txt": True, "great_file.txt": True, "fail.txt": False},
        ),
        (
            [
                {
                    "project": "Dud",
                    "filename": "bad_file.txt",
                    "url": "https://osf.io/bad",
                    "md5": "badmd5",
                },
                {
                    "project": "Dud",
                    "filename": "also_very_bad.txt",
                    "url": "https://osf.io/also_very_bad",
                    "md5": "badmd5",
                },
            ],
            [
                "Could not download file from https://osf.io/bad: BOOM!",
                "Could not download file from https://osf.io/also_very_bad: BOOM!",
            ],
            {"bad_file.txt": False, "also_very_bad.txt": False},
        ),
    ],
)
@pytest.mark.parametrize("use_destination", VALID_DESTINATIONS)
def test_osf_file_downloader_error_handling(
    atb_file_list: list[dict[str, Any]],
    expected_exceptions: list[str],
    expected_paths: dict[str, bool],
    use_destination: str,
    caplog: pytest.LogCaptureFixture,
    request: pytest.FixtureRequest,
) -> None:
    """Ensure that errors during file download are handled correctly."""
    settings = request.getfixturevalue("test_s3_settings" if use_destination == "s3" else "test_settings")

    def file_downloader_boom(**args) -> None:  # noqa: ANN003
        if "url" in args:
            if "bad" in args["url"]:
                msg = "BOOM!"
                raise ValueError(msg)
            if "fail" in args["url"]:
                msg_0 = "Loser!"
                raise RuntimeError(msg_0)

    mock_requests = MagicMock()
    mock_download_client = MagicMock()
    mock_download_client.download.side_effect = file_downloader_boom
    mock_stream_to_s3 = MagicMock(side_effect=file_downloader_boom)

    with (
        patch(
            "cdm_data_loaders.pipelines.all_the_bacteria.FileDownloader",
            return_value=mock_download_client,
        ),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.stream_to_s3", mock_stream_to_s3),
        patch("cdm_data_loaders.pipelines.all_the_bacteria.requests", mock_requests),
    ):
        list(osf_file_downloader(settings, atb_file_list))

    if use_destination == "s3":
        mock_download_client.download.assert_not_called()
        call_args = [c.kwargs for c in mock_stream_to_s3.call_args_list]
        assert call_args == [
            {
                "url": f["url"],
                "s3_path": f"{settings.raw_data_dir}/{f['project'].replace('AllTheBacteria/', '')}/{f['filename']}",
                "requests": mock_requests,
            }
            for f in atb_file_list
        ]
    else:
        mock_stream_to_s3.assert_not_called()

    expected_file_names = {
        "Two/good_file.txt": f"{settings.raw_data_dir}/One/Two/good_file.txt",
        "great_file.txt": f"{settings.raw_data_dir}/One/Two/great_file.txt",
    }

    for item in atb_file_list:
        if item["filename"] in expected_paths and expected_paths[item["filename"]]:
            assert item["path"] == expected_file_names[item["filename"]]
        else:
            assert "path" not in item

    log_messages = [r.message for r in caplog.records]
    assert expected_exceptions == log_messages


def test_run_atb_pipeline(
    test_settings: AtbSettings,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """run_atb_pipeline binds settings and delegates to run_pipeline with correct args."""
    mock_atb_file_list = MagicMock()
    mock_file_downloader = MagicMock()
    mock_run_pipeline = MagicMock()

    monkeypatch.setattr(all_the_bacteria, "atb_file_list", mock_atb_file_list)
    monkeypatch.setattr(all_the_bacteria, "file_downloader", mock_file_downloader)
    monkeypatch.setattr(all_the_bacteria, "run_pipeline", mock_run_pipeline)

    output = run_atb_pipeline(test_settings)
    assert output is None

    mock_atb_file_list.bind.assert_called_once_with(test_settings)
    mock_file_downloader.bind.assert_called_once_with(test_settings)
    mock_run_pipeline.assert_called_once_with(
        settings=test_settings,
        resource=mock_file_downloader,
        destination_kwargs={"max_table_nesting": 0},
        pipeline_kwargs={
            "pipeline_name": DATASET_NAME,
            "dataset_name": DATASET_NAME,
        },
        pipeline_run_kwargs=None,
    )


def test_run_atb_pipeline_bind_order(
    test_settings: AtbSettings,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both resources are bound before run_pipeline is called."""
    call_order: list[str] = []

    mock_atb_file_list = MagicMock()
    mock_atb_file_list.bind.side_effect = lambda *_: call_order.append("atb_file_list.bind")

    mock_file_downloader = MagicMock()
    mock_file_downloader.bind.side_effect = lambda *_: call_order.append("file_downloader.bind")

    def track_run_pipeline(**_: Any) -> None:
        call_order.append("run_pipeline")

    monkeypatch.setattr(all_the_bacteria, "atb_file_list", mock_atb_file_list)
    monkeypatch.setattr(all_the_bacteria, "file_downloader", mock_file_downloader)
    monkeypatch.setattr(all_the_bacteria, "run_pipeline", track_run_pipeline)

    run_atb_pipeline(test_settings)

    assert call_order.index("atb_file_list.bind") < call_order.index("run_pipeline")
    assert call_order.index("file_downloader.bind") < call_order.index("run_pipeline")


def test_run_atb_pipeline_resource_passed_to_run_pipeline_is_file_downloader(
    test_settings: AtbSettings,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The resource argument to run_pipeline is the (possibly mutated) file_downloader object."""
    sentinel = MagicMock(name="file_downloader_sentinel")
    mock_run_pipeline = MagicMock()

    monkeypatch.setattr(all_the_bacteria, "atb_file_list", MagicMock())
    monkeypatch.setattr(all_the_bacteria, "file_downloader", sentinel)
    monkeypatch.setattr(all_the_bacteria, "run_pipeline", mock_run_pipeline)

    run_atb_pipeline(test_settings)

    assert mock_run_pipeline.call_args.kwargs["resource"] is sentinel


@pytest.mark.parametrize("use_output_dir_for_pipeline_metadata", [True, False])
@pytest.mark.parametrize("dev_mode", [True, False])
def test_run_atb_pipeline_pipeline_dir_present_or_absent(
    use_output_dir_for_pipeline_metadata: bool,
    dev_mode: bool,
    mock_dlt: MagicMock,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Check that the appropriate args are passed as pipeline_kwargs."""
    settings = AtbSettings(
        output_dir="/my/output",
        use_destination=VALID_DESTINATIONS[0],
        use_output_dir_for_pipeline_metadata=use_output_dir_for_pipeline_metadata,
        dev_mode=dev_mode,
    )
    mock_file_downloader = MagicMock(name="file_downloader")
    if use_output_dir_for_pipeline_metadata:
        assert settings.pipeline_dir == "/my/output/.dlt_conf"
    else:
        assert settings.pipeline_dir is None

    monkeypatch.setattr(all_the_bacteria, "atb_file_list", MagicMock())
    monkeypatch.setattr(all_the_bacteria, "file_downloader", mock_file_downloader)

    run_atb_pipeline(settings)

    assert mock_dlt.destination.call_args_list == [call(VALID_DESTINATIONS[0], max_table_nesting=0)]
    mock_pipeline = mock_dlt.pipeline.return_value
    assert mock_pipeline.run.call_args_list == [call(mock_file_downloader)]

    pipeline_kwargs: dict[str, Any] = {"pipeline_name": DATASET_NAME, "dataset_name": DATASET_NAME}
    if use_output_dir_for_pipeline_metadata:
        pipeline_kwargs["pipelines_dir"] = settings.pipeline_dir
    if dev_mode:
        pipeline_kwargs["dev_mode"] = True

    assert mock_dlt.pipeline.call_args_list == [call(destination=mock_dlt.destination.return_value, **pipeline_kwargs)]
