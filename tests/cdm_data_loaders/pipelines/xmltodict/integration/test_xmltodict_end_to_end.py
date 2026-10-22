"""End-to-end tests for the xmltodict pipeline."""

from collections.abc import Callable
from pathlib import Path

import pytest
from dlt.common.pipeline import LoadInfo
from pandas import DataFrame

from cdm_data_loaders.pipelines.xmltodict.settings import PIPELINE_NAME, XmlToDictSettings
from tests.cdm_data_loaders.pipelines.conftest import DEFAULT_DLT_TABLES
from tests.cdm_data_loaders.pipelines.framework import PipelineTestRunner
from tests.cdm_data_loaders.pipelines.helpers import (
    LOADER_FILE_FORMATS,
    assert_pipeline_report_written,
    assert_run_ok,
)
from tests.helpers import settings_to_argv
from tests.xml_samples import SIMPLE_LIBRARY_XML


def check_book_list_results(load_info: LoadInfo | None, table_name: str | None = None) -> DataFrame:
    """Check over the output of parsing SIMPLE_LIBRARY_XML."""
    load_info = assert_run_ok(load_info)
    dataset = load_info.pipeline.dataset()
    table_name = table_name or "book"
    assert set(dataset.tables) == {table_name, *DEFAULT_DLT_TABLES}
    book_df = dataset.table(table_name).df()
    assert set(book_df.columns.tolist()) >= {"book___id", "book__title", "_dlt_id", "_dlt_load_id"}
    return book_df


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_run_xmltodict_pipeline_pass_writes_expected_output(
    settings_factory: Callable[..., XmlToDictSettings],
    loader_file_format: str,
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Running the pipeline on xml files loads one row per matching element into the configured table."""
    settings = settings_factory(loader_file_format=loader_file_format)
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    load_info = pipeline_runner.run(PIPELINE_NAME, settings=settings)
    book_df = check_book_list_results(load_info)
    assert sorted(book_df["book___id"].tolist()) == ["1", "2", "3"]
    assert_pipeline_report_written(Path(settings.output_dir), load_info)


def test_run_xmltodict_pipeline_pass_gzip_files_are_loaded(
    settings_factory: Callable[..., XmlToDictSettings],
    write_gzip_file: Callable[[Path, str, str], Path],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Gzip-compressed xml files matching the glob are decompressed and loaded."""
    settings = settings_factory()
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "plain.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")
    write_gzip_file(input_dir, "library.xml.gz", SIMPLE_LIBRARY_XML)

    load_info = pipeline_runner.run(PIPELINE_NAME, settings=settings)
    book_df = check_book_list_results(load_info)
    # one copy of the three books from the plain file, one from the gzipped copy
    assert sorted(book_df["book___id"].tolist()) == ["1", "1", "2", "2", "3", "3"]


def test_run_xmltodict_pipeline_pass_no_matching_files_yields_no_data_table(
    settings_factory: Callable[..., XmlToDictSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """An input dir with no matching xml files does not fail the run and creates no data table."""
    settings = settings_factory()
    Path(settings.input_dir).mkdir(exist_ok=True)

    load_info = assert_run_ok(pipeline_runner.run(PIPELINE_NAME, settings=settings))
    dataset = load_info.pipeline.dataset()
    # only dlt's own metadata tables exist; the data table is never created
    assert set(dataset.tables) == DEFAULT_DLT_TABLES
    assert_pipeline_report_written(Path(settings.output_dir), load_info)


def test_run_xmltodict_pipeline_pass_custom_table_name_is_respected(
    settings_factory: Callable[..., XmlToDictSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Rows land in the table named by table_name, not in a name derived from the xml tag."""
    settings = settings_factory(table_name="library_entries")
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    load_info = pipeline_runner.run(PIPELINE_NAME, settings=settings)
    book_df = check_book_list_results(load_info, "library_entries")
    # the record dict is keyed by the xml tag, so the flattened column prefix is `book`
    assert sorted(book_df["book___id"].tolist()) == ["1", "2", "3"]


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    settings_factory: Callable[..., XmlToDictSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """cli() reads command-line arguments and runs the pipeline, producing the expected output."""
    settings = settings_factory(dataset_name="cli_dataset")
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    load_info = assert_run_ok(pipeline_runner.run(PIPELINE_NAME, cli_args=settings_to_argv(settings)))
    book_df = check_book_list_results(load_info)
    assert sorted(book_df["book___id"].tolist()) == ["1", "2", "3"]
    assert_pipeline_report_written(Path(settings.output_dir), load_info)
