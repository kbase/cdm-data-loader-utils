"""End-to-end tests for the xmltopyarrow pipeline."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from dlt.common.pipeline import LoadInfo
from pandas import DataFrame

from cdm_data_loaders.pipelines.xmltopyarrow.settings import PIPELINE_NAME, XmlToPyArrowSettings
from tests.cdm_data_loaders.pipelines.conftest import DEFAULT_DLT_TABLES
from tests.cdm_data_loaders.pipelines.framework import PipelineTestRunner
from tests.cdm_data_loaders.pipelines.helpers import (
    LOADER_FILE_FORMATS,
    assert_pipeline_report_written,
    assert_run_ok,
    settings_to_argv,
    sorted_json,
)
from tests.cdm_data_loaders.pipelines.xmltopyarrow.integration.conftest import (
    normalize_uniref_entry_for_xsd,
    strip_xsd_arrow_row,
)
from tests.conftest import N_REFERENCE_XML_ENTRIES
from tests.xml_samples import SIMPLE_LIBRARY_XML


def check_book_list_results(load_info: LoadInfo | None, table_name: str | None = None) -> DataFrame:
    """Check over the output of parsing SIMPLE_LIBRARY_XML with LIBRARY_XSD."""
    load_info = assert_run_ok(load_info)
    dataset = load_info.pipeline.dataset()
    table_name = table_name or "book"
    assert set(dataset.tables) == {table_name, *DEFAULT_DLT_TABLES}
    book_df = dataset.table(table_name).df()
    assert set(book_df.columns.tolist()) >= {"_id", "title"}
    return book_df


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_run_xmltopyarrow_pipeline_pass_writes_expected_output(
    settings_factory: Callable[..., XmlToPyArrowSettings],
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
    assert sorted(book_df["_id"].tolist()) == ["1", "2", "3"]
    assert_pipeline_report_written(Path(settings.output_dir), load_info)


def test_run_xmltopyarrow_pipeline_pass_gzip_files_are_loaded(
    settings_factory: Callable[..., XmlToPyArrowSettings],
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
    assert sorted(book_df["_id"].tolist()) == ["1", "1", "2", "2", "3", "3"]


def test_run_xmltopyarrow_pipeline_pass_no_matching_files_yields_no_data_table(
    settings_factory: Callable[..., XmlToPyArrowSettings],
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


def test_run_xmltopyarrow_pipeline_pass_custom_table_name_is_respected(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Rows land in the table named by table_name, not in a name derived from the xml tag."""
    settings = settings_factory(table_name="library_entries")
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    load_info = pipeline_runner.run(PIPELINE_NAME, settings=settings)
    book_df = check_book_list_results(load_info, "library_entries")
    assert sorted(book_df["_id"].tolist()) == ["1", "2", "3"]


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """cli() reads command-line arguments and runs the pipeline, producing the expected output."""
    settings = settings_factory(dataset_name="cli_dataset")
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    load_info = assert_run_ok(pipeline_runner.run(PIPELINE_NAME, cli_args=settings_to_argv(settings)))
    book_df = check_book_list_results(load_info)
    assert sorted(book_df["_id"].tolist()) == ["1", "2", "3"]
    assert_pipeline_report_written(Path(settings.output_dir), load_info)


def test_run_xmltopyarrow_pipeline_pass_uniref_chunk_5_el(
    settings_factory: Callable[..., XmlToPyArrowSettings],
    reference_xml_data_dir: Path,
    reference_xml_entries: list[dict[str, Any]],
    pipeline_runner: PipelineTestRunner,
) -> None:
    """Running the pipeline on uniref chunk_5_el with uniref.xsd loads all 100 entries with Arrow types."""
    settings = settings_factory(
        input_dir=str(reference_xml_data_dir / "chunk_5_el"),
        xsd_file=str(reference_xml_data_dir / "uniref.xsd"),
        xml_tag="entry",
        table_name="entry",
    )

    load_info = assert_run_ok(pipeline_runner.run(PIPELINE_NAME, settings=settings))
    dataset = load_info.pipeline.dataset()
    assert set(dataset.tables) == {"entry", *DEFAULT_DLT_TABLES}

    loaded_rows = [strip_xsd_arrow_row(row) for row in dataset.table("entry").arrow().to_pylist()]
    expected_rows = [normalize_uniref_entry_for_xsd(entry) for entry in reference_xml_entries]
    assert len(loaded_rows) == N_REFERENCE_XML_ENTRIES
    assert sorted_json(loaded_rows) == sorted_json(expected_rows)
    assert_pipeline_report_written(Path(settings.output_dir), load_info)
