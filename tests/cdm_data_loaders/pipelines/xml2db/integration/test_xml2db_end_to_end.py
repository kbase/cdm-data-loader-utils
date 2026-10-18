"""End-to-end tests for the xml2db pipeline, using a small library.xsd schema."""

import sys
from collections.abc import Callable
from pathlib import Path

import pytest
from dlt.common.pipeline import LoadInfo
from pandas import DataFrame

from cdm_data_loaders.pipelines.xml2db.pipeline import cli, run_xml2db_pipeline
from cdm_data_loaders.pipelines.xml2db.settings import PIPELINE_NAME, Xml2DbSettings
from tests.cdm_data_loaders.pipelines.conftest import DEFAULT_DLT_TABLES
from tests.xml_samples import SIMPLE_LIBRARY_XML

EXPECTED_BOOK_COUNT = 3
EXPECTED_LIBRARY_ROOT_COUNT_TWO_FILES = 2


def check_book_list_results(load_info: LoadInfo | None) -> tuple[DataFrame, DataFrame]:
    """Check over the output of parsing SIMPLE_LIBRARY_XML: one root row, three book rows."""
    assert load_info is not None
    assert load_info.has_failed_jobs is False
    dataset = load_info.pipeline.dataset()
    assert {"library", "book", "library_book", *DEFAULT_DLT_TABLES} == set(dataset.tables)
    book_df = dataset.table("book").df()
    library_df = dataset.table("library").df()
    assert set(book_df.columns.tolist()) >= {"pk_book", "id", "title"}
    return library_df, book_df


def test_run_xml2db_pipeline_pass_writes_expected_output(
    settings_factory: Callable[..., Xml2DbSettings],
) -> None:
    """Running the pipeline loads one row per library root and one row per book."""
    settings = settings_factory()
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    load_info = run_xml2db_pipeline(settings)
    library_df, book_df = check_book_list_results(load_info)
    assert len(library_df) == 1
    assert sorted(book_df["id"].tolist()) == ["1", "2", "3"]
    assert sorted(book_df["title"].tolist()) == ["The Shining", "The Stand", "The Tommyknockers"]


def test_run_xml2db_pipeline_pass_gzip_files_are_loaded(
    settings_factory: Callable[..., Xml2DbSettings],
    write_gzip_file: Callable[[Path, str, str], Path],
) -> None:
    """Gzip-compressed xml files matching the glob are decompressed and loaded.

    Both files contain the same three books. `book` is a "reused" table, so post-load
    compaction (on by default) collapses the two files' identical book rows down to one copy of
    each -- unlike `library`, the schema's root table, which is exempt from content-addressing
    (see `is_xml2db_content_addressed_table`) and keeps one row per source file regardless.
    """
    settings = settings_factory()
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "plain.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")
    write_gzip_file(input_dir, "library.xml.gz", SIMPLE_LIBRARY_XML)

    load_info = run_xml2db_pipeline(settings)
    library_df, book_df = check_book_list_results(load_info)
    # one root row from each of the two files: the root table is never compacted
    assert len(library_df) == EXPECTED_LIBRARY_ROOT_COUNT_TWO_FILES
    # the two files' identical books are compacted down to one copy of each
    assert sorted(book_df["id"].tolist()) == ["1", "2", "3"]


def test_run_xml2db_pipeline_pass_no_matching_files_yields_no_data_table(
    settings_factory: Callable[..., Xml2DbSettings],
) -> None:
    """An input dir with no matching xml files does not fail the run and creates no data table."""
    settings = settings_factory()
    Path(settings.input_dir).mkdir(exist_ok=True)

    load_info = run_xml2db_pipeline(settings)
    assert load_info is not None
    assert load_info.has_failed_jobs is False
    dataset = load_info.pipeline.dataset()
    # only dlt's own metadata tables exist; the data tables are never created
    assert set(dataset.tables) == DEFAULT_DLT_TABLES


def test_run_xml2db_pipeline_pass_custom_short_name_used_for_root_table(
    settings_factory: Callable[..., Xml2DbSettings],
) -> None:
    """The `short_name` setting only affects the virtual root table name for multi-root schemas.

    library.xsd has a single root element, so xml2db names the root table after that element
    ("library") regardless of `short_name`; this just confirms the setting is accepted and does
    not break a single-root schema.
    """
    settings = settings_factory(short_name="my_schema")
    input_dir = Path(settings.input_dir)
    input_dir.mkdir(exist_ok=True)
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")

    load_info = run_xml2db_pipeline(settings)
    library_df, book_df = check_book_list_results(load_info)
    assert len(library_df) == 1
    assert sorted(book_df["id"].tolist()) == ["1", "2", "3"]


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    dlt_destination_config: str,
    library_xsd: Path,
) -> None:
    """cli() reads command-line arguments and runs the pipeline, producing the expected output."""
    input_dir = tmp_path / "cli_input"
    input_dir.mkdir()
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")
    output_dir = tmp_path / "bucket"
    output_dir.mkdir(exist_ok=True)
    log_config_file = tmp_path / "logging.json"
    log_config_file.write_text('{"version": 1}')

    argv = [
        PIPELINE_NAME,
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--log-config-file",
        str(log_config_file),
        "--use-destination",
        dlt_destination_config,
        "--dataset-name",
        "cli_dataset",
        "--xsd-file",
        str(library_xsd),
        "--use-output-dir-for-pipeline-metadata",
        "false",
        "--dlt-dev-mode",
        "false",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    load_info = cli.run()
    library_df, book_df = check_book_list_results(load_info)
    assert len(library_df) == 1
    assert sorted(book_df["id"].tolist()) == ["1", "2", "3"]
