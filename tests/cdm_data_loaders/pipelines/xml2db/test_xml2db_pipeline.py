"""Unit tests for run_xml2db_pipeline, _read_items, _maybe_compact, and cli wiring."""

import gzip
import sys
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock, patch
from uuid import uuid4

import pendulum
import pytest
from dlt.common.pipeline import LoadInfo
from dlt.common.storages.fsspec_filesystem import FileItemDict
from dlt.extract import DltResource
from xml2db import DataModel

import cdm_data_loaders.pipelines.xml2db.pipeline as xml2db_pipeline_module
from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.pipelines.xml2db.pipeline import _read_items, cli, run_xml2db_pipeline
from cdm_data_loaders.pipelines.xml2db.settings import PIPELINE_NAME, Xml2DbSettings
from cdm_data_loaders.readers.xml2db_doc import build_xml2db_model
from tests.cdm_data_loaders.pipelines.conftest import duckdb_pipeline

SIMPLE_LIBRARY_XML = """<?xml version="1.0"?>
<library>
    <book id="1"><title>The Shining</title></book>
    <book id="2"><title>The Stand</title></book>
</library>
"""

LIBRARY_XSD = """<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:element name="library">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="book" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="title" type="xs:string"/>
                        </xs:sequence>
                        <xs:attribute name="id" type="xs:string" use="required"/>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
"""


@pytest.fixture
def library_xsd(tmp_path: Path) -> Path:
    """Write the small library.xsd schema used for pipeline-mechanics tests."""
    xsd_path = tmp_path / "library.xsd"
    xsd_path.write_text(LIBRARY_XSD, encoding="utf-8")
    return xsd_path


@pytest.fixture
def library_model(library_xsd: Path) -> DataModel:
    """Build the xml2db DataModel from library.xsd."""
    return build_xml2db_model(library_xsd, short_name="library_test")


def make_file_item(file_path: Path) -> FileItemDict:
    """Build a dlt FileItemDict for a local file, as the filesystem source would."""
    return FileItemDict(
        {
            "file_url": file_path.as_uri(),
            "file_name": file_path.name,
            "relative_path": file_path.name,
            "mime_type": "application/xml",
            "modification_date": pendulum.now(),
            "size_in_bytes": file_path.stat().st_size,
        }
    )


def fake_settings(
    tmp_path: Path,
    xsd_file: Path,
    **overrides: Any,  # noqa: ANN401
) -> Xml2DbSettings:
    """Build a valid Xml2DbSettings with fields open to override."""
    log_config_file = tmp_path / "logging.json"
    log_config_file.write_text('{"version": 1}', encoding="utf-8")
    input_dir = tmp_path / "input"
    input_dir.mkdir(exist_ok=True)
    output_dir = tmp_path / "output"
    output_dir.mkdir(exist_ok=True)

    kwargs: dict[str, Any] = {
        "buffer_size": 100,
        "log_interval": 1000,
        "dataset_name": "xml2db_test_dataset",
        "loader_file_format": "parquet",
        "log_config_file": str(log_config_file),
        "input_dir": str(input_dir),
        "output_dir": str(output_dir),
        "dlt_dev_mode": False,
        "use_destination": "local_fs",
        "use_output_dir_for_pipeline_metadata": False,
        "xsd_file": str(xsd_file),
    }
    kwargs.update(overrides)
    return Xml2DbSettings(**kwargs)


# _read_items


def test_read_items_pass_yields_pages_of_rows_per_table(tmp_path: Path, library_model: DataModel) -> None:
    """_read_items parses each file item and yields table-tagged pages of rows."""
    xml_file = tmp_path / "data.xml"
    xml_file.write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")
    settings = fake_settings(tmp_path, tmp_path / "library.xsd")

    items = list(_read_items(iter([make_file_item(xml_file)]), settings, library_model))

    assert items, "no pages yielded"
    for item in items:
        assert isinstance(item.data, list)
    table_names = [item.meta.table_name for item in items]
    assert "book" in table_names


def test_read_items_pass_reads_every_file_in_items(tmp_path: Path, library_model: DataModel) -> None:
    """One page set is produced per file, and every file's rows are present."""
    settings = fake_settings(tmp_path, tmp_path / "library.xsd")
    file_paths = []
    for index in range(3):
        xml_file = tmp_path / f"data_{index}.xml"
        xml_file.write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")
        file_paths.append(xml_file)

    items = list(_read_items(iter(make_file_item(path) for path in file_paths), settings, library_model))

    book_rows = [row for item in items if item.meta.table_name == "book" for row in item.data]
    assert len(book_rows) == 6  # noqa: PLR2004 -- 2 books per file x 3 files
    ids = {row["id"] for row in book_rows}
    assert ids == {"1", "2"}


def test_read_items_pass_gzip_files_are_decompressed(tmp_path: Path, library_model: DataModel) -> None:
    """A gzipped input file parses the same as its decompressed copy."""
    xml_file = tmp_path / "data.xml.gz"
    with gzip.open(xml_file, "wb") as fh:
        fh.write(SIMPLE_LIBRARY_XML.encode("utf-8"))
    settings = fake_settings(tmp_path, tmp_path / "library.xsd")

    items = list(_read_items(iter([make_file_item(xml_file)]), settings, library_model))

    book_rows = [row for item in items if item.meta.table_name == "book" for row in item.data]
    assert len(book_rows) == 2  # noqa: PLR2004


# run_xml2db_pipeline


@pytest.mark.parametrize("loader_file_format", ["jsonl", "parquet"])
def test_run_xml2db_pipeline_pass_sets_core_run_pipeline_args_correctly(
    tmp_path: Path, library_xsd: Path, loader_file_format: str
) -> None:
    """run_xml2db_pipeline builds the reader and delegates to run_pipeline with correct args."""
    settings = fake_settings(tmp_path, library_xsd, loader_file_format=loader_file_format)

    with patch.object(xml2db_pipeline_module, "run_pipeline") as mock_run_pipeline:
        run_xml2db_pipeline(settings)

    assert mock_run_pipeline.call_count == 1
    _, kwargs = mock_run_pipeline.call_args
    assert kwargs.keys() == {"settings", "resource", "pipeline_kwargs", "pipeline_run_kwargs"}
    assert kwargs["settings"] == settings
    assert kwargs["pipeline_kwargs"] == {
        "pipeline_name": PIPELINE_NAME,
        "dataset_name": settings.dataset_name,
    }
    assert kwargs["pipeline_run_kwargs"] == {"loader_file_format": loader_file_format}
    assert isinstance(kwargs["resource"], DltResource)


def test_run_xml2db_pipeline_pass_binds_reader_before_run_pipeline(tmp_path: Path, library_xsd: Path) -> None:
    """The resource passed to run_pipeline is the bound reader piped into the filesystem source."""
    settings = fake_settings(tmp_path, library_xsd)

    with patch.object(xml2db_pipeline_module, "run_pipeline") as mock_run_pipeline:
        run_xml2db_pipeline(settings)

    mock_run_pipeline.assert_called_once()
    resource = mock_run_pipeline.call_args.kwargs["resource"]
    assert resource.name == "xml2db_reader"
    assert resource._pipe.parent.name == "filesystem"  # noqa: SLF001


def test_run_xml2db_pipeline_pass_builds_model_with_settings_schema_and_short_name(
    tmp_path: Path, library_xsd: Path
) -> None:
    """The DataModel is built from settings.xsd_file and settings.short_name."""
    settings = fake_settings(tmp_path, library_xsd, short_name="my_short_name")

    with (
        patch.object(
            xml2db_pipeline_module, "build_xml2db_model", wraps=xml2db_pipeline_module.build_xml2db_model
        ) as mock_build,
        patch.object(xml2db_pipeline_module, "run_pipeline"),
    ):
        run_xml2db_pipeline(settings)

    mock_build.assert_called_once()
    args, kwargs = mock_build.call_args
    assert args[0] == settings.xsd_file
    assert kwargs["short_name"] == "my_short_name"
    assert kwargs["model_config"] is None


def test_run_xml2db_pipeline_pass_loads_model_config_from_xml2db_config_file(tmp_path: Path, library_xsd: Path) -> None:
    """A settings.xml2db_config_file is loaded via xml2db.load_config and passed to the model."""
    config_file = tmp_path / "model_config.yaml"
    config_file.write_text("record_hash_size: 16\n", encoding="utf-8")
    settings = fake_settings(tmp_path, library_xsd, xml2db_config_file=str(config_file))

    with (
        patch.object(
            xml2db_pipeline_module, "build_xml2db_model", wraps=xml2db_pipeline_module.build_xml2db_model
        ) as mock_build,
        patch.object(xml2db_pipeline_module, "run_pipeline"),
    ):
        run_xml2db_pipeline(settings)

    mock_build.assert_called_once()
    _, kwargs = mock_build.call_args
    assert kwargs["model_config"] == {"record_hash_size": 16}


# _maybe_compact


def make_load_info(dataset_name: str, has_failed_jobs: bool = False) -> LoadInfo:
    """Build a real LoadInfo-shaped object for _maybe_compact tests."""
    load_info = MagicMock(spec=LoadInfo)
    load_info.dataset_name = dataset_name
    load_info.has_failed_jobs = has_failed_jobs
    return load_info


def test_maybe_compact_pass_runs_when_chunked_and_enabled(tmp_path: Path, library_model: DataModel) -> None:
    """Compaction runs for a chunked, local_fs run with compact_reused_tables set."""
    settings = fake_settings(tmp_path, tmp_path / "library.xsd", xml_tag="book")
    load_info = make_load_info(settings.dataset_name)

    with patch.object(xml2db_pipeline_module, "compact_reused_tables") as mock_compact:
        xml2db_pipeline_module._maybe_compact(settings, library_model, load_info)  # noqa: SLF001

    mock_compact.assert_called_once()
    args, _kwargs = mock_compact.call_args
    assert args[0] == Path(settings.output_dir) / settings.dataset_name
    assert args[1] is library_model
    assert args[2] == "parquet"


def test_maybe_compact_pass_runs_when_not_chunked_and_enabled(tmp_path: Path, library_model: DataModel) -> None:
    """Compaction also runs for an unchunked, local_fs run with compact_reused_tables set.

    A multi-file unchunked run can leave the same kind of literal duplicate rows behind as a
    chunked one (one xml2db Document per file, same as one per chunk), so compaction is not
    gated on xml_tag being set.
    """
    settings = fake_settings(tmp_path, tmp_path / "library.xsd")
    load_info = make_load_info(settings.dataset_name)

    with patch.object(xml2db_pipeline_module, "compact_reused_tables") as mock_compact:
        xml2db_pipeline_module._maybe_compact(settings, library_model, load_info)  # noqa: SLF001

    mock_compact.assert_called_once()
    args, _kwargs = mock_compact.call_args
    assert args[0] == Path(settings.output_dir) / settings.dataset_name
    assert args[1] is library_model
    assert args[2] == "parquet"


def test_maybe_compact_pass_skipped_when_compact_reused_tables_false(tmp_path: Path, library_model: DataModel) -> None:
    """No compaction when compact_reused_tables is off, even for a chunked run."""
    settings = fake_settings(tmp_path, tmp_path / "library.xsd", xml_tag="book", compact_reused_tables=False)
    load_info = make_load_info(settings.dataset_name)

    with patch.object(xml2db_pipeline_module, "compact_reused_tables") as mock_compact:
        xml2db_pipeline_module._maybe_compact(settings, library_model, load_info)  # noqa: SLF001

    mock_compact.assert_not_called()


def test_maybe_compact_pass_skipped_when_load_info_is_none(tmp_path: Path, library_model: DataModel) -> None:
    """No compaction when the pipeline run failed to start (load_info is None)."""
    settings = fake_settings(tmp_path, tmp_path / "library.xsd", xml_tag="book")

    with patch.object(xml2db_pipeline_module, "compact_reused_tables") as mock_compact:
        xml2db_pipeline_module._maybe_compact(settings, library_model, None)  # noqa: SLF001

    mock_compact.assert_not_called()


def test_maybe_compact_pass_skipped_when_run_has_failed_jobs(tmp_path: Path, library_model: DataModel) -> None:
    """No compaction when the load had failed jobs."""
    settings = fake_settings(tmp_path, tmp_path / "library.xsd", xml_tag="book")
    load_info = make_load_info(settings.dataset_name, has_failed_jobs=True)

    with patch.object(xml2db_pipeline_module, "compact_reused_tables") as mock_compact:
        xml2db_pipeline_module._maybe_compact(settings, library_model, load_info)  # noqa: SLF001

    mock_compact.assert_not_called()


def test_maybe_compact_pass_logs_when_duplicates_removed(
    tmp_path: Path, library_model: DataModel, caplog: pytest.LogCaptureFixture
) -> None:
    """A non-empty compaction result is logged."""
    settings = fake_settings(tmp_path, tmp_path / "library.xsd", xml_tag="book")
    load_info = make_load_info(settings.dataset_name)

    with patch.object(xml2db_pipeline_module, "compact_reused_tables", return_value={"book": 3}):
        xml2db_pipeline_module._maybe_compact(settings, library_model, load_info)  # noqa: SLF001

    assert any("compaction" in record.getMessage().lower() for record in caplog.records)


# cli


def test_cli_pass_calls_run_cli_with_settings_class_and_pipeline_fn() -> None:
    """Cli is a PipelineCli wired to Xml2DbSettings and run_xml2db_pipeline."""
    assert isinstance(cli, PipelineCli)
    assert cli.settings_cls is Xml2DbSettings
    assert cli.pipeline_fn is run_xml2db_pipeline


def test_cli_pass_runs_end_to_end_from_command_line_arguments(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """cli() reads command-line arguments and loads rows via a real DuckDB pipeline."""
    input_dir = tmp_path / "cli_input"
    input_dir.mkdir()
    (input_dir / "library.xml").write_text(SIMPLE_LIBRARY_XML, encoding="utf-8")
    output_dir = tmp_path / "cli_output"
    output_dir.mkdir()
    xsd_file = tmp_path / "library.xsd"
    xsd_file.write_text(LIBRARY_XSD, encoding="utf-8")
    log_config_file = tmp_path / "logging.json"
    log_config_file.write_text('{"version": 1}', encoding="utf-8")

    argv = [
        "xml2db",
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--use-destination",
        "local_fs",
        "--log-config-file",
        str(log_config_file),
        "--dataset-name",
        "cli_dataset",
        "--xsd-file",
        str(xsd_file),
        "--use-output-dir-for-pipeline-metadata",
        "false",
        "--dlt-dev-mode",
        "false",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    # unique pipeline name: the duckdb database file is derived from it, so
    # repeated runs of the suite do not see rows from previous runs
    pipeline_name = f"test_xml2db_cli_pipeline_{uuid4().hex}"
    captured: dict[str, Any] = {}

    def fake_run_pipeline(
        *,
        resource: Any,  # noqa: ANN401
        pipeline_kwargs: dict[str, Any],
        **_: Any,  # noqa: ANN401
    ) -> None:
        assert pipeline_kwargs == {"pipeline_name": PIPELINE_NAME, "dataset_name": "cli_dataset"}
        pipeline = duckdb_pipeline(tmp_path, pipeline_name, "cli_dataset")
        captured["load_info"] = pipeline.run(resource)

    with patch.object(xml2db_pipeline_module, "run_pipeline", fake_run_pipeline):
        cli()

    load_info = captured["load_info"]
    assert not load_info.has_failed_jobs

    pipeline = duckdb_pipeline(tmp_path, pipeline_name, "cli_dataset")
    with (
        pipeline.sql_client() as client,
        client.execute_query("SELECT COUNT(*) FROM book") as cur,
    ):
        (row_count,) = cur.fetchone()
    assert row_count == 2  # noqa: PLR2004
