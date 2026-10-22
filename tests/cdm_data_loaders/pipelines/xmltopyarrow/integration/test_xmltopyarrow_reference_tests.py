"""Integration tests for the xmltopyarrow pipeline against the uniref reference dataset."""

from collections.abc import Callable
from typing import Any

import dlt
import duckdb
import pytest

from cdm_data_loaders.core.fields import PARQUET
from tests.cdm_data_loaders.pipelines.helpers import (
    BUFFER_SIZES,
    LOADER_FILE_FORMATS,
    WORKER_CONFIGS,
    assert_dataset_matches_reference,
    assert_datasets_equal,
    read_pipeline_tables,
    reference_entries_from_duckdb,
    run_and_compare_variations,
    run_and_read,
    sorted_json,
)
from tests.cdm_data_loaders.pipelines.xmltopyarrow.integration.conftest import (
    normalize_uniref_entry_for_xsd,
    reconstruct_xsd_entries,
    strip_xsd_arrow_row,
)

CHUNK_DIRS = ("chunk_5_el", "chunk_20_el", "chunk_100_el")
EXPECTED_ENTRY_COUNT = 100


def test_xmltopyarrow_pass_diff_format_same_results(
    run_xmltopyarrow_pipeline_fxt: Callable[..., Any],
    sorted_reference_xml_entries: list[str],
) -> None:
    """Ensure that the save format does not affect the reconstructed output."""
    reconstructed_by_format: dict[str, list[str]] = {}
    for output_format in LOADER_FILE_FORMATS:
        fmt = str(output_format)
        load_info, output_dir = run_and_read(
            run_xmltopyarrow_pipeline_fxt,
            "chunk_20_el",
            loader_file_format=output_format,
        )
        tables = read_pipeline_tables(output_dir / load_info.dataset_name, fmt)
        reconstructed_by_format[fmt] = sorted_json(reconstruct_xsd_entries(tables))

    formats = list(reconstructed_by_format.keys())
    for fmt in formats[1:]:
        assert reconstructed_by_format[fmt] == reconstructed_by_format[formats[0]]
    assert reconstructed_by_format[formats[0]] == sorted_reference_xml_entries


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xmltopyarrow_pass_output_identical_across_chunk_dirs(
    run_xmltopyarrow_pipeline_fxt: Callable[..., Any],
    sorted_reference_xml_entries: list[str],
    loader_file_format: str,
) -> None:
    """Pipeline output is identical across all three chunk dirs and matches the reference."""
    run_and_compare_variations(
        run_xmltopyarrow_pipeline_fxt,
        {f"{key}_{loader_file_format}": (key, {"loader_file_format": loader_file_format}) for key in CHUNK_DIRS},
        sorted_reference_json=sorted_reference_xml_entries,
        output_format=loader_file_format,
        reconstruct=reconstruct_xsd_entries,
    )


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xmltopyarrow_pass_output_identical_across_buffer_sizes(
    run_xmltopyarrow_pipeline_fxt: Callable[..., Any],
    sorted_reference_xml_entries: list[str],
    loader_file_format: str,
) -> None:
    """ArrowListBuffer size does not affect pipeline output for any chunk dir."""
    for chunk_dir in CHUNK_DIRS:
        run_and_compare_variations(
            run_xmltopyarrow_pipeline_fxt,
            {
                f"{buffer_size}_{loader_file_format}": (
                    chunk_dir,
                    {"buffer_size": buffer_size, "loader_file_format": loader_file_format},
                )
                for buffer_size in BUFFER_SIZES
            },
            sorted_reference_json=sorted_reference_xml_entries,
            output_format=loader_file_format,
            reconstruct=reconstruct_xsd_entries,
        )


@pytest.mark.parametrize(
    ("extract_workers", "normalize_workers", "load_workers"),
    WORKER_CONFIGS,
    ids=[f"extract_{e}-normalize_{n}-load_{l}" for e, n, l in WORKER_CONFIGS],  # noqa: E741
)
@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS, ids=str)
def test_xmltopyarrow_pass_no_data_lost_across_worker_configs(
    *,
    run_xmltopyarrow_pipeline_fxt: Callable[..., Any],
    sorted_reference_xml_entries: list[str],
    extract_workers: int,
    normalize_workers: int,
    load_workers: int,
    loader_file_format: str,
) -> None:
    """No data is lost regardless of how the job is split across extract/normalize/load workers."""
    load_info_and_dir = {}
    for worker_key in ("baseline", "current"):
        worker_config = (
            {}
            if worker_key == "baseline"
            else {
                "extract.workers": extract_workers,
                "normalize.workers": normalize_workers,
                "load.workers": load_workers,
            }
        )
        with dlt.config.values(worker_config):
            load_info_and_dir[worker_key] = run_and_read(
                run_xmltopyarrow_pipeline_fxt, "chunk_5_el", loader_file_format=loader_file_format
            )

    assert_datasets_equal(load_info_and_dir)
    assert_dataset_matches_reference(
        load_info_and_dir,
        sorted_reference_xml_entries,
        loader_file_format,
        reconstruct=reconstruct_xsd_entries,
    )


def test_xmltopyarrow_pass_uniref_chunk_5_el_matches_reference(
    run_xmltopyarrow_pipeline_fxt: Callable[..., Any],
    reference_xml_entries: list[dict[str, Any]],
) -> None:
    """Running on chunk_5_el with uniref.xsd produces a single Arrow-typed entry table matching the 100 reference entries."""
    (load_info, output_dir) = run_and_read(
        run_xmltopyarrow_pipeline_fxt,
        "chunk_5_el",
    )
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, PARQUET)
    assert set(tables) == {"entry"}
    loaded_rows = [strip_xsd_arrow_row(row) for row in tables["entry"]]
    expected_rows = [normalize_uniref_entry_for_xsd(entry) for entry in reference_xml_entries]
    assert len(loaded_rows) == EXPECTED_ENTRY_COUNT
    assert sorted_json(loaded_rows) == sorted_json(expected_rows)


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xmltopyarrow_pass_reconstructed_matches_duckdb_reference(
    run_xmltopyarrow_pipeline_fxt: Callable[..., Any],
    reference_xml_dataset: duckdb.DuckDBPyConnection,
    loader_file_format: str,
) -> None:
    """Reconstructed pipeline entries equal the entries loaded from the DuckDB reference."""
    (load_info, output_dir) = run_and_read(
        run_xmltopyarrow_pipeline_fxt,
        "chunk_100_el",
        loader_file_format=loader_file_format,
    )
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, loader_file_format)
    reconstructed = reconstruct_xsd_entries(tables)
    from_duckdb = [
        normalize_uniref_entry_for_xsd(entry) for entry in reference_entries_from_duckdb(reference_xml_dataset)
    ]
    assert len(from_duckdb) == EXPECTED_ENTRY_COUNT
    assert sorted_json(reconstructed) == sorted_json(from_duckdb)
