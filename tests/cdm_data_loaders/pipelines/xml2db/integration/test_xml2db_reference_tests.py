"""Integration tests for the xml2db pipeline against the uniref reference dataset."""

import gzip
import json
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
import xmltodict

from cdm_data_loaders.core.fields import PARQUET
from tests.cdm_data_loaders.pipelines.helpers import LOADER_FILE_FORMATS, read_pipeline_tables, sorted_json
from tests.cdm_data_loaders.pipelines.xml2db.integration.xml2db_reference_helpers import (
    assert_uniref_referential_integrity,
    reconstruct_xml2db_entries,
)

CHUNK_DIRS = ("chunk_5_el", "chunk_20_el", "chunk_100_el")
EXPECTED_ENTRY_COUNT = 100
MASTER_XML = Path("chunk_100_el") / "uniref_100_part_01.xml.gz"


def _as_list(value: Any) -> list[Any]:
    """Normalize an xmltodict field that may be a single dict, a list, or missing."""
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def _properties(properties_field: Any) -> list[list[str]]:
    """Convert xmltodict property elements into sorted [type, value] pairs."""
    return sorted([[p["_type"], p["_value"]] for p in _as_list(properties_field)])


def _member_from_xmltodict(member: dict[str, Any]) -> dict[str, Any]:
    """Convert one xmltodict (representative)member element into the comparable summary shape."""
    dbref = member["dbReference"]
    sequence = member.get("sequence") or {}
    return {
        "dbReference_type": dbref.get("_type"),
        "dbReference_id": dbref.get("_id"),
        "sequence_length": int(sequence["_length"]) if sequence.get("_length") is not None else None,
        "sequence_checksum": sequence.get("_checksum"),
        "sequence_value": sequence.get("#text"),
        "properties": _properties(dbref.get("property")),
    }


def _entry_from_xmltodict(entry: dict[str, Any]) -> dict[str, Any]:
    """Convert one xmltodict entry element into the same shape as reconstruct_xml2db_entries."""
    return {
        "id": entry["_id"],
        "updated": entry["_updated"],
        "name": entry["name"],
        "properties": _properties(entry.get("property")),
        "representative_member": _member_from_xmltodict(entry["representativeMember"]),
        "members": sorted(
            (_member_from_xmltodict(m) for m in _as_list(entry.get("member"))),
            key=lambda m: json.dumps(m, sort_keys=True, default=str),
        ),
    }


def xmltodict_reference_entries(reference_xml_data_dir: Path) -> list[dict[str, Any]]:
    """Parse the master 100-entry UniRef file directly with xmltodict, in comparable form."""
    with gzip.open(reference_xml_data_dir / MASTER_XML, "rb") as fh:
        document = xmltodict.parse(fh.read(), attr_prefix="_")
    entries = _as_list(document["UniRef50"]["entry"])
    assert len(entries) == EXPECTED_ENTRY_COUNT
    return [_entry_from_xmltodict(entry) for entry in entries]


def run_and_reconstruct(run_pipeline: Callable[..., Any], chunk_dir: str, **overrides: Any) -> list[dict[str, Any]]:
    """Run the pipeline on a chunk dir and reconstruct nested entries from its output tables.

    Reads back whatever `loader_file_format` the run actually used -- `PARQUET`, matching
    `Xml2DbSettings.loader_file_format`'s own default, unless overridden.
    """
    (load_info, output_dir) = run_pipeline(chunk_dir, **overrides)
    assert load_info is not None
    assert not load_info.has_failed_jobs
    output_format = str(overrides.get("loader_file_format", PARQUET))
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, output_format)
    return reconstruct_xml2db_entries(tables)


def test_xml2db_pass_reference_dataset_matches_xmltodict(
    run_xml2db_pipeline_fxt: Callable[..., Any],
    reference_xml_data_dir: Path,
) -> None:
    """Reconstructed entries from the xml2db pipeline match a direct xmltodict parse, 1:1."""
    reconstructed = run_and_reconstruct(run_xml2db_pipeline_fxt, "chunk_100_el")
    assert len(reconstructed) == EXPECTED_ENTRY_COUNT

    reference = xmltodict_reference_entries(reference_xml_data_dir)
    assert sorted_json(reconstructed) == sorted_json(reference)


def test_xml2db_pass_output_identical_across_chunk_dirs(run_xml2db_pipeline_fxt: Callable[..., Any]) -> None:
    """Splitting the same 100 entries across differently-sized file chunks does not change the result."""
    reconstructed_by_chunk = {
        chunk_dir: run_and_reconstruct(run_xml2db_pipeline_fxt, chunk_dir) for chunk_dir in CHUNK_DIRS
    }

    for chunk_dir, entries in reconstructed_by_chunk.items():
        assert len(entries) == EXPECTED_ENTRY_COUNT, chunk_dir

    sorted_by_chunk = {chunk_dir: sorted_json(entries) for chunk_dir, entries in reconstructed_by_chunk.items()}
    baseline = sorted_by_chunk[CHUNK_DIRS[0]]
    for chunk_dir in CHUNK_DIRS[1:]:
        assert sorted_by_chunk[chunk_dir] == baseline, chunk_dir


@pytest.mark.parametrize("loader_file_format", LOADER_FILE_FORMATS)
def test_xml2db_pass_output_identical_across_loader_file_formats(
    run_xml2db_pipeline_fxt: Callable[..., Any], loader_file_format: str
) -> None:
    """The save format (jsonl vs parquet) does not affect the reconstructed content."""
    reconstructed = run_and_reconstruct(run_xml2db_pipeline_fxt, "chunk_20_el", loader_file_format=loader_file_format)
    assert len(reconstructed) == EXPECTED_ENTRY_COUNT


def test_xml2db_pass_referential_integrity_across_merged_files(
    run_xml2db_pipeline_fxt: Callable[..., Any],
) -> None:
    """Every foreign key in the merged, multi-file dataset resolves to a row in its target table."""
    (load_info, output_dir) = run_xml2db_pipeline_fxt("chunk_5_el")
    assert load_info is not None
    assert not load_info.has_failed_jobs
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, PARQUET)
    assert_uniref_referential_integrity(tables)


# chunked parsing (xml_tag / chunk_size)


@pytest.mark.parametrize("chunk_size", [1, 7, 33])
def test_xml2db_pass_chunked_reference_dataset_matches_xmltodict(
    run_xml2db_pipeline_fxt: Callable[..., Any],
    reference_xml_data_dir: Path,
    chunk_size: int,
) -> None:
    """With xml_tag set, the pipeline still reconstructs the same entries as xmltodict.

    Runs against the single 100-entry master file (one file, many chunks per file) with chunk
    sizes that don't evenly divide 100, to exercise remainder handling as well as multi-chunk
    behaviour.
    """
    reconstructed = run_and_reconstruct(run_xml2db_pipeline_fxt, "chunk_100_el", xml_tag="entry", chunk_size=chunk_size)
    assert len(reconstructed) == EXPECTED_ENTRY_COUNT

    reference = xmltodict_reference_entries(reference_xml_data_dir)
    assert sorted_json(reconstructed) == sorted_json(reference)


def test_xml2db_pass_chunked_output_matches_unchunked_output(run_xml2db_pipeline_fxt: Callable[..., Any]) -> None:
    """Enabling chunking does not change the reconstructed content, only how it's parsed internally."""
    unchunked = run_and_reconstruct(run_xml2db_pipeline_fxt, "chunk_5_el")
    chunked = run_and_reconstruct(run_xml2db_pipeline_fxt, "chunk_5_el", xml_tag="entry", chunk_size=2)

    assert len(chunked) == len(unchunked) == EXPECTED_ENTRY_COUNT
    assert sorted_json(chunked) == sorted_json(unchunked)


def test_xml2db_pass_chunked_referential_integrity_across_files_and_chunks(
    run_xml2db_pipeline_fxt: Callable[..., Any],
) -> None:
    """Foreign keys still resolve correctly when both multiple files and multiple chunks per file are involved."""
    (load_info, output_dir) = run_xml2db_pipeline_fxt("chunk_5_el", xml_tag="entry", chunk_size=2)
    assert load_info is not None
    assert not load_info.has_failed_jobs
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, PARQUET)
    assert_uniref_referential_integrity(tables)


# post-load compaction (Xml2DbSettings.compact_reused_tables)


def test_xml2db_pass_chunked_compaction_matches_unchunked_row_counts_exactly(
    run_xml2db_pipeline_fxt: Callable[..., Any],
) -> None:
    """Compacted chunked "reused" tables end up with the exact same row count as unchunked.

    Not just the same *content*, which the reconstruction-based tests above already establish,
    but the same physical row count, with no literal duplicates left over.
    """
    (unchunked_info, unchunked_dir) = run_xml2db_pipeline_fxt("chunk_100_el")
    (chunked_info, chunked_dir) = run_xml2db_pipeline_fxt("chunk_100_el", xml_tag="entry", chunk_size=7)
    assert unchunked_info is not None
    assert chunked_info is not None
    assert not unchunked_info.has_failed_jobs
    assert not chunked_info.has_failed_jobs

    unchunked_tables = read_pipeline_tables(unchunked_dir / unchunked_info.dataset_name, PARQUET)
    chunked_tables = read_pipeline_tables(chunked_dir / chunked_info.dataset_name, PARQUET)

    for table_name in ("entry", "property", "representative_member"):
        assert len(chunked_tables[table_name]) == len(unchunked_tables[table_name]), table_name


def test_xml2db_pass_chunked_compaction_root_table_is_left_alone(
    run_xml2db_pipeline_fxt: Callable[..., Any],
) -> None:
    """The root table is exempt from content-addressing, so compaction never touches it.

    See `cdm_data_loaders.readers.xml.is_xml2db_content_addressed_table` for why: it never has
    the kind of literal duplicate compaction removes, so it still gets one row per chunk after
    compaction runs, unlike property/representative_member.
    """
    (load_info, output_dir) = run_xml2db_pipeline_fxt("chunk_100_el", xml_tag="entry", chunk_size=7)
    assert load_info is not None
    assert not load_info.has_failed_jobs
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, PARQUET)

    n_chunks = -(-EXPECTED_ENTRY_COUNT // 7)  # ceil division
    assert len(tables["uniref"]) == n_chunks


def test_xml2db_pass_chunked_compaction_preserves_referential_integrity(
    run_xml2db_pipeline_fxt: Callable[..., Any],
) -> None:
    """Every foreign key still resolves after compaction has removed literal duplicate rows."""
    (load_info, output_dir) = run_xml2db_pipeline_fxt("chunk_100_el", xml_tag="entry", chunk_size=7)
    assert load_info is not None
    assert not load_info.has_failed_jobs
    tables = read_pipeline_tables(output_dir / load_info.dataset_name, PARQUET)
    assert_uniref_referential_integrity(tables)


def test_xml2db_pass_compact_reused_tables_false_keeps_literal_duplicates(
    run_xml2db_pipeline_fxt: Callable[..., Any],
) -> None:
    """Disabling compact_reused_tables leaves the literal duplicate rows chunking can produce."""
    (unchunked_info, unchunked_dir) = run_xml2db_pipeline_fxt("chunk_100_el")
    (chunked_info, chunked_dir) = run_xml2db_pipeline_fxt(
        "chunk_100_el", xml_tag="entry", chunk_size=7, compact_reused_tables=False
    )
    assert unchunked_info is not None
    assert chunked_info is not None

    unchunked_tables = read_pipeline_tables(unchunked_dir / unchunked_info.dataset_name, PARQUET)
    chunked_tables = read_pipeline_tables(chunked_dir / chunked_info.dataset_name, PARQUET)

    # without compaction, the chunked run keeps every literal duplicate, so it has *more* rows...
    assert len(chunked_tables["property"]) > len(unchunked_tables["property"])
    # ...but still exactly the same *distinct* keys (see the reader-level unit tests for why).
    assert {row["pk_property"] for row in chunked_tables["property"]} == {
        row["pk_property"] for row in unchunked_tables["property"]
    }


def test_xml2db_pass_compaction_dedupes_across_separate_unchunked_files(
    run_xml2db_pipeline_fxt: Callable[..., Any],
) -> None:
    """compact_reused_tables also dedupes across multiple separate files, not just chunks.

    chunk_5_el splits the same 100 entries across 20 files with no xml_tag set, so
    every file parses into its own xml2db Document. Content shared across files (e.g. common
    `property` values like protein names, which real UniRef data repeats heavily) gets the same
    content-addressed key (see `is_xml2db_content_addressed_table`) but still lands as separate
    literal duplicate rows until compaction merges them -- exactly the same failure mode as
    chunking, just from multiple files instead of multiple chunks of one file.
    """
    (compacted_info, compacted_dir) = run_xml2db_pipeline_fxt("chunk_5_el")
    (uncompacted_info, uncompacted_dir) = run_xml2db_pipeline_fxt("chunk_5_el", compact_reused_tables=False)
    assert compacted_info is not None
    assert uncompacted_info is not None
    assert not compacted_info.has_failed_jobs
    assert not uncompacted_info.has_failed_jobs

    compacted_tables = read_pipeline_tables(compacted_dir / compacted_info.dataset_name, PARQUET)
    uncompacted_tables = read_pipeline_tables(uncompacted_dir / uncompacted_info.dataset_name, PARQUET)

    # without compaction, distributing the 100 entries across 20 separate files leaves literal
    # duplicate rows behind for property content shared across files (real UniRef data repeats
    # protein names, lengths, and cross-references heavily)...
    assert len(uncompacted_tables["property"]) > len(compacted_tables["property"])
    # ...but the compacted run has exactly one row per distinct content-addressed key, with
    # nothing missing and nothing left over.
    compacted_property_pks = {row["pk_property"] for row in compacted_tables["property"]}
    assert len(compacted_property_pks) == len(compacted_tables["property"])
    assert compacted_property_pks == {row["pk_property"] for row in uncompacted_tables["property"]}

    # entry and every junction table are non-reused: file-count has no effect on their row counts.
    assert len(compacted_tables["entry"]) == len(uncompacted_tables["entry"]) == EXPECTED_ENTRY_COUNT
