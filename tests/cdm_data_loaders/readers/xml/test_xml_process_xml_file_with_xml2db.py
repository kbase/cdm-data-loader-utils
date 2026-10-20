"""Unit tests for the xml2db-based transformation in `cdm_data_loaders.readers.xml`.

These exercise `build_xml2db_model`, `flatten_xml2db_document`, and
`process_xml_file_with_xml2db` against the UniRef reference fixtures (`uniref.xsd` and the
`chunk_5_el` XML files) used by the xmltodict-based reference tests.
"""

import gzip
import json
import math
from collections.abc import Iterable
from pathlib import Path
from typing import Any, Final
from unittest.mock import MagicMock

import pytest
import xmltodict
from dlt.extract.items import DataItemWithMeta, TableNameMeta
from lxml.etree import fromstring
from xml2db import DataModel, Document

from cdm_data_loaders.pipelines.xml2db.settings import Xml2DbSettings
from cdm_data_loaders.readers.xml import DEFAULT_XMLTODICT_ARGS
from cdm_data_loaders.readers.xml2db_doc import (
    build_xml2db_model,
    flatten_xml2db_document,
    iter_xml2db_chunk_fragments,
    process_xml_file_with_xml2db,
)

REFERENCE_XML_DIR: Final[Path] = Path("tests") / "data" / "uniprot" / "uniref"
UNIREF_XSD: Final[Path] = REFERENCE_XML_DIR / "uniref.xsd"
CHUNK_5_EL_DIR: Final[Path] = REFERENCE_XML_DIR / "chunk_5_el"
PART_01: Final[Path] = CHUNK_5_EL_DIR / "uniref_100_part_01.xml.gz"
PART_02: Final[Path] = CHUNK_5_EL_DIR / "uniref_100_part_02.xml.gz"
MASTER_100: Final[Path] = REFERENCE_XML_DIR / "chunk_100_el" / "uniref_100_part_01.xml.gz"
ENTRIES_PER_PART: Final[int] = 5
ENTRIES_IN_MASTER: Final[int] = 100
ENTRY_LOCAL_NAME: Final[str] = "entry"
ENTRY_EXPANDED_TAG: Final[str] = "{http://uniprot.org/uniref}entry"

# Table/type names as xml2db derives them from uniref.xsd.
ENTRY_TABLE: Final[str] = "entry"
PROPERTY_TABLE: Final[str] = "property"
MEMBER_TABLE: Final[str] = "representativeMember"
ENTRY_PROPERTY_JUNCTION: Final[str] = "entry_property"
ENTRY_MEMBER_JUNCTION: Final[str] = "entry_member_representativeMember"
MEMBER_PROPERTY_JUNCTION: Final[str] = "representativeMember_dbReference_property"

RECORD_HASH_COLUMN: Final[str] = "xml2db_record_hash"
# sha1 (xml2db's default hash) is 20 bytes -> 40 hex characters.
SHA1_HEX_LENGTH: Final[int] = 40


def _table_and_data(items: Iterable[DataItemWithMeta]) -> list[tuple[str, Any]]:
    """Flatten a sequence of DataItemWithMeta into (table_name, rows) pairs."""
    result: list[tuple[str, Any]] = []
    for item in items:
        assert isinstance(item, DataItemWithMeta)
        assert isinstance(item.meta, TableNameMeta)
        result.append((item.meta.table_name, item.data))
    return result


def _rows_by_table(items: Iterable[DataItemWithMeta]) -> dict[str, list[dict[str, Any]]]:
    """Merge all pages for the same table into a single list of rows, keyed by table name."""
    merged: dict[str, list[dict[str, Any]]] = {}
    for table_name, rows in _table_and_data(items):
        merged.setdefault(table_name, []).extend(rows)
    return merged


def fake_xml2db_settings(
    buffer_size: int = 100,
    skip_xml_validation: bool = True,
    xml_tag: str | None = None,
    chunk_size: int = 1000,
) -> Xml2DbSettings:
    """Build a MagicMock stand-in for Xml2DbSettings exposing just the attributes the reader uses."""
    settings = MagicMock(spec=Xml2DbSettings)
    settings.buffer_size = buffer_size
    settings.skip_xml_validation = skip_xml_validation
    settings.xml_tag = xml_tag
    settings.chunk_size = chunk_size
    return settings


@pytest.fixture(scope="module")
def uniref_model() -> DataModel:
    """Build the xml2db DataModel from uniref.xsd once for the whole test module."""
    return build_xml2db_model(UNIREF_XSD, short_name="uniref_test")


def _reference_entries(file_path: Path) -> dict[str, dict[str, Any]]:
    """Parse a UniRef file with xmltodict directly, keyed by entry id, for cross-checking."""
    with gzip.open(file_path, "rb") as fh:
        document = xmltodict.parse(fh.read(), **DEFAULT_XMLTODICT_ARGS)
    entries = document["UniRef50"]["entry"]
    if isinstance(entries, dict):
        entries = [entries]
    return {entry["_id"]: entry for entry in entries}


def test_build_xml2db_model_pass_derives_expected_tables(uniref_model: DataModel) -> None:
    """The model built from uniref.xsd has one table per distinct XSD type, named by element."""
    table_names = {table.name for table in uniref_model.tables.values()}
    assert {ENTRY_TABLE, PROPERTY_TABLE, MEMBER_TABLE, "uniref_test"} <= table_names
    # "member" and "representativeMember" share memberType, so they collapse into one table.
    assert "member" not in table_names


def test_process_xml_file_with_xml2db_pass_entry_rows_match_xmltodict_reference(
    uniref_model: DataModel,
) -> None:
    """Every entry row's id/name/updated matches the same entry parsed directly by xmltodict."""
    reference = _reference_entries(PART_01)
    assert len(reference) == ENTRIES_PER_PART

    items = list(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_01))
    rows = _rows_by_table(items)

    entry_rows = rows[ENTRY_TABLE]
    assert len(entry_rows) == ENTRIES_PER_PART
    assert {row["id"] for row in entry_rows} == set(reference)

    for row in entry_rows:
        ref_entry = reference[row["id"]]
        assert row["name"] == ref_entry["name"]
        assert row["updated"] == ref_entry["_updated"]
        # entry and representativeMember are both "reused" tables in this schema, so their keys
        # are content-addressed (see _content_addressed_key), not namespaced by source file.
        assert row["pk_entry"] == f"{ENTRY_TABLE}:{row[RECORD_HASH_COLUMN]}"
        assert row["fk_representativeMember"].startswith(f"{MEMBER_TABLE}:")


def test_process_xml_file_with_xml2db_pass_record_hash_is_hex_string(uniref_model: DataModel) -> None:
    """xml2db's raw-bytes record hash is converted to a hex string for downstream serialisation."""
    items = list(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_01))
    rows = _rows_by_table(items)

    for table_name in (ENTRY_TABLE, PROPERTY_TABLE, MEMBER_TABLE):
        for row in rows[table_name]:
            record_hash = row[RECORD_HASH_COLUMN]
            assert isinstance(record_hash, str)
            assert len(record_hash) == SHA1_HEX_LENGTH
            int(record_hash, 16)  # raises ValueError if not valid hex


def test_process_xml_file_with_xml2db_pass_referential_integrity(uniref_model: DataModel) -> None:
    """Every foreign key produced by the reader points at an existing primary key in its table."""
    items = list(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_01))
    rows = _rows_by_table(items)

    entry_pks = {row["pk_entry"] for row in rows[ENTRY_TABLE]}
    property_pks = {row["pk_property"] for row in rows[PROPERTY_TABLE]}
    member_pks = {row["pk_representativeMember"] for row in rows[MEMBER_TABLE]}

    # rel1: every entry's representative member exists.
    assert {row["fk_representativeMember"] for row in rows[ENTRY_TABLE]} <= member_pks

    # n-n junctions: both sides of every relation exist in their own tables.
    for row in rows[ENTRY_PROPERTY_JUNCTION]:
        assert row["fk_entry"] in entry_pks
        assert row["fk_property"] in property_pks

    for row in rows.get(ENTRY_MEMBER_JUNCTION, []):
        assert row["fk_entry"] in entry_pks
        assert row["fk_representativeMember"] in member_pks

    for row in rows[MEMBER_PROPERTY_JUNCTION]:
        assert row["fk_representativeMember"] in member_pks
        assert row["fk_property"] in property_pks


def test_process_xml_file_with_xml2db_pass_keys_are_unique_across_files(uniref_model: DataModel) -> None:
    """Distinct entries across two files get distinct keys, even though entry is a reused table.

    entry's key is content-addressed like any other reused table's (see
    `test_process_xml_file_with_xml2db_pass_entry_rows_match_xmltodict_reference`), but every
    entry's unique `id` attribute is part of that content, so two files' *distinct* entries never
    collide in practice. See
    `test_process_xml_file_with_xml2db_pass_reused_table_keys_collide_across_files_when_content_matches`
    for the same key scheme deliberately colliding two files' rows that share identical content.
    """
    items_1 = list(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_01))
    items_2 = list(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_02))

    entry_pks_1 = {row["pk_entry"] for row in _rows_by_table(items_1)[ENTRY_TABLE]}
    entry_pks_2 = {row["pk_entry"] for row in _rows_by_table(items_2)[ENTRY_TABLE]}

    assert len(entry_pks_1) == ENTRIES_PER_PART
    assert len(entry_pks_2) == ENTRIES_PER_PART
    assert entry_pks_1.isdisjoint(entry_pks_2)


def test_process_xml_file_with_xml2db_pass_reused_table_keys_collide_across_files_when_content_matches(
    uniref_model: DataModel,
) -> None:
    """Identical content in two different files now resolves to the very same key.

    Every UniRef representative member carries an `isSeed=true` property, present in both
    PART_01 and PART_02; under the old file-namespaced key scheme these could never share a key
    even though their content is identical. Content-addressing fixes that: this is precisely what
    makes the rows produced from two different files (or, in chunked mode, two different chunks
    of one file) safe to collapse later without touching any foreign key -- see
    `cdm_data_loaders.pipelines.xml2db.compaction`.
    """
    rows_1 = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_01))
    rows_2 = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_02))

    def _seed_property_pk(rows: dict[str, list[dict[str, Any]]]) -> str:
        for row in rows[PROPERTY_TABLE]:
            if row["type"] == "isSeed" and row["value"] == "true":
                return row["pk_property"]
        pytest.fail("no isSeed=true property found in fixture data")
        return ""  # unreachable; keeps type-checkers happy

    assert _seed_property_pk(rows_1) == _seed_property_pk(rows_2)


def test_process_xml_file_with_xml2db_pass_gzip_and_plain_produce_equivalent_rows(
    uniref_model: DataModel, tmp_path: Path
) -> None:
    """A gzip-compressed file and its decompressed copy parse to the same content."""
    plain_path = tmp_path / "uniref_100_part_01.xml"
    with gzip.open(PART_01, "rb") as fh:
        plain_path.write_bytes(fh.read())

    gz_rows = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, PART_01))
    plain_rows = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, plain_path))

    def _strip_file_prefix(value: Any, file_name: str) -> Any:  # noqa: ANN401
        if isinstance(value, str) and value.startswith(f"{file_name}:"):
            return value.removeprefix(f"{file_name}:")
        return value

    def _normalize(rows: list[dict[str, Any]], file_name: str) -> list[dict[str, Any]]:
        return sorted(
            ({key: _strip_file_prefix(value, file_name) for key, value in row.items()} for row in rows),
            key=lambda r: str(r.get("pk_entry") or r.get("pk_property") or r.get("pk_representativeMember") or r),
        )

    assert set(gz_rows) == set(plain_rows)
    for table_name in gz_rows:
        assert _normalize(gz_rows[table_name], PART_01.name) == _normalize(plain_rows[table_name], plain_path.name)


def test_process_xml_file_with_xml2db_pass_buffer_size_does_not_change_total_content(
    uniref_model: DataModel,
) -> None:
    """Every row is still produced regardless of buffer_size, only the number of pages may differ."""
    small_buffer_items = list(process_xml_file_with_xml2db(fake_xml2db_settings(buffer_size=1), uniref_model, PART_01))
    large_buffer_items = list(
        process_xml_file_with_xml2db(fake_xml2db_settings(buffer_size=100_000), uniref_model, PART_01)
    )

    small_rows = _rows_by_table(small_buffer_items)
    large_rows = _rows_by_table(large_buffer_items)

    assert set(small_rows) == set(large_rows)
    for table_name in small_rows:
        assert sorted(small_rows[table_name], key=str) == sorted(large_rows[table_name], key=str)


def test_flatten_xml2db_document_pass_unknown_table_key_maps_via_model_name(uniref_model: DataModel) -> None:
    """flatten_xml2db_document keys output by table .name, not by the internal XSD type name."""
    document = Document(uniref_model)
    with gzip.open(PART_01, "rb") as fh:
        document.parse_xml(fh, skip_validation=True, iterparse=True)

    tables = flatten_xml2db_document(uniref_model, document, file_key=PART_01.name)

    # "entryType" is the XSD type name; "entry" is the table name that should appear instead.
    assert "entryType" not in tables
    assert ENTRY_TABLE in tables
    assert len(tables[ENTRY_TABLE]) == ENTRIES_PER_PART


# iter_xml2db_chunk_fragments


def _entry_count(fragment: bytes) -> int:
    """Count local-name `entry` children in a serialized fragment, ignoring namespaces."""
    root = fromstring(fragment)
    return sum(1 for child in root if child.tag.rsplit("}", 1)[-1] == ENTRY_LOCAL_NAME)


@pytest.mark.parametrize("chunk_size", [1, 2, 3, 5, 100])
def test_iter_xml2db_chunk_fragments_pass_splits_into_expected_chunk_sizes(chunk_size: int) -> None:
    """Fragments cover every entry exactly once, in order, in chunks of at most `chunk_size`."""
    fragments = list(iter_xml2db_chunk_fragments(PART_01, ENTRY_LOCAL_NAME, chunk_size))
    sizes = [_entry_count(fragment) for fragment in fragments]

    assert sum(sizes) == ENTRIES_PER_PART
    assert len(fragments) == math.ceil(ENTRIES_PER_PART / chunk_size)
    assert all(size <= chunk_size for size in sizes)
    assert all(size == chunk_size for size in sizes[:-1])  # only the last fragment may be a remainder


@pytest.mark.parametrize("child_tag", [ENTRY_LOCAL_NAME, ENTRY_EXPANDED_TAG], ids=["bare-local-name", "namespaced-tag"])
def test_iter_xml2db_chunk_fragments_pass_bare_and_namespaced_tag_are_equivalent(child_tag: str) -> None:
    """A bare local name and its Clark-notation expanded form select the same elements."""
    fragments = list(iter_xml2db_chunk_fragments(PART_01, child_tag, chunk_size=2))
    assert sum(_entry_count(fragment) for fragment in fragments) == ENTRIES_PER_PART


def test_iter_xml2db_chunk_fragments_pass_preserves_root_tag_and_attributes() -> None:
    """Every fragment's root is a faithful copy of the source root, tag/attributes included."""
    with gzip.open(PART_01, "rb") as fh:
        original_root = fromstring(fh.read())

    for fragment in iter_xml2db_chunk_fragments(PART_01, ENTRY_LOCAL_NAME, chunk_size=2):
        fragment_root = fromstring(fragment)
        assert fragment_root.tag == original_root.tag
        assert dict(fragment_root.attrib) == dict(original_root.attrib)


def test_iter_xml2db_chunk_fragments_pass_no_matching_children_yields_nothing(tmp_path: Path) -> None:
    """A well-formed document with no elements matching child_tag produces no fragments."""
    xml_path = tmp_path / "empty.xml"
    xml_path.write_text("<root><other/><other/></root>", encoding="utf-8")

    assert list(iter_xml2db_chunk_fragments(xml_path, "missing", chunk_size=10)) == []


def test_iter_xml2db_chunk_fragments_pass_empty_file_yields_nothing(tmp_path: Path) -> None:
    """A zero-byte file (no root element at all) produces no fragments rather than raising."""
    xml_path = tmp_path / "truly_empty.xml"
    xml_path.write_bytes(b"")

    assert list(iter_xml2db_chunk_fragments(xml_path, ENTRY_LOCAL_NAME, chunk_size=10)) == []


#  process_xml_file_with_xml2db, chunked mode


def _data_columns(row: dict[str, Any]) -> dict[str, Any]:
    """Drop pk/fk/hash columns, leaving only a row's own data.

    Chunked parsing uses a fresh xml2db `Document` per chunk, so temp pk numbering restarts for
    every chunk and is *not* comparable to the single-pass numbering used by whole-file parsing
    (unlike, say, gzip vs. plain re-encodings of the same single-pass document, whose pks match
    exactly because both are one Document parsing the same entries in the same order). Content
    equivalence therefore has to be checked on data columns only; row counts are checked
    separately for the pk/fk-only junction tables, which have no data columns of their own.
    """
    return {
        key: value for key, value in row.items() if not key.startswith(("pk_", "fk_")) and key != RECORD_HASH_COLUMN
    }


def _normalize_data_rows(rows: list[dict[str, Any]]) -> list[str]:
    """Sorted JSON strings of each row's data columns, for order-insensitive comparison."""
    return sorted(json.dumps(_data_columns(row), sort_keys=True, default=str) for row in rows)


#: Tables xml2db deduplicates ("reused") by content hash. Whole-file parsing dedups across the
#: entire file (one hashmap); chunked parsing only dedups *within* each chunk (a fresh Document,
#: and therefore a fresh hashmap, per chunk), so these tables may legitimately end up with *more*
#: rows when chunked (duplicate content that would have collapsed to one row in a single pass can
#: now appear once per chunk it occurs in). This is the documented trade-off of chunking: bounded
#: memory in exchange for weaker (chunk-local, not file-wide) deduplication.
DEDUPLICATED_TABLES: Final[frozenset[str]] = frozenset({PROPERTY_TABLE, MEMBER_TABLE})


def test_process_xml_file_with_xml2db_pass_chunked_matches_unchunked_row_counts(uniref_model: DataModel) -> None:
    """Row counts match exactly for non-deduplicated tables; may only grow (never shrink) for reused ones."""
    unchunked_rows = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, MASTER_100))
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=7), uniref_model, MASTER_100
        )
    )
    # the virtual root table gets one (namespaced, never content-addressed -- see
    # is_xml2db_content_addressed_table) row per chunk, so its row count grows with the number of
    # chunks just like property/representativeMember's, albeit for a different underlying reason.
    # `entry` is technically "reused" too, but every one of the 100 fixture entries has a unique
    # `id`, so none of them ever collide and its row count matches exactly either way.
    root_table_name = uniref_model.tables[uniref_model.root_table].name
    deduplicated_tables = DEDUPLICATED_TABLES | {root_table_name}

    assert set(unchunked_rows) == set(chunked_rows)
    for table_name in unchunked_rows:
        if table_name in deduplicated_tables:
            assert len(chunked_rows[table_name]) >= len(unchunked_rows[table_name]), table_name
        else:
            assert len(chunked_rows[table_name]) == len(unchunked_rows[table_name]), table_name


def test_process_xml_file_with_xml2db_pass_chunked_entry_content_matches_unchunked_exactly(
    uniref_model: DataModel,
) -> None:
    """Entry rows are never deduplicated, so chunking must not change their content at all."""
    unchunked_rows = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, MASTER_100))
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=7), uniref_model, MASTER_100
        )
    )

    assert _normalize_data_rows(unchunked_rows[ENTRY_TABLE]) == _normalize_data_rows(chunked_rows[ENTRY_TABLE])


def test_process_xml_file_with_xml2db_pass_chunked_preserves_same_distinct_values_for_reused_tables(
    uniref_model: DataModel,
) -> None:
    """No distinct property/member value is lost or invented by chunking, only its deduplication."""
    unchunked_rows = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, MASTER_100))
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=7), uniref_model, MASTER_100
        )
    )

    for table_name in DEDUPLICATED_TABLES:
        assert set(_normalize_data_rows(unchunked_rows[table_name])) == set(
            _normalize_data_rows(chunked_rows[table_name])
        )


def test_process_xml_file_with_xml2db_pass_chunked_entry_count_matches(uniref_model: DataModel) -> None:
    """Chunked parsing yields exactly one entry row per source <entry> element, none lost or duplicated."""
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=7), uniref_model, MASTER_100
        )
    )
    assert len(chunked_rows[ENTRY_TABLE]) == ENTRIES_IN_MASTER
    assert len({row["id"] for row in chunked_rows[ENTRY_TABLE]}) == ENTRIES_IN_MASTER


def test_process_xml_file_with_xml2db_pass_chunked_keys_unique_across_chunks(uniref_model: DataModel) -> None:
    """Entry keys are content-addressed, and the 100 fixture entries remain distinct across chunks."""
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=7), uniref_model, MASTER_100
        )
    )
    entry_pks = [row["pk_entry"] for row in chunked_rows[ENTRY_TABLE]]
    assert len(entry_pks) == len(set(entry_pks)) == ENTRIES_IN_MASTER
    for pk in entry_pks:
        assert pk.startswith(f"{ENTRY_TABLE}:")
        assert len(pk) == len(ENTRY_TABLE) + 1 + SHA1_HEX_LENGTH


def test_process_xml_file_with_xml2db_pass_chunked_reused_table_distinct_key_count_matches_unchunked(
    uniref_model: DataModel,
) -> None:
    """Distinct keys for reused tables match exactly between chunked and unchunked parsing.

    Literal row counts may still differ until the output is compacted (see
    `cdm_data_loaders.pipelines.xml2db.compaction`), but the *set* of keys never does: identical
    content always produces the same content-addressed key, in any `Document`. Since unchunked
    (whole-file, single-`Document`) parsing already fully deduplicates reused tables, this distinct
    key count is also that table's true minimal row count. This excludes the virtual root table --
    see `test_process_xml_file_with_xml2db_pass_chunked_root_table_keys_are_not_content_addressed`.
    """
    unchunked_rows = _rows_by_table(process_xml_file_with_xml2db(fake_xml2db_settings(), uniref_model, MASTER_100))
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=7), uniref_model, MASTER_100
        )
    )

    for table_name in DEDUPLICATED_TABLES:
        pk_column = f"pk_{table_name}"
        unchunked_distinct = {row[pk_column] for row in unchunked_rows[table_name]}
        chunked_distinct = {row[pk_column] for row in chunked_rows[table_name]}
        assert chunked_distinct == unchunked_distinct, table_name
        assert len(unchunked_rows[table_name]) == len(unchunked_distinct), table_name


def test_process_xml_file_with_xml2db_pass_chunked_root_table_keys_are_not_content_addressed(
    uniref_model: DataModel,
) -> None:
    """The virtual root table is exempt from content-addressing (see `is_xml2db_content_addressed_table`).

    Its own `xml2db_record_hash` incorporates the hash of every entry attached to it; since a
    chunk's synthetic root only ever has *that chunk's* entries attached, every chunk's root row
    gets a different hash even though the underlying `releaseDate`/`version` values are identical
    -- so, unlike `property`/`representativeMember`, root rows never collide across chunks, and
    `compact_reused_tables` deliberately leaves this table alone (see
    `cdm_data_loaders.pipelines.xml2db.compaction.reused_table_names`).
    """
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=7), uniref_model, MASTER_100
        )
    )
    root_table_name = uniref_model.tables[uniref_model.root_table].name
    root_pks = [row[f"pk_{root_table_name}"] for row in chunked_rows[root_table_name]]

    # one root row per chunk, all distinct, none of them content-addressed (no ":" + 40 hex chars).
    assert len(root_pks) == len(set(root_pks)) > 1
    assert all(not pk.startswith(f"{root_table_name}:") for pk in root_pks)


def test_process_xml_file_with_xml2db_pass_chunk_size_one_still_produces_all_rows(uniref_model: DataModel) -> None:
    """The degenerate case of one entry per chunk (many tiny Documents) still yields correct output."""
    chunked_rows = _rows_by_table(
        process_xml_file_with_xml2db(
            fake_xml2db_settings(xml_tag=ENTRY_LOCAL_NAME, chunk_size=1), uniref_model, PART_01
        )
    )
    assert len(chunked_rows[ENTRY_TABLE]) == ENTRIES_PER_PART
