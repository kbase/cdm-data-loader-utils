"""Reconstruction helpers for the xml2db-based UniRef reference tests.

xml2db (unlike xmltodict) produces a normalized relational model: one row per entry/member/
property plus many-to-many junction tables, with primary/foreign keys namespaced per source file
by `process_xml_file_with_xml2db`. dlt further normalizes table and column names to snake_case
on load. These helpers join everything back into one nested dict per entry, comparable across
pipeline runs and against the plain xmltodict-based reference.
"""

from collections import defaultdict
from functools import partial
from json import dumps
from typing import Any

from tests.cdm_data_loaders.pipelines.helpers import strip_metadata_columns

ENTRY_TABLE = "entry"
PROPERTY_TABLE = "property"
MEMBER_TABLE = "representative_member"
ENTRY_PROPERTY_JUNCTION = "entry_property"
ENTRY_MEMBER_JUNCTION = "entry_member_representative_member"
MEMBER_PROPERTY_JUNCTION = "representative_member_db_reference_property"
EXPECTED_ENTRY_COUNT = 100


def _index_by(rows: list[dict[str, Any]], key: str) -> dict[Any, dict[str, Any]]:
    """Index a list of rows by a unique column value."""
    return {row[key]: row for row in rows}


def _group_by(rows: list[dict[str, Any]], key: str) -> dict[Any, list[dict[str, Any]]]:
    """Group a list of rows by a (non-unique) column value."""
    grouped: dict[Any, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[row[key]].append(row)
    return grouped


def _property_pairs(
    member_or_entry_pk: Any,
    junction_rows: dict[Any, list[dict[str, Any]]],
    properties_by_pk: dict[Any, dict[str, Any]],
) -> list[list[str]]:
    """Resolve a set of junction-table links into sorted [type, value] property pairs."""
    pairs = [
        [properties_by_pk[link["fk_property"]]["type"], properties_by_pk[link["fk_property"]]["value"]]
        for link in junction_rows.get(member_or_entry_pk, [])
    ]
    return sorted(pairs)


def _member_summary(
    member_row: dict[str, Any],
    member_property_junction: dict[Any, list[dict[str, Any]]],
    properties_by_pk: dict[Any, dict[str, Any]],
) -> dict[str, Any]:
    """Build a comparable summary of a (representative or plain) member row."""
    return {
        "dbReference_type": member_row.get("db_reference_type"),
        "dbReference_id": member_row.get("db_reference_id"),
        "sequence_length": member_row.get("sequence_length"),
        "sequence_checksum": member_row.get("sequence_checksum"),
        "sequence_value": member_row.get("sequence_value"),
        "properties": _property_pairs(
            member_row["pk_representative_member"], member_property_junction, properties_by_pk
        ),
    }


def reconstruct_xml2db_entries(tables: dict[str, list[dict[str, Any]]]) -> list[dict[str, Any]]:
    """Rebuild one nested dict per UniRef entry from the xml2db pipeline's relational tables.

    :param tables: table name -> list of row dicts, as read from a pipeline's output dataset.
    :type tables: dict[str, list[dict[str, Any]]]
    :return: one dict per entry with id/name/properties/representative_member/members, comparable
        across pipeline runs and (after normalizing key ordering) against an xmltodict reference.
    :rtype: list[dict[str, Any]]
    """
    entry_rows = [strip_metadata_columns(r) for r in tables.get(ENTRY_TABLE, [])]
    property_rows = [strip_metadata_columns(r) for r in tables.get(PROPERTY_TABLE, [])]
    member_rows = [strip_metadata_columns(r) for r in tables.get(MEMBER_TABLE, [])]
    entry_property_rows = [strip_metadata_columns(r) for r in tables.get(ENTRY_PROPERTY_JUNCTION, [])]
    entry_member_rows = [strip_metadata_columns(r) for r in tables.get(ENTRY_MEMBER_JUNCTION, [])]
    member_property_rows = [strip_metadata_columns(r) for r in tables.get(MEMBER_PROPERTY_JUNCTION, [])]

    properties_by_pk = _index_by(property_rows, "pk_property")
    members_by_pk = _index_by(member_rows, "pk_representative_member")
    entry_property_junction = _group_by(entry_property_rows, "fk_entry")
    entry_member_junction = _group_by(entry_member_rows, "fk_entry")
    member_property_junction = _group_by(member_property_rows, "fk_representative_member")

    entries = []
    for entry in entry_rows:
        representative_member = members_by_pk[entry["fk_representative_member"]]
        members = [
            _member_summary(members_by_pk[link["fk_representative_member"]], member_property_junction, properties_by_pk)
            for link in entry_member_junction.get(entry["pk_entry"], [])
        ]
        entries.append(
            {
                "id": entry["id"],
                "updated": entry["updated"],
                "name": entry["name"],
                "properties": _property_pairs(entry["pk_entry"], entry_property_junction, properties_by_pk),
                "representative_member": _member_summary(
                    representative_member, member_property_junction, properties_by_pk
                ),
                "members": sorted(members, key=partial(dumps, sort_keys=True, default=str)),
            }
        )
    return entries


def assert_uniref_referential_integrity(tables: dict[str, list[dict[str, Any]]]) -> None:
    """Every foreign key in a merged xml2db dataset resolves to a row in its target table."""
    entry_pks = {row["pk_entry"] for row in tables["entry"]}
    property_pks = {row["pk_property"] for row in tables["property"]}
    member_pks = {row["pk_representative_member"] for row in tables["representative_member"]}

    assert len(entry_pks) == len(tables["entry"]) == EXPECTED_ENTRY_COUNT
    assert {row["fk_representative_member"] for row in tables["entry"]} <= member_pks

    for row in tables["entry_property"]:
        assert row["fk_entry"] in entry_pks
        assert row["fk_property"] in property_pks

    for row in tables.get("entry_member_representative_member", []):
        assert row["fk_entry"] in entry_pks
        assert row["fk_representative_member"] in member_pks

    for row in tables["representative_member_db_reference_property"]:
        assert row["fk_representative_member"] in member_pks
        assert row["fk_property"] in property_pks
