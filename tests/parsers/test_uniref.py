"""Tests for the UniRef importer."""

import os
import sys

sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import gzip
import tempfile
import textwrap
import xml.etree.ElementTree as ET
from datetime import datetime
from types import SimpleNamespace

import pytest

from cdm_data_loader_utils.parsers.uniref import (
    add_cluster_members,
    cdm_entity_id,
    extract_cluster,
    extract_cross_refs,
    get_accession_and_seed,
    get_timestamps,
    load_existing_created,
    parse_uniref_entry,
    parse_uniref_xml,
)

NS = {"ns": "http://uniprot.org/uniref"}


class FakeSparkDF:
    """A fake DataFrame returned by spark.read.format().load().select()."""

    def __init__(self, rows) -> None:
        self._rows = rows

    def collect(self):
        return self._rows


class FakeSparkReader:
    """Mock spark.read.format('delta').load().select() chain."""

    def __init__(self, rows=None, fail=False) -> None:
        self._rows = rows
        self._fail = fail

    def format(self, fmt):
        assert fmt == "delta"
        return self

    def load(self, path):
        if self._fail:
            msg = "Table does not exist"
            raise Exception(msg)
        return self

    def select(self, *cols):
        return FakeSparkDF(self._rows)


@pytest.mark.parametrize(
    ("entity_table", "expected"),
    [
        (None, {}),  # no path
        ("", {}),  # empty path
    ],
)
def test_load_existing_created_no_path(entity_table, expected) -> None:
    """Should return empty dict when entity_table path is missing."""
    fake_spark = SimpleNamespace()
    assert load_existing_created(fake_spark, entity_table) == expected


def test_load_existing_created_success(monkeypatch) -> None:
    """Delta table exists: should return dict of id → created timestamp."""
    rows = [
        {"data_source_entity_id": "UniRef100_A", "created": "2024-01-01T00:00:00"},
        {"data_source_entity_id": "UniRef100_B", "created": "2024-01-02T00:00:00"},
    ]

    fake_reader = FakeSparkReader(rows=rows)

    # Patch spark.read to our fake reader
    fake_spark = SimpleNamespace(read=fake_reader)

    result = load_existing_created(fake_spark, "/fake/path/entity")

    assert result == {
        "UniRef100_A": "2024-01-01T00:00:00",
        "UniRef100_B": "2024-01-02T00:00:00",
    }


def test_load_existing_created_missing_table(monkeypatch) -> None:
    """If Delta table does not exist (load fails), return empty dict."""
    fake_reader = FakeSparkReader(fail=True)

    fake_spark = SimpleNamespace(read=fake_reader)

    result = load_existing_created(fake_spark, "/fake/path/entity")

    assert result == {}


@pytest.mark.parametrize(
    ("accession", "expected_prefix"),
    [
        ("A0B0123456", "CDM:"),
        ("P01234", "CDM:"),
        ("", None),
        (None, None),
    ],
)
def test_cdm_entity_id(accession, expected_prefix) -> None:
    result = cdm_entity_id(accession)
    if expected_prefix is None:
        assert result is None
    else:
        assert result.startswith(expected_prefix)


@pytest.mark.parametrize(
    ("xml_str", "expected_name"),
    [
        (
            "<entry xmlns='http://uniprot.org/uniref' id='UniRef100_A0A009GP46' updated='2016-10-05'>"
            "<name>TestName</name></entry>",
            "TestName",
        ),
        (
            "<entry xmlns='http://uniprot.org/uniref' id='UniRef100_XYZ' updated='2024-01-01'/>",
            "UNKNOWN",
        ),
    ],
)
def test_extract_cluster(xml_str, expected_name) -> None:
    ns = {"ns": "http://uniprot.org/uniref"}
    elem = ET.fromstring(xml_str)

    uniref_id = elem.attrib.get("id")
    cluster_id, name = extract_cluster(elem, ns, uniref_id)

    assert cluster_id.startswith("cdm_ccol_")
    assert isinstance(cluster_id, str)
    assert name == expected_name


@pytest.mark.parametrize(
    ("uniref_id", "existing_created", "now", "expected"),
    [
        # Has existing_created
        (
            "UniRef100_A",
            {"UniRef100_A": "2024-01-01T00:00:00"},
            datetime(2025, 1, 1, 0, 0, 0),
            ("2025-01-01T00:00:00", "2024-01-01T00:00:00"),
        ),
        # There is no existing_created
        (
            "UniRef100_B",
            {"UniRef100_A": "2024-01-01T00:00:00"},
            datetime(2025, 1, 1, 0, 0, 0),
            ("2025-01-01T00:00:00", "2025-01-01T00:00:00"),
        ),
        # There is no existing_created，also not provide "now"
        (
            "UniRef100_C",
            {},
            None,  # The system automatically use the current time
            None,  # Only assert that the return is a string and they are equal
        ),
    ],
)
def test_get_timestamps(uniref_id, existing_created, now, expected) -> None:
    result = get_timestamps(uniref_id, existing_created, now)
    if expected is not None:
        assert result == expected
    else:
        formatted_now, created_time = result
        assert formatted_now == created_time
        assert isinstance(formatted_now, str)
        assert len(formatted_now) == 19  # "YYYY-MM-DDTHH:MM:SS"


@pytest.mark.parametrize(
    ("xml_str", "expected_acc", "expected_is_seed"),
    [
        # Have accession and isSeed
        (
            """
        <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="A0A009HJL9_ACIB9">
            <property type="UniProtKB accession" value="A0A009HJL9"/>
            <property type="isSeed" value="true"/>
        </dbReference>
        """,
            "A0A009HJL9",
            True,
        ),
        # Only accession，No isSeed
        (
            """
        <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="A0A241V597_9GAMM">
            <property type="UniProtKB accession" value="A0A241V597"/>
        </dbReference>
        """,
            "A0A241V597",
            False,
        ),
        # No accession，only id
        (
            """
        <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="ID_ONLY"></dbReference>
        """,
            "ID_ONLY",
            False,
        ),
        # None
        (
            None,
            None,
            False,
        ),
    ],
)
def test_get_accession_and_seed(xml_str, expected_acc, expected_is_seed) -> None:
    ns = {"ns": "http://uniprot.org/uniref"}
    dbref = ET.fromstring(xml_str) if xml_str else None
    acc, is_seed = get_accession_and_seed(dbref, ns)
    assert acc == expected_acc
    assert is_seed == expected_is_seed


def make_entry_with_members(member_xmls, ns_uri="http://uniprot.org/uniref"):
    """
    Receives a list of xml strings from dbReference,
    generates an <entry> element with <member> child nodes.
    """
    entry_elem = ET.Element(f"{{{ns_uri}}}entry")
    for dbref_xml in member_xmls:
        dbref_elem = ET.fromstring(dbref_xml)
        member_elem = ET.SubElement(entry_elem, f"{{{ns_uri}}}member")
        member_elem.append(dbref_elem)
    return entry_elem


@pytest.mark.parametrize(
    ("repr_xml", "member_xmls", "expected"),
    [
        pytest.param(
            # representative member
            textwrap.dedent("""
                <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="REP_ID">
                    <property type="UniProtKB accession" value="REP_ACC"/>
                    <property type="isSeed" value="true"/>
                </dbReference>
            """),
            [
                textwrap.dedent("""
                    <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="MEM1_ID">
                        <property type="UniProtKB accession" value="MEM1_ACC"/>
                    </dbReference>
                """),
                textwrap.dedent("""
                    <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="MEM2_ID">
                        <property type="UniProtKB accession" value="MEM2_ACC"/>
                        <property type="isSeed" value="true"/>
                    </dbReference>
                """),
            ],
            [
                # cdm_prot_
                ("CLUSTER_X", "cdm_prot_", "true", "true", "1.0"),
                ("CLUSTER_X", "cdm_prot_", "false", "false", "1.0"),
                ("CLUSTER_X", "cdm_prot_", "false", "true", "1.0"),
            ],
            id="with-representative-and-members",
        ),
        pytest.param(
            # Only members, no representative member
            None,
            [
                textwrap.dedent("""
                    <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="MEM_ID">
                        <property type="UniProtKB accession" value="MEM_ACC"/>
                    </dbReference>
                """),
            ],
            [
                # cdm_prot_
                ("CLUSTER_X", "cdm_prot_", "false", "false", "1.0"),
            ],
            id="members-only",
        ),
        pytest.param(
            # No members, no representative member
            None,
            [],
            [],
            id="no-members",
        ),
    ],
)
def test_add_cluster_members(repr_xml, member_xmls, expected) -> None:
    """Test add_cluster_members with various representative/member combinations."""
    ns = {"ns": "http://uniprot.org/uniref"}
    cluster_id = "CLUSTER_X"

    # Structure (representative members) dbReference if it exists
    repr_db = ET.fromstring(repr_xml) if repr_xml else None

    # Structure <entry> nodes, and add <member>
    elem = make_entry_with_members(member_xmls)

    # Calling the function under test
    cluster_member_data = []
    add_cluster_members(cluster_id, repr_db, elem, cluster_member_data, ns)

    assert len(cluster_member_data) == len(expected)
    for i, (clu_id, cdm_prefix, is_repr, is_seed, score) in enumerate(expected):
        out = cluster_member_data[i]
        assert out[0] == clu_id, f"Wrong cluster_id at idx {i}: {out[0]}"
        assert out[1].startswith(cdm_prefix), f"Wrong entity_id at idx {i}: {out[1]}"
        assert out[2] == is_repr, f"Wrong is_representative at idx {i}: {out[2]}"
        assert out[3] == is_seed, f"Wrong is_seed at idx {i}: {out[3]}"
        assert out[4] == score, f"Wrong score at idx {i}: {out[4]}"


@pytest.mark.parametrize(
    ("dbref_props", "expected_xrefs"),
    [
        (
            # all cross-ref fields present
            [
                ("UniRef90 ID", "UniRef90_N8Q6C0"),
                ("UniRef50 ID", "UniRef50_A0A7Z7LP76"),
                ("UniParc ID", "UPI00044F6C4F"),
                ("protein name", "foo"),
            ],
            [
                ("UniRef90 ID", "UniRef90_N8Q6C0"),
                ("UniRef50 ID", "UniRef50_A0A7Z7LP76"),
                ("UniParc ID", "UPI00044F6C4F"),
            ],
        ),
        (
            # partial cross-ref
            [
                ("UniRef90 ID", "UniRef90_ABC"),
                ("protein name", "bar"),
            ],
            [
                ("UniRef90 ID", "UniRef90_ABC"),
            ],
        ),
        (
            # No cross-ref
            [
                ("protein name", "baz"),
            ],
            [],
        ),
    ],
)
def test_extract_cross_refs_param(dbref_props, expected_xrefs) -> None:
    """
    Test that extract_cross_refs correctly extracts all UniRef cross-reference fields.
    """
    dbref = ET.Element(
        "{http://uniprot.org/uniref}dbReference",
        type="UniProtKB ID",
        id="TEST_ID",
    )

    for t, v in dbref_props:
        ET.SubElement(
            dbref,
            "{http://uniprot.org/uniref}property",
            type=t,
            value=v,
        )

    ns = {"ns": "http://uniprot.org/uniref"}
    cross_reference_data = []
    extract_cross_refs(dbref, cross_reference_data, ns)

    entity_id = cdm_entity_id("TEST_ID")
    expected = {(entity_id, typ, val) for typ, val in expected_xrefs}
    got = set(cross_reference_data)
    assert got == expected


@pytest.mark.parametrize(
    ("xml_str", "existing_created", "expected_created", "expect_member_count", "expect_xref_count"),
    [
        pytest.param(
            # CASE 1:
            # The old creation time exists, it should be retained.
            """
            <entry xmlns="http://uniprot.org/uniref"
                id="UniRef100_TEST"
                updated="2024-01-01">
                <name>Cluster: Example protein</name>

                <representativeMember>
                    <dbReference type="UniProtKB ID" id="REP_ID">
                        <property type="UniProtKB accession" value="REP_ACC"/>
                        <property type="isSeed" value="true"/>
                        <property type="UniParc ID" value="UPI000001"/>
                    </dbReference>
                </representativeMember>

                <member>
                    <dbReference type="UniProtKB ID" id="MEM_ID">
                        <property type="UniProtKB accession" value="MEM_ACC"/>
                    </dbReference>
                </member>
            </entry>
            """,
            {"UniRef100_TEST": "2020-01-01T00:00:00"},
            "2020-01-01T00:00:00",  # expected created
            2,  # 1 repr + 1 member
            1,  # one UniParc xref
            id="with_existing_created",
        ),
        pytest.param(
            # CASE 2:
            # No existing_created → created == updated
            """
            <entry xmlns="http://uniprot.org/uniref"
                id="UniRef100_NEW"
                updated="2024-03-01">
                <name>Cluster: New protein</name>

                <representativeMember>
                    <dbReference type="UniProtKB ID" id="REP2">
                        <property type="UniProtKB accession" value="REP2_ACC"/>
                        <property type="isSeed" value="true"/>
                    </dbReference>
                </representativeMember>
            </entry>
            """,
            {},
            None,  # meaning created == updated
            1,  # only representative
            0,  # no xrefs
            id="no_existing_created",
        ),
    ],
)
def test_parse_uniref_entry_param(xml_str, existing_created, expected_created, expect_member_count, expect_xref_count) -> None:
    elem = ET.fromstring(xml_str)

    result = parse_uniref_entry(elem, existing_created, NS)

    cluster_rows = result["cluster_data"]
    entity_rows = result["entity_data"]
    member_rows = result["cluster_member_data"]
    xref_rows = result["cross_reference_data"]

    # -----------------------------
    # Validate cluster
    # -----------------------------
    assert len(cluster_rows) == 1
    cluster_id = cluster_rows[0][0]
    assert cluster_id.startswith("cdm_ccol_")

    # -----------------------------
    # Validate Entity
    # -----------------------------
    assert len(entity_rows) == 1
    (
        ent_entity_id,
        _data_source_entity_id,
        ent_type,
        data_source,
        updated,
        created,
    ) = entity_rows[0]

    assert ent_entity_id == cluster_id
    assert data_source.startswith("UniRef")
    assert ent_type == "Cluster"

    if expected_created is not None:
        assert created == expected_created
        assert updated != created
    else:
        # created == updated when no existing_created
        assert created == updated

    # -----------------------------
    # Validate members
    # -----------------------------
    assert len(member_rows) == expect_member_count
    for row in member_rows:
        cid, entity_id, _, _, _ = row
        assert cid == cluster_id
        assert entity_id.startswith("cdm_prot_")

    # -----------------------------
    # Validate xrefs
    # -----------------------------
    assert len(xref_rows) == expect_xref_count
    for e_id, _x_type, _x_val in xref_rows:
        # xrefs use CDM: prefix
        assert e_id.startswith("CDM:")


def make_fake_uniref_xml(num_entries=2):
    """
    Create a minimal UniRef XML with N <entry> elements.
    """
    entries = []
    for i in range(num_entries):
        entries.append(f"""
        <entry xmlns="http://uniprot.org/uniref" id="UniRef100_FAKE{i}">
            <name>Cluster FAKE {i}</name>
            <representativeMember>
                <dbReference type="UniProtKB ID" id="FAKE{i}_REP">
                    <property type="UniProtKB accession" value="FAKE{i}_REP_ACC"/>
                    <property type="isSeed" value="true"/>
                </dbReference>
            </representativeMember>
        </entry>
        """)

    xml = f"""
    <root xmlns="http://uniprot.org/uniref">
        {"".join(entries)}
    </root>
    """
    return xml.strip().encode("utf-8")


@pytest.mark.parametrize("batch_size", [1, 2])
def test_parse_uniref_xml_batch(batch_size) -> None:
    # Prepare fake XML inside a gzipped temp file
    with tempfile.NamedTemporaryFile(suffix=".xml.gz", delete=True) as tmp:
        xml_bytes = make_fake_uniref_xml(num_entries=2)

        with gzip.open(tmp.name, "wb") as gz:
            gz.write(xml_bytes)

        # No existing created timestamps
        existing_created = {}

        # Call the parser
        result = parse_uniref_xml(tmp.name, batch_size, existing_created)

        # Validate the number of parsed entries matches batch_size
        assert len(result["cluster_data"]) == batch_size
        assert len(result["entity_data"]) == batch_size

        # Member data: each entry has exactly 1 representative member → 1 row per entry
        assert len(result["cluster_member_data"]) == batch_size

        # Cross references: none included in fake XML
        assert len(result["cross_reference_data"]) == 0

        # Validate cluster_id prefix
        cluster_id = result["cluster_data"][0][0]
        assert cluster_id.startswith("cdm_ccol_")
