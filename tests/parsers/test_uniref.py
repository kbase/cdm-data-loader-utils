import os
import sys

sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

import gzip
import tempfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
import pytest

from cdm_data_loader_utils.parsers.uniref import (
    cdm_entity_id,
    get_timestamps,
    extract_cluster,
    get_accession_and_seed,
    add_cluster_members,
    extract_cross_refs,
    parse_uniref_xml,
)

NS = {"ns": "http://uniprot.org/uniref"}


# ---------------------------------------------------------
# cdm_entity_id
# ---------------------------------------------------------
@pytest.mark.parametrize(
    "value, should_raise",
    [
        ("A0A009HJL9", False),
        ("UniRef100_A0A009HJL9", False),
        ("", True),
        (None, True),
    ],
)
def test_cdm_entity_id(value, should_raise):
    if should_raise:
        with pytest.raises(ValueError):
            cdm_entity_id(value)
    else:
        out = cdm_entity_id(value)
        assert isinstance(out, str)
        assert out.startswith("CDM:")


# ---------------------------------------------------------
# get_timestamps
# ---------------------------------------------------------
@pytest.mark.parametrize(
    "uniref_id, existing, now, expect_created_same_as_updated",
    [
        (
            "UniRef100_A",
            {"UniRef100_A": "2024-01-01T00:00:00+00:00"},
            datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc),
            False,
        ),
        (
            "UniRef100_B",
            {},
            datetime(2025, 1, 1, 0, 0, 0, tzinfo=timezone.utc),
            True,
        ),
        (
            "UniRef100_C",
            {},
            None,
            True,
        ),
    ],
)
def test_get_timestamps(uniref_id, existing, now, expect_created_same_as_updated):
    updated, created = get_timestamps(uniref_id, existing, now)

    assert isinstance(updated, str)
    assert isinstance(created, str)
    assert updated.endswith("+00:00")

    if expect_created_same_as_updated:
        assert updated == created
    else:
        assert updated != created


@pytest.mark.parametrize("bad_id", ["", None])
def test_get_timestamps_rejects_empty_uniref_id(bad_id):
    with pytest.raises(ValueError):
        get_timestamps(bad_id, {}, None)


# ---------------------------------------------------------
# add_cluster_members
# ---------------------------------------------------------
@pytest.mark.parametrize(
    "repr_xml, member_xmls, expected_count",
    [
        (
            """
            <dbReference xmlns="http://uniprot.org/uniref">
                <property type="UniProtKB accession" value="REP_ACC"/>
                <property type="isSeed" value="true"/>
            </dbReference>
            """,
            [
                """
                <dbReference xmlns="http://uniprot.org/uniref">
                    <property type="UniProtKB accession" value="MEM1_ACC"/>
                </dbReference>
                """,
                """
                <dbReference xmlns="http://uniprot.org/uniref">
                    <property type="UniProtKB accession" value="MEM2_ACC"/>
                </dbReference>
                """,
            ],
            3,
        ),
        (
            None,
            [
                """
                <dbReference xmlns="http://uniprot.org/uniref">
                    <property type="UniProtKB accession" value="ONLY_ACC"/>
                </dbReference>
                """,
            ],
            1,
        ),
        (None, [], 0),
    ],
)
def test_add_cluster_members(repr_xml, member_xmls, expected_count):
    cluster_id = "CDM_CLUSTER"
    repr_db = ET.fromstring(repr_xml) if repr_xml else None

    entry = ET.Element("{http://uniprot.org/uniref}entry")
    for m in member_xmls:
        mem = ET.SubElement(entry, "{http://uniprot.org/uniref}member")
        mem.append(ET.fromstring(m))

    rows = []
    add_cluster_members(cluster_id, repr_db, entry, rows, NS)

    assert len(rows) == expected_count
    for r in rows:
        assert r[0] == cluster_id
        assert r[1].startswith("CDM:")
        assert r[4] == "1.0"


# ---------------------------------------------------------
# extract_cluster
# ---------------------------------------------------------
@pytest.mark.parametrize(
    "xml_str, uniref_id, expected_name",
    [
        (
            "<entry xmlns='http://uniprot.org/uniref' id='UniRef100_A'><name>Test Cluster Name</name></entry>",
            "UniRef100_A",
            "Test Cluster Name",
        ),
        (
            "<entry xmlns='http://uniprot.org/uniref' id='UniRef100_B'/>",
            "UniRef100_B",
            "UNKNOWN",
        ),
    ],
)
def test_extract_cluster(xml_str, uniref_id, expected_name):
    elem = ET.fromstring(xml_str)

    cluster_id, name = extract_cluster(elem, NS, uniref_id)

    # ---- cluster_id checks ----
    assert isinstance(cluster_id, str)
    assert cluster_id.startswith("CDM:")

    # ---- name checks ----
    assert name == expected_name


@pytest.mark.parametrize(
    "xml_str, expected_acc, expected_is_seed",
    [
        # accession + isSeed=true
        (
            """
            <dbReference xmlns="http://uniprot.org/uniref">
                <property type="UniProtKB accession" value="A0A009HJL9"/>
                <property type="isSeed" value="true"/>
            </dbReference>
            """,
            "A0A009HJL9",
            True,
        ),
        # accession only
        (
            """
            <dbReference xmlns="http://uniprot.org/uniref">
                <property type="UniProtKB accession" value="A0A241V597"/>
            </dbReference>
            """,
            "A0A241V597",
            False,
        ),
        # no accession
        (
            """
            <dbReference xmlns="http://uniprot.org/uniref">
                <property type="protein name" value="Some protein"/>
            </dbReference>
            """,
            None,
            False,
        ),
        # dbref is None
        (
            None,
            None,
            False,
        ),
    ],
)
def test_get_accession_and_seed(xml_str, expected_acc, expected_is_seed):
    dbref = ET.fromstring(xml_str) if xml_str else None

    acc, is_seed = get_accession_and_seed(dbref, NS)

    assert acc == expected_acc
    assert is_seed == expected_is_seed


# ---------------------------------------------------------
# extract_cross_refs
# ---------------------------------------------------------
@pytest.mark.parametrize(
    "props, expected",
    [
        (
            [
                ("UniProtKB accession", "A0A1"),
                ("UniRef90 ID", "UniRef90_X"),
                ("UniParc ID", "UPI0001"),
            ],
            {
                ("UniRef90 ID", "UniRef90_X"),
                ("UniParc ID", "UPI0001"),
            },
        ),
        (
            [
                ("UniProtKB accession", "A0A2"),
            ],
            set(),
        ),
    ],
)
def test_extract_cross_refs(props, expected):
    dbref = ET.Element("{http://uniprot.org/uniref}dbReference", id="UniProtKB:A0A1")

    for k, v in props:
        ET.SubElement(
            dbref,
            "{http://uniprot.org/uniref}property",
            type=k,
            value=v,
        )

    rows = []
    extract_cross_refs(dbref, rows, NS)

    got = {(t, v) for _, t, v in rows}
    assert got == expected

    for entity_id, _, _ in rows:
        assert entity_id is not None
        assert isinstance(entity_id, str)


# ---------------------------------------------------------
# parse_uniref_xml
# ---------------------------------------------------------
@pytest.mark.parametrize("batch_size", [1, 2])
def test_parse_uniref_xml_batch(batch_size):
    xml = """
    <root xmlns="http://uniprot.org/uniref">
        <entry id="UniRef100_A">
            <name>A</name>
            <representativeMember>
                <dbReference id="UniProtKB:A1">
                    <property type="UniProtKB accession" value="A1"/>
                </dbReference>
            </representativeMember>
        </entry>

        <entry id="UniRef100_B">
            <name>B</name>
            <representativeMember>
                <dbReference id="UniProtKB:B1">
                    <property type="UniProtKB accession" value="B1"/>
                </dbReference>
            </representativeMember>
        </entry>
    </root>
    """.strip()

    with tempfile.TemporaryDirectory() as tmpdir:
        gz_path = f"{tmpdir}/uniref_test.xml.gz"
        with gzip.open(gz_path, "wb") as gz:
            gz.write(xml.encode("utf-8"))

        result = parse_uniref_xml(gz_path, batch_size, {})

    assert len(result["cluster_data"]) == batch_size
    assert len(result["entity_data"]) == batch_size
    assert len(result["cluster_member_data"]) == batch_size
    assert len(result["cross_reference_data"]) in (0, batch_size)
