import pytest
import textwrap
from datetime import datetime
import xml.etree.ElementTree as ET
from uniref_revise import extract_cluster,cdm_entity_id,get_timestamps,get_accession_and_seed,add_cluster_members,extract_cross_refs


@pytest.mark.parametrize("accession,expected_prefix", [
    ("A0B0123456", "CDM:"),
    ("P01234", "CDM:"),
    ("", None),
    (None, None)])
def test_cdm_entity_id(accession, expected_prefix):
    result = cdm_entity_id(accession)
    if expected_prefix is None:
        assert result is None
    else:
        assert result.startswith(expected_prefix)


@pytest.mark.parametrize("xml_str, expected_name", [
        (
        "<entry xmlns='http://uniprot.org/uniref' id='UniRef100_A0A009GP46' updated='2016-10-05'>"
        "<name>TestName</name></entry>", 
        "TestName"
    ),
    (
        "<entry xmlns='http://uniprot.org/uniref' id='UniRef100_XYZ' updated='2024-01-01'/>",
        "UNKNOWN"
    )
])

def test_extract_cluster(xml_str, expected_name):
    ns = {"ns": "http://uniprot.org/uniref"}
    elem = ET.fromstring(xml_str)
    cluster_id, name = extract_cluster(elem, ns)
    assert cluster_id.startswith("CDM:")
    assert isinstance(cluster_id, str)
    assert name == expected_name


@pytest.mark.parametrize(
    "uniref_id, existing_created, now, expected",
    [
        # Has existing_created 
        (
            "UniRef100_A",
            {"UniRef100_A": "2024-01-01T00:00:00"},
            datetime(2025, 1, 1, 0, 0, 0),
            ("2025-01-01T00:00:00", "2024-01-01T00:00:00")
        ),
        # There is no existing_created
        (
            "UniRef100_B",
            {"UniRef100_A": "2024-01-01T00:00:00"},
            datetime(2025, 1, 1, 0, 0, 0),
            ("2025-01-01T00:00:00", "2025-01-01T00:00:00")
        ),
        # There is no existing_created，also not provide "now"
        (
            "UniRef100_C",
            {},
            None,  # The system automatically use the current time
            None   # Only assert that the return is a string and they are equal
        ),
    ]
)

def test_get_timestamps(uniref_id, existing_created, now, expected):
    result = get_timestamps(uniref_id, existing_created, now)
    if expected is not None:
        assert result == expected
    else:
        formatted_now, created_time = result
        assert formatted_now == created_time
        assert isinstance(formatted_now, str)
        assert len(formatted_now) == 19  # "YYYY-MM-DDTHH:MM:SS" ---> 19 bites 


@pytest.mark.parametrize("xml_str, expected_acc, expected_is_seed", [
    # Have accession and isSeed
    (
        '''
        <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="A0A009HJL9_ACIB9">
            <property type="UniProtKB accession" value="A0A009HJL9"/>
            <property type="isSeed" value="true"/>
        </dbReference>
        ''',
        "A0A009HJL9", True
    ),
    # Only accession，No isSeed
    (
        '''
        <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="A0A241V597_9GAMM">
            <property type="UniProtKB accession" value="A0A241V597"/>
        </dbReference>
        ''',
        "A0A241V597", False
    ),
    # No accession，only id
    (
        '''
        <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="ID_ONLY"></dbReference>
        ''',
        "ID_ONLY", False
    ),
    # None
    (
        None, None, False
    ),
])

def test_get_accession_and_seed(xml_str, expected_acc, expected_is_seed):
    ns = {"ns": "http://uniprot.org/uniref"}
    dbref = ET.fromstring(xml_str) if xml_str else None
    acc, is_seed = get_accession_and_seed(dbref, ns)
    assert acc == expected_acc
    assert is_seed == expected_is_seed



def make_entry_with_members(member_xmls, ns_uri="http://uniprot.org/uniref"):
    """
    receives a list of xml strings from dbReference, generates an <entry> element with <member> child nodes

    """
    entry_elem = ET.Element(f"{{{ns_uri}}}entry")
    for dbref_xml in member_xmls:
        dbref_elem = ET.fromstring(dbref_xml)
        member_elem = ET.SubElement(entry_elem, f"{{{ns_uri}}}member")
        member_elem.append(dbref_elem)
    return entry_elem

@pytest.mark.parametrize(
    "repr_xml, member_xmls, expected",
    [
        pytest.param(
            # representative member, two members
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
                """)
            ],
            [
                ("CLUSTER_X", "CDM:", "true", "true", "1.0"),
                ("CLUSTER_X", "CDM:", "false", "false", "1.0"),
                ("CLUSTER_X", "CDM:", "false", "true", "1.0"),
            ],
            id="with-representative-and-members"
        ),
        pytest.param(
            # Only memebers, no representative member
            None,
            [
                textwrap.dedent("""
                    <dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="MEM_ID">
                        <property type="UniProtKB accession" value="MEM_ACC"/>
                    </dbReference>
                """)
            ],
            [
                ("CLUSTER_X", "CDM:", "false", "false", "1.0")
            ],
            id="members-only"
        ),
        pytest.param(
            # No members, no representative member
            None,
            [],
            [],
            id="no-members"
        ),
    ]
)
def test_add_cluster_members(repr_xml, member_xmls, expected):
    """ Test add_cluster_members with various representative/member combinations """

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
   
 
XREF_TYPES = ["UniRef90 ID", "UniRef50 ID", "UniParc ID"]
@pytest.mark.parametrize(
    "dbref_props, expected_xrefs",
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
            ]
        ),
        (
            # partial cross-ref
            [
                ("UniRef90 ID", "UniRef90_ABC"),
                ("protein name", "bar"),
            ],
            [
                ("UniRef90 ID", "UniRef90_ABC"),
            ]
        ),
        (
            # No cross-ref
            [
                ("protein name", "baz"),
            ],
            []
        ),
    ]
)
def test_extract_cross_refs_param(dbref_props, expected_xrefs):
    """
    Test that extract_cross_refs correctly extracts all UniRef cross-reference fields
    """
    dbref = ET.Element("{http://uniprot.org/uniref}dbReference", type="UniProtKB ID", id="TEST_ID")

    for t, v in dbref_props:
        ET.SubElement(dbref, "{http://uniprot.org/uniref}property", type=t, value=v)

    ns = {"ns": "http://uniprot.org/uniref"}
    cross_reference_data = []
    extract_cross_refs(dbref, cross_reference_data, ns)

    entity_id = cdm_entity_id("TEST_ID")
    expected = set((entity_id, typ, val) for typ, val in expected_xrefs)
    got = set(cross_reference_data)
    assert got == expected

    