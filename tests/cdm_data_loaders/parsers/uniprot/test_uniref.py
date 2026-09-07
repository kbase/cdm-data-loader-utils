"""Tests for the UniRef Parser."""

import datetime
import logging
from pathlib import Path
from typing import Any

import pytest
from lxml.etree import XMLParser, fromstring

from cdm_data_loaders.parsers.uniprot.uniref import (
    ENTRY_XML_TAG,
    UNIREF_URL,
    dump_xml_element,
    extract_cluster,
    extract_cross_refs,
    generate_dbxref,
    parse_uniref_entry,
)
from cdm_data_loaders.readers.xml import stream_xml_file

NOW = datetime.datetime.now(tz=datetime.UTC)
PROTOCOL = "uniref_9000"
FILE_PATH = "/path/to/file"


def entry_wrap(xml: str) -> str:
    """Wrap XML text in an `entry` element.

    :param xml: the XML to wrap
    :type xml: str
    :return: XML string with an entry around it.
    :rtype: str
    """
    return f"<entry xmlns='{UNIREF_URL}'>{xml}</entry>"


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        (
            """  <entry xmlns='http://uniprot.org/uniref' id="UniRef50_A0AB34IYJ6" updated="2025-02-05">
    <name>Cluster: PKZILLA-1 giant polyketide synthase</name>
    <property type="member count" value="1"/>
    <property type="common taxon" value="Prymnesium parvum"/>
    <property type="common taxon ID" value="97485"/>
    <property type="GO Molecular Function" value="GO:0004315"/>
    <property type="GO Molecular Function" value="GO:0004312"/>
    <property type="GO Molecular Function" value="GO:0016491"/>
    <property type="GO Molecular Function" value="GO:0031177"/>
    <property type="GO Biological Process" value="GO:0006633"/>
    <property type="GO Biological Process" value="GO:0044550"/>
    <representativeMember>
      <dbReference type="UniProtKB ID" id="A0AB34IYJ6_PRYPA">
        <property type="UniProtKB accession" value="A0AB34IYJ6"/>
        <property type="UniParc ID" value="UPI003596FD28"/>
        <property type="UniRef100 ID" value="UniRef100_A0AB34IYJ6"/>
        <property type="UniRef90 ID" value="UniRef90_A0AB34IYJ6"/>
        <property type="protein name" value="PKZILLA-1 giant polyketide synthase"/>
        <property type="source organism" value="Prymnesium parvum (Golden alga)"/>
        <property type="NCBI taxonomy" value="97485"/>
        <property type="length" value="45212"/>
        <property type="isSeed" value="true"/>
      </dbReference>
      <sequence length="45212" checksum="385E596EB33D030E">SNIP</sequence>
    </representativeMember>
  </entry>
""",
            (
                {
                    "cluster_id": "uniref:UniRef50_A0AB34IYJ6",
                    "name": "Cluster: PKZILLA-1 giant polyketide synthase",
                    "cluster_type": "Protein",
                    "description": None,
                },
                {
                    "entity_id": "uniref:UniRef50_A0AB34IYJ6",
                    "entity_type": "Cluster",
                    "data_source": "UniRef",
                    "data_source_entity_id": "UniRef50_A0AB34IYJ6",
                    "data_source_updated": "2025-02-05",
                },
            ),
        ),
        (
            """
    <entry xmlns='http://uniprot.org/uniref' id="UniRef50_UPI00358F51CD" updated="2024-11-27">
        <name>Cluster: LOW QUALITY PROTEIN: titin</name>
        <property type="member count" value="1"/>
        <property type="common taxon" value="Myxine glutinosa"/>
        <property type="common taxon ID" value="7769"/>
        <representativeMember>
        <dbReference type="UniParc ID" id="UPI00358F51CD">
            <property type="UniRef100 ID" value="UniRef100_UPI00358F51CD"/>
            <property type="UniRef90 ID" value="UniRef90_UPI00358F51CD"/>
            <property type="protein name" value="LOW QUALITY PROTEIN: titin"/>
            <property type="source organism" value="Myxine glutinosa"/>
            <property type="NCBI taxonomy" value="7769"/>
            <property type="length" value="47063"/>
            <property type="isSeed" value="true"/>
        </dbReference>
        <sequence length="47063" checksum="48729625616C010E">SNIP</sequence>
        </representativeMember>
    </entry>
    """,
            (
                {
                    "cluster_id": "uniref:UniRef50_UPI00358F51CD",
                    "name": "Cluster: LOW QUALITY PROTEIN: titin",
                    "cluster_type": "Protein",
                    "description": None,
                },
                {
                    "entity_id": "uniref:UniRef50_UPI00358F51CD",
                    "entity_type": "Cluster",
                    "data_source": "UniRef",
                    "data_source_entity_id": "UniRef50_UPI00358F51CD",
                    "data_source_updated": "2024-11-27",
                },
            ),
        ),
        (
            """
    <entry xmlns='http://uniprot.org/uniref' id="UniRef50_UPI0035" updated="2022-11-27">
        <property type="member count" value="1"/>
        <property type="common taxon" value="Myxine glutinosa"/>
        <property type="common taxon ID" value="7769"/>
    </entry>
    """,
            (
                {
                    "cluster_id": "uniref:UniRef50_UPI0035",
                    "name": None,
                    "cluster_type": "Protein",
                    "description": None,
                },
                {
                    "entity_id": "uniref:UniRef50_UPI0035",
                    "entity_type": "Cluster",
                    "data_source": "UniRef",
                    "data_source_entity_id": "UniRef50_UPI0035",
                    "data_source_updated": "2022-11-27",
                },
            ),
        ),
    ],
)
def test_extract_cluster(xml: str, expected: tuple[dict[str, Any], dict[str, Any]]) -> None:
    """Test the extraction of cluster information from an entry."""
    elem = fromstring(xml)
    (cluster, entity) = extract_cluster(elem)
    assert (cluster, entity) == expected


ENTITY_ID = "uniprot:A2ASS6"
UNIPROT_NAME_ID = "uniprot_name:TITIN_MOUSE"
CLUSTER_DATA = {
    "entity_id": ENTITY_ID,
    "cluster_id": "cluster:id",
    "is_seed": False,
}
CLUSTER_SEED_DATA = {**CLUSTER_DATA, "is_seed": True}


@pytest.mark.parametrize(
    ("xml", "expected", "warning"),
    [
        (
            # all cross-ref fields present
            """<dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="TITIN_MOUSE">
        <property type="UniProtKB accession" value="A2ASS6"/>
        <property type="UniProtKB accession" value="A2ASS5"/>
        <property type="UniProtKB accession" value="A2ASS7"/>
        <property type="UniParc ID" value="UPI0000D77B45"/>
        <property type="UniRef100 ID" value="UniRef100_A2ASS6"/>
        <property type="UniRef90 ID" value="UniRef90_A2ASS6"/>
        <property type="protein name" value="Titin"/>
        <property type="source organism" value="Mus musculus (Mouse)"/>
        <property type="NCBI taxonomy" value="10090"/>
        <property type="length" value="35213"/>
      </dbReference>
            """,
            CLUSTER_DATA,
            None,
        ),
        (
            # xrefs only
            """<dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="TITIN_MOUSE">
        <property type="UniProtKB accession" value="A2ASS6"/>
      </dbReference>
            """,
            CLUSTER_DATA,
            None,
        ),
        (
            # protein and name info, isSeed is true
            """<dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="TITIN_MOUSE">
        <property type="UniProtKB accession" value="A2ASS6"/>
        <property type="protein name" value="Titin"/>
        <property type="length" value="35213"/>
        <property type="isSeed" value="true"/>
      </dbReference>
            """,
            CLUSTER_SEED_DATA,
            None,
        ),
        (
            # no uniprot acc for a uniprot ID
            """<dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID" id="TITIN_MOUSE">
        <property type="source organism" value="Mus musculus (Mouse)"/>
        <property type="fave breakfast cereal" value="Cheerios"/>
        <property type="height in inches" value="3"/>
      </dbReference>
            """,
            {**CLUSTER_DATA, "entity_id": UNIPROT_NAME_ID},
            "Could not find UniProt accession for UniProtKB ID TITIN_MOUSE",
        ),
        (  # incomplete
            """<dbReference xmlns="http://uniprot.org/uniref" type="UniProtKB ID"></dbReference>""",
            {},
            "Incomplete dbReference in cluster:id",
        ),
        (  # incomplete
            """<dbReference xmlns="http://uniprot.org/uniref" id="TITIN_MOUSE" />""",
            {},
            "Incomplete dbReference in cluster:id",
        ),
        (  # incomplete
            """<dbReference xmlns="http://uniprot.org/uniref" />""",
            {},
            "Incomplete dbReference in cluster:id",
        ),
    ],
)
@pytest.mark.parametrize("is_representative", [True, False])
def test_extract_cross_refs_param(
    xml: str, expected: dict[str, Any], warning: str | None, is_representative: bool, caplog: pytest.LogCaptureFixture
) -> None:
    """
    Test that extract_cross_refs correctly extracts all UniRef cross-reference fields.
    """
    # edit the cluster member as appropriate
    if expected:
        expected["is_representative"] = is_representative

    db_ref = extract_cross_refs(fromstring(xml), "cluster:id", is_representative)
    assert db_ref == expected

    if warning:
        assert len(caplog.records) == 1
        assert caplog.records[0].levelno == logging.WARNING
        assert caplog.records[0].message == warning
    else:
        assert not caplog.records


@pytest.mark.parametrize(
    "entry",
    [
        pytest.param(
            {
                "xml": """  <entry xmlns="http://uniprot.org/uniref" id="UniRef50_A0AB34IYJ6" updated="2025-02-05">
    <name>Cluster: PKZILLA-1 giant polyketide synthase</name>
    <property type="member count" value="1"/>
    <property type="common taxon" value="Prymnesium parvum"/>
    <property type="common taxon ID" value="97485"/>
    <property type="GO Molecular Function" value="GO:0004315"/>
    <property type="GO Molecular Function" value="GO:0004312"/>
    <property type="GO Molecular Function" value="GO:0016491"/>
    <property type="GO Molecular Function" value="GO:0031177"/>
    <property type="GO Biological Process" value="GO:0006633"/>
    <property type="GO Biological Process" value="GO:0044550"/>
    <representativeMember>
      <dbReference type="UniProtKB ID" id="A0AB34IYJ6_PRYPA">
        <property type="UniProtKB accession" value="A0AB34IYJ6"/>
        <property type="UniParc ID" value="UPI003596FD28"/>
        <property type="UniRef100 ID" value="UniRef100_A0AB34IYJ6"/>
        <property type="UniRef90 ID" value="UniRef90_A0AB34IYJ6"/>
        <property type="protein name" value="PKZILLA-1 giant polyketide synthase"/>
        <property type="source organism" value="Prymnesium parvum (Golden alga)"/>
        <property type="NCBI taxonomy" value="97485"/>
        <property type="length" value="45212"/>
        <property type="isSeed" value="true"/>
      </dbReference>
      <sequence length="45212" checksum="385E596EB33D030E">SNIP</sequence>
    </representativeMember>
  </entry>
""",
                "exp": {
                    "cluster": [
                        {
                            "cluster_id": "uniref:UniRef50_A0AB34IYJ6",
                            "name": "Cluster: PKZILLA-1 giant polyketide synthase",
                            "cluster_type": "Protein",
                            "protocol": PROTOCOL,
                            "description": None,
                        }
                    ],
                    "entity": [
                        {
                            "entity_id": "uniref:UniRef50_A0AB34IYJ6",
                            "entity_type": "Cluster",
                            "data_source": "UniRef",
                            "data_source_entity_id": "UniRef50_A0AB34IYJ6",
                            "data_source_updated": "2025-02-05",
                            "updated": NOW,
                        }
                    ],
                    "clustermember": [
                        {
                            "entity_id": "uniprot:A0AB34IYJ6",
                            "cluster_id": "uniref:UniRef50_A0AB34IYJ6",
                            "is_representative": True,
                            "is_seed": True,
                        },
                    ],
                },
            },
            id="Single UniProt rep",
        ),
        pytest.param(
            {
                "xml": """  <entry xmlns="http://uniprot.org/uniref" id="UniRef50_A0A410P257" updated="2024-07-24">
    <name>Cluster: Glycogen synthase</name>
    <property type="member count" value="2"/>
    <property type="common taxon" value="Velamenicoccus archaeovorus"/>
    <property type="common taxon ID" value="1930593"/>
    <representativeMember>
      <dbReference type="UniProtKB ID" id="A0A410P257_VELA1">
        <property type="UniProtKB accession" value="A0A410P257"/>
        <property type="UniParc ID" value="UPI000FFCD83E"/>
        <property type="UniRef100 ID" value="UniRef100_A0A410P257"/>
        <property type="UniRef90 ID" value="UniRef90_A0A410P257"/>
        <property type="protein name" value="Glycogen synthase"/>
        <property type="source organism" value="Velamenicoccus archaeovorus"/>
        <property type="NCBI taxonomy" value="1930593"/>
        <property type="length" value="39677"/>
      </dbReference>
      <sequence length="39677" checksum="976AC558673E4C32">SNIP</sequence>
    </representativeMember>
    <member>
      <dbReference type="UniParc ID" id="UPI000FFEDAD7">
        <property type="UniRef100 ID" value="UniRef100_UPI000FFEDAD7"/>
        <property type="UniRef90 ID" value="UniRef90_A0A410P257"/>
        <property type="protein name" value="UPF0489 family protein"/>
        <property type="source organism" value="Candidatus Velamenicoccus archaeovorus"/>
        <property type="NCBI taxonomy" value="1930593"/>
        <property type="length" value="39686"/>
        <property type="isSeed" value="true"/>
      </dbReference>
    </member>
  </entry>
""",
                "exp": {
                    "cluster": [
                        {
                            "cluster_id": "uniref:UniRef50_A0A410P257",
                            "name": "Cluster: Glycogen synthase",
                            "cluster_type": "Protein",
                            "protocol": PROTOCOL,
                            "description": None,
                        }
                    ],
                    "entity": [
                        {
                            "entity_id": "uniref:UniRef50_A0A410P257",
                            "entity_type": "Cluster",
                            "data_source": "UniRef",
                            "data_source_entity_id": "UniRef50_A0A410P257",
                            "data_source_updated": "2024-07-24",
                            "updated": NOW,
                        }
                    ],
                    "clustermember": [
                        {
                            "entity_id": "uniprot:A0A410P257",
                            "cluster_id": "uniref:UniRef50_A0A410P257",
                            "is_representative": True,
                            "is_seed": False,
                        },
                        {
                            "entity_id": "uniparc:UPI000FFEDAD7",
                            "cluster_id": "uniref:UniRef50_A0A410P257",
                            "is_representative": False,
                            "is_seed": True,
                        },
                    ],
                },
            },
            id="Rep different to seed",
        ),
        pytest.param(
            {
                "xml": """  <entry xmlns="http://uniprot.org/uniref" id="UniRef50_Q8WZ42" updated="2025-10-08">
    <name>Cluster: Titin</name>
    <property type="member count" value="4624"/>
    <property type="common taxon" value="Vertebrata"/>
    <property type="common taxon ID" value="7742"/>
    <property type="GO Biological Process" value="GO:0009987"/>
    <property type="GO Cellular Component" value="GO:0110165"/>
    <representativeMember>
      <dbReference type="UniProtKB ID" id="TITIN_HUMAN">
        <property type="UniProtKB accession" value="Q8WZ42"/>
        <property type="UniProtKB accession" value="A0A0A0MTS7"/>
        <property type="UniProtKB accession" value="A6NKB1"/>
        <property type="UniProtKB accession" value="C9JQJ2"/>
        <property type="UniProtKB accession" value="E7EQE6"/>
        <property type="UniProtKB accession" value="E7ET18"/>
        <property type="UniProtKB accession" value="K7ENY1"/>
        <property type="UniProtKB accession" value="Q10465"/>
        <property type="UniProtKB accession" value="Q10466"/>
        <property type="UniProtKB accession" value="Q15598"/>
        <property type="UniProtKB accession" value="Q2XUS3"/>
        <property type="UniProtKB accession" value="Q32Q60"/>
        <property type="UniProtKB accession" value="Q4U1Z6"/>
        <property type="UniProtKB accession" value="Q4ZG20"/>
        <property type="UniProtKB accession" value="Q6NSG0"/>
        <property type="UniParc ID" value="UPI00025287CD"/>
        <property type="UniRef100 ID" value="UniRef100_Q8WZ42"/>
        <property type="UniRef90 ID" value="UniRef90_Q8WZ42"/>
        <property type="protein name" value="Titin"/>
        <property type="source organism" value="Homo sapiens (Human)"/>
        <property type="NCBI taxonomy" value="9606"/>
        <property type="length" value="34350"/>
      </dbReference>
      <sequence length="34350" checksum="DEB216410AD560D9">SNIP</sequence>
    </representativeMember>
    <member>
      <dbReference type="UniProtKB ID" id="TITIN_MOUSE">
        <property type="UniProtKB accession" value="A2ASS6"/>
        <property type="UniProtKB accession" value="A2ASS5"/>
        <property type="UniProtKB accession" value="A2ASS7"/>
        <property type="UniParc ID" value="UPI0000D77B45"/>
        <property type="UniRef100 ID" value="UniRef100_A2ASS6"/>
        <property type="UniRef90 ID" value="UniRef90_A2ASS6"/>
        <property type="protein name" value="Titin"/>
        <property type="source organism" value="Mus musculus (Mouse)"/>
        <property type="NCBI taxonomy" value="10090"/>
        <property type="length" value="35213"/>
      </dbReference>
    </member>
    <member>
      <dbReference type="UniProtKB ID" id="Q8WZ42-8">
        <property type="UniProtKB accession" value="Q8WZ42-8"/>
        <property type="UniParc ID" value="UPI000255FF54"/>
        <property type="UniRef100 ID" value="UniRef100_Q8WZ42-8"/>
        <property type="UniRef90 ID" value="UniRef90_Q8WZ42"/>
        <property type="protein name" value="Isoform 8 of Titin"/>
        <property type="source organism" value="Homo sapiens (Human)"/>
        <property type="NCBI taxonomy" value="9606"/>
        <property type="length" value="34475"/>
      </dbReference>
    </member>
    <member>
      <dbReference type="UniProtKB ID" id="A0A2J8PRH0_PANTR">
        <property type="UniProtKB accession" value="A0A2J8PRH0"/>
        <property type="UniParc ID" value="UPI000CB538BE"/>
        <property type="UniRef100 ID" value="UniRef100_A0A2J8PRH0"/>
        <property type="UniRef90 ID" value="UniRef90_Q8WZ42"/>
        <property type="protein name" value="Titin"/>
        <property type="source organism" value="Pan troglodytes (Chimpanzee)"/>
        <property type="NCBI taxonomy" value="9598"/>
        <property type="length" value="34356"/>
      </dbReference>
    </member>
    <member>
      <dbReference type="UniProtKB ID" id="SOME_FAKE_ID_WITH_NO_ACC" />
    </member>
    <!-- many more omitted -->
  </entry>
""",
                "exp": {
                    "cluster": [
                        {
                            "cluster_id": "uniref:UniRef50_Q8WZ42",
                            "name": "Cluster: Titin",
                            "cluster_type": "Protein",
                            "protocol": PROTOCOL,
                            "description": None,
                        }
                    ],
                    "entity": [
                        {
                            "entity_id": "uniref:UniRef50_Q8WZ42",
                            "entity_type": "Cluster",
                            "data_source": "UniRef",
                            "data_source_entity_id": "UniRef50_Q8WZ42",
                            "data_source_updated": "2025-10-08",
                            "updated": NOW,
                        }
                    ],
                    "clustermember": [
                        {
                            "entity_id": "uniprot:Q8WZ42",
                            "cluster_id": "uniref:UniRef50_Q8WZ42",
                            "is_representative": True,
                            "is_seed": False,
                        },
                        {
                            "entity_id": "uniprot:A2ASS6",
                            "cluster_id": "uniref:UniRef50_Q8WZ42",
                            "is_representative": False,
                            "is_seed": False,
                        },
                        {
                            "entity_id": "uniprot:Q8WZ42-8",
                            "cluster_id": "uniref:UniRef50_Q8WZ42",
                            "is_representative": False,
                            "is_seed": False,
                        },
                        {
                            "entity_id": "uniprot:A0A2J8PRH0",
                            "cluster_id": "uniref:UniRef50_Q8WZ42",
                            "is_representative": False,
                            "is_seed": False,
                        },
                        {
                            "entity_id": "uniprot_name:SOME_FAKE_ID_WITH_NO_ACC",
                            "cluster_id": "uniref:UniRef50_Q8WZ42",
                            "is_representative": False,
                            "is_seed": False,
                        },
                    ],
                },
            },
            id="multiple UniProt accs, one UniProt name",
        ),
        pytest.param(
            {
                "xml": """
  <entry xmlns="http://uniprot.org/uniref" id="UniRef50_UPI00358F51CD" updated="2024-11-27">
    <name>Cluster: LOW QUALITY PROTEIN: titin</name>
    <property type="member count" value="1"/>
    <property type="common taxon" value="Myxine glutinosa"/>
    <property type="common taxon ID" value="7769"/>
    <representativeMember>
      <dbReference type="UniParc ID" id="UPI00358F51CD">
        <property type="UniRef100 ID" value="UniRef100_UPI00358F51CD"/>
        <property type="UniRef90 ID" value="UniRef90_UPI00358F51CD"/>
        <property type="protein name" value="LOW QUALITY PROTEIN: titin"/>
        <property type="source organism" value="Myxine glutinosa"/>
        <property type="NCBI taxonomy" value="7769"/>
        <property type="length" value="47063"/>
        <property type="isSeed" value="true"/>
      </dbReference>
      <sequence length="47063" checksum="48729625616C010E">SNIP</sequence>
    </representativeMember>
  </entry>""",
                "exp": {
                    "cluster": [
                        {
                            "cluster_id": "uniref:UniRef50_UPI00358F51CD",
                            "name": "Cluster: LOW QUALITY PROTEIN: titin",
                            "cluster_type": "Protein",
                            "protocol": PROTOCOL,
                            "description": None,
                        }
                    ],
                    "entity": [
                        {
                            "entity_id": "uniref:UniRef50_UPI00358F51CD",
                            "entity_type": "Cluster",
                            "data_source": "UniRef",
                            "data_source_entity_id": "UniRef50_UPI00358F51CD",
                            "data_source_updated": "2024-11-27",
                            "updated": NOW,
                        }
                    ],
                    "clustermember": [
                        {
                            "entity_id": "uniparc:UPI00358F51CD",
                            "cluster_id": "uniref:UniRef50_UPI00358F51CD",
                            "is_representative": True,
                            "is_seed": True,
                        },
                    ],
                },
            },
            id="single UniParc rep",
        ),
        pytest.param(
            {
                "xml": """
  <entry xmlns="http://uniprot.org/uniref" id="UniRef20_CRAP" updated="1066-10-10">
    <property type="member count" value="0"/>
    <representativeMember>
      <dbReference type="UniParc ID" id="UPI00358F51CD">
        <property type="length" value="2"/>
        <property type="isSeed" value="true"/>
      </dbReference>
      <sequence length="2" checksum="48729625616C010E">SNIP</sequence>
    </representativeMember>
  </entry>
""",
                "exp": {
                    "cluster": [
                        {
                            "cluster_id": "uniref:UniRef20_CRAP",
                            "name": None,
                            "cluster_type": "Protein",
                            "protocol": PROTOCOL,
                            "description": None,
                        }
                    ],
                    "entity": [
                        {
                            "entity_id": "uniref:UniRef20_CRAP",
                            "entity_type": "Cluster",
                            "data_source": "UniRef",
                            "data_source_entity_id": "UniRef20_CRAP",
                            "data_source_updated": "1066-10-10",
                            "updated": NOW,
                        }
                    ],
                    "clustermember": [
                        {
                            "entity_id": "uniparc:UPI00358F51CD",
                            "cluster_id": "uniref:UniRef20_CRAP",
                            "is_representative": True,
                            "is_seed": True,
                        },
                    ],
                },
            },
            id="single UniParc rep, no name",
        ),
    ],
)
def test_parse_whole_uniref_entry(entry: dict[str, Any]) -> None:
    """Ensure that parsing the whole entry results in the expected tables being yielded."""
    parsed = parse_uniref_entry(fromstring(entry["xml"]), NOW, "/path/to/file", PROTOCOL)
    assert parsed["entity"] == entry["exp"]["entity"]
    assert parsed["cluster"] == entry["exp"]["cluster"]
    assert parsed["clustermember"] == entry["exp"]["clustermember"]
    assert parsed["entity_x_source_file"] == [
        {
            "source_file": "/path/to/file",
            **{k: entry["exp"]["entity"][0][k] for k in ["entity_id", "data_source"]},
        }
    ]


@pytest.mark.parametrize(
    ("db", "acc", "expected"),
    [
        # known prefixes are translated
        ("UniProtKB ID", "TITIN_MOUSE", "uniprot_name:TITIN_MOUSE"),
        ("UniProtKB accession", "A2ASS6", "uniprot:A2ASS6"),
        ("UniParc ID", "UPI0000D77B45", "uniparc:UPI0000D77B45"),
        ("UniRef90 ID", "UniRef90_A2ASS6", "uniref:UniRef90_A2ASS6"),
        ("UniRef50 ID", "UniRef50_A2ASS6", "uniref:UniRef50_A2ASS6"),
        ("UniRef100 ID", "UniRef100_A2ASS6", "uniref:UniRef100_A2ASS6"),
        ("NCBI taxonomy", "10090", "NCBITaxon:10090"),
        # unknown prefixes fall through unchanged
        ("SomeUnknownDB", "XYZ123", "SomeUnknownDB:XYZ123"),
    ],
)
def test_generate_dbxref(db: str, acc: str, expected: str) -> None:
    """Test that dbxrefs use the BioRegistry-style prefixes, falling back to the raw db name."""
    assert generate_dbxref(db, acc) == expected


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        # the http uniref namespace declaration is stripped
        (
            f"<entry xmlns='{UNIREF_URL}'><name>Cluster: x</name></entry>",
            "<entry><name>Cluster: x</name></entry>",
        ),
        # the https variant of the namespace is also stripped
        (
            "<entry xmlns='https://uniprot.org/uniref'><name>Cluster: y</name></entry>",
            "<entry><name>Cluster: y</name></entry>",
        ),
        # the xsi namespace declaration is stripped too
        (
            f'<entry xmlns="{UNIREF_URL}" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            "<name>Cluster: z</name></entry>",
            "<entry><name>Cluster: z</name></entry>",
        ),
        # self-closing element
        (f"<entry xmlns='{UNIREF_URL}'/>", "<entry/>"),
    ],
)
def test_dump_xml_element(xml: str, expected: str) -> None:
    """Test that namespace declarations are stripped when dumping XML to a string."""
    dumped = dump_xml_element(fromstring(xml, parser=XMLParser(remove_blank_text=True)))
    assert dumped == expected


def test_extract_cluster_no_id() -> None:
    """An entry with no id attribute yields empty cluster and entity dicts."""
    xml = f"<entry xmlns='{UNIREF_URL}'><name>Cluster: no id</name></entry>"
    assert extract_cluster(fromstring(xml)) == ({}, {})


def test_parse_uniref_entry_no_cluster_id() -> None:
    """An entry with no id attribute is reported as a _parse_error rather than parsed."""
    xml = f"<entry xmlns='{UNIREF_URL}'><name>Cluster: no id</name></entry>"
    result = parse_uniref_entry(
        fromstring(xml, parser=XMLParser(remove_blank_text=True)),
        NOW,
        FILE_PATH,
        PROTOCOL,
    )

    assert list(result.keys()) == ["_parse_error"]
    error = result["_parse_error"][0]
    assert error["error"] == "No cluster ID found in entry"
    assert error["source_file"] == FILE_PATH
    # the raw XML (namespace-stripped) is retained for later inspection
    assert error["xml"] == "<entry><name>Cluster: no id</name></entry>"


FIXTURE_CHUNK = Path("tests/fixtures/uniref/uniref_chunk_00001.xml")


@pytest.mark.skipif(not FIXTURE_CHUNK.exists(), reason="UniRef XML fixture not available")
def test_parse_uniref_entry_integration() -> None:
    """Integration test: stream real UniRef entries from a fixture file and parse each one.

    Exercises the parser end-to-end against genuine UniRef XML (as produced upstream),
    ensuring every entry parses without error and produces the expected schema tables.
    """
    # NOTE: stream_xml_file clears each element right after yielding it, so entries must be
    # parsed inside the streaming loop rather than collected into a list first.
    entry_count = 0
    seen_cluster_ids: set[str] = set()
    for entry in stream_xml_file(FIXTURE_CHUNK, ENTRY_XML_TAG):
        entry_count += 1
        parsed = parse_uniref_entry(entry, NOW, str(FIXTURE_CHUNK), PROTOCOL)

        # every real entry should parse cleanly
        assert "_parse_error" not in parsed, f"unexpected parse error: {parsed.get('_parse_error')}"

        # cluster + entity rows present and consistent
        assert len(parsed["cluster"]) == 1
        assert len(parsed["entity"]) == 1
        cluster = parsed["cluster"][0]
        entity = parsed["entity"][0]
        cluster_id = cluster["cluster_id"]
        seen_cluster_ids.add(cluster_id)

        assert cluster_id.startswith("uniref:")
        assert cluster["cluster_type"] == "Protein"
        assert cluster["protocol"] == PROTOCOL
        assert entity["entity_id"] == cluster_id
        assert entity["entity_type"] == "Cluster"
        assert entity["data_source"] == "UniRef"
        assert entity["updated"] == NOW

        # at least the representative member is present, and all members reference this cluster
        assert parsed["clustermember"], f"no cluster members for {cluster_id}"
        assert any(m["is_representative"] for m in parsed["clustermember"])
        assert all(m["cluster_id"] == cluster_id for m in parsed["clustermember"])

        # source-file linkage present
        assert parsed["entity_x_source_file"] == [
            {"entity_id": cluster_id, "data_source": "UniRef", "source_file": str(FIXTURE_CHUNK)}
        ]

    assert entry_count > 0, "expected at least one entry in the fixture"
    # the three distinct clusters in the fixture should each be parsed exactly once
    assert seen_cluster_ids == {
        "uniref:UniRef50_A0AB34IYJ6",
        "uniref:UniRef50_A0A410P257",
        "uniref:UniRef50_UPI00358F51CD",
    }
