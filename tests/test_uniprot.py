"""

This file uses pytest to provide parameterized and functional tests for all major
UniProt parsing utility functions, ensuring correct parsing and transformation of
UniProt XML into structured CDM records.

Coverage:
    - generate_cdm_id:     Stable CDM entity ID from accession
    - build_datasource_record: Datasource provenance and metadata
    - parse_identifiers:   UniProt accessions to identifier records
    - parse_names:         Protein names (top-level, recommended, alternative)
    - parse_protein_info:  EC numbers, existence evidence, sequences, etc.
    - parse_evidence_map:  Evidence element mapping and supporting objects
    - parse_associations:  Biological and database associations (taxonomy, PDB, Rhea, ChEBI)
    - parse_publications:  Supported literature references (PMID, DOI, etc.)
    - parse_uniprot_entry: Full record parsing, all fields together

How to run in the terminal:
    PYTHONPATH=src pytest tests/test_uniprot.py

"""

import re
import pytest
import xml.etree.ElementTree as ET
import datetime
from parsers.uniprot import (
    generate_cdm_id,
    build_datasource_record,
    parse_identifiers,
    parse_names,
    parse_protein_info,
    parse_evidence_map,
    parse_associations,
    parse_publications,
    parse_uniprot_entry
)

# Regular expression to validate UUID format
UUID_PATTERN = re.compile(
    r"^[a-f0-9]{8}-[a-f0-9]{4}-[1-5][a-f0-9]{3}-[89ab][a-f0-9]{3}-[a-f0-9]{12}$",
    re.IGNORECASE
)

@pytest.mark.parametrize("n", range(5))
def test_generate_cdm_id_format(n):
    uuid = generate_cdm_id()
    assert uuid.startswith("CDM:")
    uuid_str = uuid[4:]
    assert UUID_PATTERN.match(uuid_str), f"{uuid_str} is not a valid UUID"


## build_datasource_record ##
def test_build_datasource_record():
    url = "https://example.com/uniprot.xml.gz"
    record = build_datasource_record(url)
    assert isinstance(record, dict)
    assert set(record.keys()) == {"name", "source", "url", "accessed", "version"}
    assert record["name"] == "UniProt import"
    assert record["source"] == "UniProt"
    assert record["url"] == url

    # check accessed
    accessed_dt = datetime.datetime.fromisoformat(record["accessed"])
    now = datetime.datetime.now(datetime.timezone.utc)
    delta = abs((now - accessed_dt).total_seconds())
    assert delta < 10
    assert record["version"] == 115


@pytest.mark.parametrize("bad_url", [None, ""])
def test_build_datasource_record_bad(bad_url):
    record = build_datasource_record(bad_url)
    assert record["url"] == bad_url


## parse_identifiers function test ##
@pytest.mark.parametrize(
    "xml_str,cdm_id,expected",
    [
        ### multiple accessions, expect two dict, every dic use the same cdm_id
        ### identifier according to <accession> number
        (
            """
            <entry xmlns="https://uniprot.org/uniprot">
                <accession>Q9V2L2</accession>
                <accession>G8ZFP4</accession>
            </entry>
            """,
            "CDM:001",
            [
                {
                    "entity_id": "CDM:001",
                    "identifier": "UniProt:Q9V2L2",
                    "source": "UniProt",
                    "description": "UniProt accession",
                },
                {
                    "entity_id": "CDM:001",
                    "identifier": "UniProt:G8ZFP4",
                    "source": "UniProt",
                    "description": "UniProt accession",
                },
            ],
        ),
        ### Use single accession
        (
            """
            <entry xmlns="https://uniprot.org/uniprot">
                <accession>X00001</accession>
            </entry>
            """,
            "CDM:002",
            [
                {
                    "entity_id": "CDM:002",
                    "identifier": "UniProt:X00001",
                    "source": "UniProt",
                    "description": "UniProt accession",
                }
            ],
        ),
        ### No accession
        (
            """
            <entry xmlns="https://uniprot.org/uniprot">
            </entry>
            """,
            "CDM:003",
            [],
        ),
    ],
)
def test_parse_identifiers(xml_str, cdm_id, expected):
    """
    This approach ensures that parse_identifiers correctly parses and structures identifier data.

    The parsed Element object and the provided CDM_id are passed to the parse_identifiers funtion.
    The function is expected to extract all relevant identifier information from the XML and return list of dict.

    The test compares the result output with the predefined expected result using an assert statement.

    """
    entry = ET.fromstring(xml_str)
    result = parse_identifiers(entry, cdm_id)
    assert result == expected


"""
    This parameterized pytest function tests the correctness of the parse_names function for various UniProt XML entry scenarios. 

    XML string representing a UniProt entry with different protein names: 
    top-level <name> 
    recommended names,
    alternative names,
    combinations, 
    no names

    cdm_id: CDM entry ID 

    Output: 
    A list of name records with their metadata 

"""


## parse_names function test ##
@pytest.mark.parametrize(
    "xml_str, cdm_id, expected",
    [
        # Only top-level <name>
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <name>MainProteinName</name>
               </entry>""",
            "CDM:001",
            [
                {
                    "entity_id": "CDM:001",
                    "name": "MainProteinName",
                    "description": "UniProt protein name",
                    "source": "UniProt",
                }
            ],
        ),
        # RecommendedName (fullName and shortName)
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <protein>
                     <recommendedName>
                        <fullName>RecFullName</fullName>
                        <shortName>RecShort</shortName>
                     </recommendedName>
                   </protein>
               </entry>""",
            "CDM:002",
            [
                {
                    "entity_id": "CDM:002",
                    "name": "RecFullName",
                    "description": "UniProt recommended full name",
                    "source": "UniProt",
                },
                {
                    "entity_id": "CDM:002",
                    "name": "RecShort",
                    "description": "UniProt recommended short name",
                    "source": "UniProt",
                },
            ],
        ),
        # AlternativeName (fullName and shortName)
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <protein>
                     <alternativeName>
                        <fullName>AltFullName1</fullName>
                        <shortName>AltShort1</shortName>
                     </alternativeName>
                     <alternativeName>
                        <fullName>AltFullName2</fullName>
                     </alternativeName>
                   </protein>
               </entry>""",
            "CDM:003",
            [
                {
                    "entity_id": "CDM:003",
                    "name": "AltFullName1",
                    "description": "UniProt alternative full name",
                    "source": "UniProt",
                },
                {
                    "entity_id": "CDM:003",
                    "name": "AltShort1",
                    "description": "UniProt alternative short name",
                    "source": "UniProt",
                },
                {
                    "entity_id": "CDM:003",
                    "name": "AltFullName2",
                    "description": "UniProt alternative full name",
                    "source": "UniProt",
                },
            ],
        ),
        # Mixed: top-level <name> and <protein>
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <name>TopLevel</name>
                   <protein>
                     <recommendedName>
                        <fullName>MixedFull</fullName>
                     </recommendedName>
                   </protein>
               </entry>""",
            "CDM:004",
            [
                {
                    "entity_id": "CDM:004",
                    "name": "TopLevel",
                    "description": "UniProt protein name",
                    "source": "UniProt",
                },
                {
                    "entity_id": "CDM:004",
                    "name": "MixedFull",
                    "description": "UniProt recommended full name",
                    "source": "UniProt",
                },
            ],
        ),
        # No names at all
        (
            """<entry xmlns="https://uniprot.org/uniprot">
               </entry>""",
            "CDM:005",
            [],
        ),
    ],
)
def test_parse_names(xml_str, cdm_id, expected):
    entry = ET.fromstring(xml_str)
    result = parse_names(entry, cdm_id)
    assert result == expected


"""

    This test ensures parse_protein_info works correctly for different combinations of data 
    Including cases with no protein info, sequence only, existence only or EC numbers
    
    This approach thoroughly validates that parse_protein_info can accurately extract, combine and structure metadata field. 
    
    Include: 
    EC Number, 
    existence evidence, 
    sequence

"""


## parse_protein_info function test ##
@pytest.mark.parametrize(
    "xml_str, cdm_id, expected",
    [
        # There are multiple ecNumbers under the recommend names
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <protein>
                  <recommendedName>
                    <ecNumber>1.2.3.4</ecNumber>
                    <ecNumber>5.6.7.8</ecNumber>
                  </recommendedName>
                </protein>
            </entry>""",
            "CDM:001",
            {"ec_numbers": ["1.2.3.4", "5.6.7.8"]},
        ),
        # alternativeName has EC Number
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <protein>
                  <alternativeName>
                    <ecNumber>3.3.3.3</ecNumber>
                  </alternativeName>
                </protein>
            </entry>""",
            "CDM:002",
            {"ec_numbers": ["3.3.3.3"]},
        ),
        # If have both proteinExistence evidence and existence
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <proteinExistence type="evidence at protein level"/>
            </entry>""",
            "CDM:003",
            {
                "protein_id": "CDM:003",
                "evidence_for_existence": "evidence at protein level",
            },
        ),
        # Sequence only
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <sequence length="357" mass="12345" checksum="ABCD" modified="2024-05-21" version="2">
                MAGNLSKVAAVSGVAAAVLGK
                </sequence>
            </entry>""",
            "CDM:004",
            {
                "length": "357",
                "mass": "12345",
                "checksum": "ABCD",
                "modified": "2024-05-21",
                "sequence_version": "2",
                "sequence": "MAGNLSKVAAVSGVAAAVLGK",
            },
        ),
        # Combine with three elements: proteinExistence, sequence and ecNumbers
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <protein>
                  <recommendedName>
                    <ecNumber>3.3.3.3</ecNumber>
                  </recommendedName>
                  <alternativeName>
                    <ecNumber>8.8.8.8</ecNumber>
                  </alternativeName>
                </protein>
                <proteinExistence type="evidence at transcript level"/>
                <sequence length="10" mass="1000" checksum="XYZ" modified="2021-12-01" version="1">
                  MKTLLTGAAT
                </sequence>
            </entry>""",
            "CDM:005",
            {
                "ec_numbers": ["3.3.3.3", "8.8.8.8"],
                "protein_id": "CDM:005",
                "evidence_for_existence": "evidence at transcript level",
                "length": "10",
                "mass": "1000",
                "checksum": "XYZ",
                "modified": "2021-12-01",
                "sequence_version": "1",
                "sequence": "MKTLLTGAAT",
            },
        ),
        # return None
        ("""<entry xmlns="https://uniprot.org/uniprot"></entry>""", "CDM:006", None),
    ],
)
def test_parse_protein_info(xml_str, cdm_id, expected):
    entry = ET.fromstring(xml_str)
    result = parse_protein_info(entry, cdm_id)
    assert result == expected


"""

    This parameterized pytest function verifies the behavior of the parse_evidence_map function
    for different UniProt XML entry structures involving evidence elements.

    xml_str: Simulates a UniProt entry with various <evidence> and <source> sub-structures, 
    including cases with multiple evidence elements, missing sources, or no evidence at all.

    expected: A dictionary mapping evidence keys to their extracted details—such as evidence type, 
    supporting objects, and publication references.

    Ensure parse_evidence_map: 
    Accurately extract evidence keys and types
    Correctly classify supporting objects and publication references
    Handle entries with absent sources or evidence elements
    Represent all relevant evidence metadata in the required structure

"""


## parse_evidence_map function test ##
@pytest.mark.parametrize(
    "xml_str, expected",
    [
        # Single evidence，include PubMed and supporting object
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <evidence key="1" type="ECO:0000255">
                  <source>
                    <dbReference type="PubMed" id="123456"/>
                    <dbReference type="Ensembl" id="ENSG00001"/>
                  </source>
                </evidence>
            </entry>""",
            {
                "1": {
                    "evidence_type": "ECO:0000255",
                    "supporting_objects": ["Ensembl:ENSG00001"],
                    "publications": ["PMID:123456"],
                }
            },
        ),
        # multiple evidences
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <evidence key="E1" type="ECO:0000313">
                  <source>
                    <dbReference type="PubMed" id="654321"/>
                  </source>
                </evidence>
                <evidence key="E2" type="ECO:0000250">
                  <source>
                    <dbReference type="PDB" id="2N7Q"/>
                  </source>
                </evidence>
            </entry>""",
            {
                "E1": {
                    "evidence_type": "ECO:0000313",
                    "supporting_objects": None,
                    "publications": ["PMID:654321"],
                },
                "E2": {
                    "evidence_type": "ECO:0000250",
                    "supporting_objects": ["PDB:2N7Q"],
                    "publications": None,
                },
            },
        ),
        # no source
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <evidence key="X1" type="ECO:9999999"/>
            </entry>""",
            {
                "X1": {
                    "evidence_type": "ECO:9999999",
                    "supporting_objects": None,
                    "publications": None,
                }
            },
        ),
        # no evidence
        (
            """<entry xmlns="https://uniprot.org/uniprot">
            </entry>""",
            {},
        ),
        # one evidence with multiple supporting objects
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <evidence key="K1" type="ECO:0000269">
                  <source>
                    <dbReference type="Ensembl" id="ENS1"/>
                    <dbReference type="RefSeq" id="RS123"/>
                  </source>
                </evidence>
            </entry>""",
            {
                "K1": {
                    "evidence_type": "ECO:0000269",
                    "supporting_objects": ["Ensembl:ENS1", "RefSeq:RS123"],
                    "publications": None,
                }
            },
        ),
    ],
)
def test_parse_evidence_map(xml_str, expected):
    entry = ET.fromstring(xml_str)
    result = parse_evidence_map(entry)
    assert result == expected


"""

    xml_strings: models a UniProt entry with different types of possible associations
    cdm_id: uniquely identifies the protein being parsed
    evidence_map:  supplies external evidence metadata for associations
    expected: list of association dictionaries 

    Arg: 
    The function correctly links proteins to organism taxonomy.
	Cross-references are properly included, evidence metadata is correctly merged.
	Associations derived from catalytic activity and cofactor comments are correctly generated.
	All combinations and edge cases are handled robustly.

"""


## parse_associations function test ##
@pytest.mark.parametrize(
    "xml_str, cdm_id, evidence_map, expected",
    [
        # organism association（NCBI Taxonomy dbReference）
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <organism>
                      <dbReference type="NCBI Taxonomy" id="9606"/>
                   </organism>
            </entry>""",
            "CDM:1",
            {},
            [{"subject": "CDM:1", "object": "NCBITaxon:9606"}],
        ),
        # dbReference with evidence key
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <dbReference type="PDB" id="2N7Q" evidence="E1"/>
            </entry>""",
            "CDM:2",
            {
                "E1": {
                    "evidence_type": "ECO:0000250",
                    "supporting_objects": ["Ensembl:ENS1"],
                    "publications": ["PMID:1234"],
                }
            },
            [
                {
                    "subject": "CDM:2",
                    "object": "PDB:2N7Q",
                    "evidence_type": "ECO:0000250",
                    "supporting_objects": ["Ensembl:ENS1"],
                    "publications": ["PMID:1234"],
                }
            ],
        ),
        # comment catalytic activity (reaction) with evidence key
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <comment type="catalytic activity">
                     <reaction evidence="E2">
                        <dbReference type="Rhea" id="12345"/>
                     </reaction>
                   </comment>
            </entry>""",
            "CDM:3",
            {
                "E2": {
                    "evidence_type": "ECO:0000313",
                    "publications": ["PMID:2222"],
                }
            },
            [
                {
                    "subject": "CDM:3",
                    "predicate": "catalyzes",
                    "object": "Rhea:12345",
                    "evidence_type": "ECO:0000313",
                    "publications": ["PMID:2222"],
                }
            ],
        ),
        # Comment cofactor without evidence
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <comment type="cofactor">
                     <cofactor>
                       <dbReference type="ChEBI" id="CHEBI:15377"/>
                     </cofactor>
                   </comment>
            </entry>""",
            "CDM:4",
            {},
            [
                {
                    "subject": "CDM:4",
                    "predicate": "requires_cofactor",
                    "object": "ChEBI:CHEBI:15377",
                }
            ],
        ),
        # Several relevant relationship（with organism and dbReference）
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                   <organism>
                      <dbReference type="NCBI Taxonomy" id="562"/>
                   </organism>
                   <dbReference type="RefSeq" id="NP_414543"/>
            </entry>""",
            "CDM:5",
            {},
            [
                {"subject": "CDM:5", "object": "NCBITaxon:562"},
                {"subject": "CDM:5", "object": "RefSeq:NP_414543"},
            ],
        ),
        # if it is empty entry, return to []
        ("""<entry xmlns="https://uniprot.org/uniprot"></entry>""", "CDM:6", {}, []),
    ],
)
def test_parse_associations(xml_str, cdm_id, evidence_map, expected):
    entry = ET.fromstring(xml_str)
    result = parse_associations(entry, cdm_id, evidence_map)
    assert result == expected


"""

    xml_str: Uniprot entry include <reference>, <citation>, 
    Refer: PubMed, DOI, GeneBank, DDBJ, EMBL 

    Output: List of publication identifier 

    Arg: 
    Extract publication of references
	Recognize and format database types ( with prefixing “PMID:”, “DOI:”)
	Handle entries with multiple or mixed publication types
	Return an empty list if no publication data.
    
"""


## parse_publications function test ##
@pytest.mark.parametrize(
    "xml_str, expected",
    [
        # Single PubMed
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <reference>
                  <citation>
                    <dbReference type="PubMed" id="12345"/>
                  </citation>
                </reference>
            </entry>""",
            ["PMID:12345"],
        ),
        # Multiple types include (PubMed, DOI, GenBank)
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <reference>
                  <citation>
                    <dbReference type="PubMed" id="55555"/>
                    <dbReference type="DOI" id="10.1000/j.jmb.2020.01.001"/>
                    <dbReference type="GenBank" id="AB123456"/>
                  </citation>
                </reference>
            </entry>""",
            ["PMID:55555", "DOI:10.1000/j.jmb.2020.01.001"],
        ),
        # Multiple references
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <reference>
                  <citation>
                    <dbReference type="DOI" id="10.1000/jmb.123456"/>
                  </citation>
                </reference>
                <reference>
                  <citation>
                    <dbReference type="PubMed" id="98765"/>
                  </citation>
                </reference>
            </entry>""",
            ["DOI:10.1000/jmb.123456", "PMID:98765"],
        ),
        # dbReference: DDBJ and EMBL
        (
            """<entry xmlns="https://uniprot.org/uniprot">
                <reference>
                  <citation>
                    <dbReference type="DDBJ" id="BA000001"/>
                    <dbReference type="EMBL" id="AB987654"/>
                  </citation>
                </reference>
            </entry>""",
            [],
        ),
        # no publication
        ("""<entry xmlns="https://uniprot.org/uniprot"></entry>""", []),
    ],
)
def test_parse_publications(xml_str, expected):
    entry = ET.fromstring(xml_str)
    result = parse_publications(entry)
    assert result == expected


## parse_uniprot_entry function test ##
@pytest.mark.parametrize(
    "xml_str, datasource_name, prev_created",
    [
        (
            """
            <entry xmlns="https://uniprot.org/uniprot" created="2020-01-01" modified="2021-01-01" version="3">
                <accession>P12345</accession>
                <name>ProteinX</name>
                <protein>
                    <recommendedName>
                        <fullName>ProteinX Full Name</fullName>
                    </recommendedName>
                </protein>
                <organism>
                    <dbReference type="NCBI Taxonomy" id="9606"/>
                </organism>
                <reference>
                    <citation>
                        <dbReference type="PubMed" id="99999"/>
                    </citation>
                </reference>
            </entry>
            """,
            "UniProt import",
            None,
        ),
    ],
)

def test_parse_uniprot_entry(xml_str, datasource_name, prev_created):
    import xml.etree.ElementTree as ET
    entry = ET.fromstring(xml_str)
    cdm_id = generate_cdm_id()   

    current_timestamp = "2024-07-17T13:00:00Z"

    record = parse_uniprot_entry(entry, cdm_id, current_timestamp, datasource_name, prev_created)

    entity = record["entity"]
    assert entity["entity_type"] == "protein"
    assert entity["data_source"] == datasource_name
    assert entity["version"] == "3"
    assert entity["uniprot_created"] == "2020-01-01"
    assert entity["uniprot_modified"] == "2021-01-01"
    assert entity["entity_id"].startswith("CDM:")

    # identifiers/names/associations/publications
    assert isinstance(record["identifiers"], list)
    assert isinstance(record["names"], list)
    assert isinstance(record["associations"], list)
    assert isinstance(record["publications"], list)
