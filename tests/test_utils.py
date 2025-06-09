import re 
from datetime import date 
import pytest 
import xml.etree.ElementTree as ET

from uniprot_to_cdm import generate_unique_id, build_datasource_record

"""
    
@pytest.mark.parametrize is pytest's parameterization mechanism. 
It can automatically call the same test function multiple times with a set of parameter data.

"""

# Test generate id, ensure it is unqiue id 
@pytest.mark.parametrize("run", range(10))

def test_generate_unqiue_id(run):
    random_id = generate_unique_id()
    # Check the format as `CDM:00000000-0000-0000-12345678`
    assert random_id.startswith("CDM:")
    # Check if there is function uuid.uuid4
    check_uuid = random_id.split(":")[1]
    assert re.fullmatch(r"[0-9a-fA-F-]{36}", check_uuid)

# Test build datasource record 
def test_build_datasource_record():
    result = build_datasource_record()
    assert isinstance(result, dict)
    assert result['name'] == 'UniProt archaea'
    assert result['source'] == 'UniProt'
    assert 'url' in result
    # check the date time 
    assert result['accessed'] == date.today().strftime('%Y-%m-%d')
    assert isinstance(result['version'], int)


from uniprot_to_cdm import parse_identifiers
@pytest.mark.parametrize(
    ## Three parameters for test this function 
    "xml_str,cdm_id,expected", [
        ### multiple accessions, expect two dict, every dic use the same cdm_id
        ### identifier according to <accession> number 
        ### List all accessions 
        (
            '''
            <entry xmlns="https://uniprot.org/uniprot">
                <accession>Q9V2L2</accession>
                <accession>G8ZFP4</accession>
            </entry>
            ''',
            "CDM:001",
            [
                {
                    'entity_id': "CDM:001",
                    'identifier': "UniProt:Q9V2L2",
                    'source': "UniProt",
                    'description': "UniProt accession"
                },
                {
                    'entity_id': "CDM:001",
                    'identifier': "UniProt:G8ZFP4",
                    'source': "UniProt",
                    'description': "UniProt accession"
                }
            ]
        ), 
        ### Use single accession 
        (
            '''
            <entry xmlns="https://uniprot.org/uniprot">
                <accession>X00001</accession>
            </entry>
            ''',
            "CDM:002",
            [
                {
                    'entity_id': "CDM:002",
                    'identifier': "UniProt:X00001",
                    'source': "UniProt",
                    'description': "UniProt accession"
                }
            ]
        ),
        ### No accession 
        (
            '''
            <entry xmlns="https://uniprot.org/uniprot">
            </entry>
            ''',
            "CDM:003",
            []
        ),
    ]
)

def test_parse_identifiers(xml_str, cdm_id, expected):
    """
    This approach ensures that parse_identifiers correctly parses and structures identifier data. 

    The parsed Element object and the provided CDM_id are passed to the parse_identifiers funtion. 
    The function is expected to extract all relevant identifier information from the XML and return list of dict. 
    
    The test compares the result output with the predefined expected result using an assert statement.
  
    """
    # parse ET.fromstring() to Element, parsing xml strings into Element objects
    entry = ET.fromstring(xml_str)
    result = parse_identifiers(entry, cdm_id)
    assert result == expected



from uniprot_to_cdm import parse_names
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
@pytest.mark.parametrize(
    "xml_str, cdm_id, expected",
    [
        # Only top-level <name>
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <name>MainProteinName</name>
               </entry>''',
            "CDM:001",
            [
                {
                    'entity_id': "CDM:001",
                    'name': "MainProteinName",
                    'description': "UniProt protein name",
                    'source': "UniProt"
                }
            ]
        ),
        # RecommendedName (fullName and shortName)
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <protein>
                     <recommendedName>
                        <fullName>RecFullName</fullName>
                        <shortName>RecShort</shortName>
                     </recommendedName>
                   </protein>
               </entry>''',
            "CDM:002",
            [
                {
                    'entity_id': "CDM:002",
                    'name': "RecFullName",
                    'description': "UniProt recommended full name",
                    'source': "UniProt"
                },
                {
                    'entity_id': "CDM:002",
                    'name': "RecShort",
                    'description': "UniProt recommended short name",
                    'source': "UniProt"
                }
            ]
        ),
        # AlternativeName (fullName and shortName)
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <protein>
                     <alternativeName>
                        <fullName>AltFullName1</fullName>
                        <shortName>AltShort1</shortName>
                     </alternativeName>
                     <alternativeName>
                        <fullName>AltFullName2</fullName>
                     </alternativeName>
                   </protein>
               </entry>''',
            "CDM:003",
            [
                {
                    'entity_id': "CDM:003",
                    'name': "AltFullName1",
                    'description': "UniProt alternative full name",
                    'source': "UniProt"
                },
                {
                    'entity_id': "CDM:003",
                    'name': "AltShort1",
                    'description': "UniProt alternative short name",
                    'source': "UniProt"
                },
                {
                    'entity_id': "CDM:003",
                    'name': "AltFullName2",
                    'description': "UniProt alternative full name",
                    'source': "UniProt"
                }
            ]
        ),
        # Mixed: top-level <name> and <protein>
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <name>TopLevel</name>
                   <protein>
                     <recommendedName>
                        <fullName>MixedFull</fullName>
                     </recommendedName>
                   </protein>
               </entry>''',
            "CDM:004",
            [
                {
                    'entity_id': "CDM:004",
                    'name': "TopLevel",
                    'description': "UniProt protein name",
                    'source': "UniProt"
                },
                {
                    'entity_id': "CDM:004",
                    'name': "MixedFull",
                    'description': "UniProt recommended full name",
                    'source': "UniProt"
                }
            ]
        ),
        # No names at all
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
               </entry>''',
            "CDM:005",
            []
        ),
    ]
)

def test_parse_names(xml_str, cdm_id, expected):
    entry = ET.fromstring(xml_str)
    result = parse_names(entry, cdm_id)
    assert result == expected



from uniprot_to_cdm import parse_protein_info
"""

    This test ensures parse_protein_info works correctly for different combinations of data 
    Including cases with no protein info, sequence only, existence only or EC numbers
    
    This approach thoroughly validates that parse_protein_info can accurately extract, combine and structure metadata field. 
    
    Include: 
    EC Number, 
    existence evidence, 
    sequence

"""

@pytest.mark.parametrize(
    "xml_str, cdm_id, expected",
    [
        # There are multiple ecNumbers under the recommend names
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <protein>
                  <recommendedName>
                    <ecNumber>1.2.3.4</ecNumber>
                    <ecNumber>5.6.7.8</ecNumber>
                  </recommendedName>
                </protein>
            </entry>''',
            "CDM:001",
            {
                "ec_numbers": ["1.2.3.4", "5.6.7.8"]
            }
        ),

        # alternativeName has EC Number
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <protein>
                  <alternativeName>
                    <ecNumber>3.3.3.3</ecNumber>
                  </alternativeName>
                </protein>
            </entry>''',
            "CDM:002",
            {
                "ec_numbers": ["3.3.3.3"]
            }
        ),

        # If have both proteinExistence evidence and existence 
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <proteinExistence type="evidence at protein level"/>
            </entry>''',
            "CDM:003",
            {
                "protein_id": "CDM:003",
                "evidence_for_existence": "evidence at protein level"
            }
        ),

        # Sequence only
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <sequence length="357" mass="12345" checksum="ABCD" modified="2024-05-21" version="2">
                MAGNLSKVAAVSGVAAAVLGK
                </sequence>
            </entry>''',
            "CDM:004",
            {
                "length": "357",
                "mass": "12345",
                "checksum": "ABCD",
                "modified": "2024-05-21",
                "sequence_version": "2",
                "sequence": "MAGNLSKVAAVSGVAAAVLGK"
            }
        ),

        # Combine with three elements: proteinExistence, sequence and ecNumbers
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
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
            </entry>''',
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
                "sequence": "MKTLLTGAAT"
            }
        ),
        # return None
        (
            '''<entry xmlns="https://uniprot.org/uniprot"></entry>''',
            "CDM:006",
            None
        ),
    ]
)

def test_parse_protein_info(xml_str, cdm_id, expected):
    entry = ET.fromstring(xml_str)
    result = parse_protein_info(entry, cdm_id)
    assert result == expected



from uniprot_to_cdm import parse_evidence_map
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

@pytest.mark.parametrize(
    "xml_str, expected",
    [
        # Single evidence，include PubMed and supporting object
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <evidence key="1" type="ECO:0000255">
                  <source>
                    <dbReference type="PubMed" id="123456"/>
                    <dbReference type="Ensembl" id="ENSG00001"/>
                  </source>
                </evidence>
            </entry>''',
            {
                "1": {
                    "evidence_type": "ECO:0000255",
                    "supporting_objects": ["Ensembl:ENSG00001"],
                    "publications": ["PMID:123456"]
                }
            }
        ),

        # multiple evidences
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
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
            </entry>''',
            {
                "E1": {
                    "evidence_type": "ECO:0000313",
                    "supporting_objects": None,
                    "publications": ["PMID:654321"]
                },
                "E2": {
                    "evidence_type": "ECO:0000250",
                    "supporting_objects": ["PDB:2N7Q"],
                    "publications": None
                }
            }
        ),

        # no source
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <evidence key="X1" type="ECO:9999999"/>
            </entry>''',
            {
                "X1": {
                    "evidence_type": "ECO:9999999",
                    "supporting_objects": None,
                    "publications": None
                }
            }
        ),

        # no evidence 
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
            </entry>''',
            {}
        ),

        # one evidence with multiple supporting objects 
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <evidence key="K1" type="ECO:0000269">
                  <source>
                    <dbReference type="Ensembl" id="ENS1"/>
                    <dbReference type="RefSeq" id="RS123"/>
                  </source>
                </evidence>
            </entry>''',
            {
                "K1": {
                    "evidence_type": "ECO:0000269",
                    "supporting_objects": ["Ensembl:ENS1", "RefSeq:RS123"],
                    "publications": None
                }
            }
        ),
    ]
)
def test_parse_evidence_map(xml_str, expected):
    entry = ET.fromstring(xml_str)
    result = parse_evidence_map(entry)
    assert result == expected



from uniprot_to_cdm import parse_associations

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

@pytest.mark.parametrize(
    "xml_str, cdm_id, evidence_map, expected",
    [
        # organism association（NCBI Taxonomy dbReference）
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <organism>
                      <dbReference type="NCBI Taxonomy" id="9606"/>
                   </organism>
            </entry>''',
            "CDM:1",
            {},
            [
                {'subject': "CDM:1", 'object': "NCBITaxon:9606"}
            ]
        ),

        # dbReference with evidence key
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <dbReference type="PDB" id="2N7Q" evidence="E1"/>
            </entry>''',
            "CDM:2",
            {
                "E1": {
                    "evidence_type": "ECO:0000250",
                    "supporting_objects": ["Ensembl:ENS1"],
                    "publications": ["PMID:1234"]
                }
            },
            [
                {
                    "subject": "CDM:2",
                    "object": "PDB:2N7Q",
                    "evidence_type": "ECO:0000250",
                    "supporting_objects": ["Ensembl:ENS1"],
                    "publications": ["PMID:1234"]
                }
            ]
        ),

        # comment catalytic activity (reaction) with evidence key
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <comment type="catalytic activity">
                     <reaction evidence="E2">
                        <dbReference type="Rhea" id="12345"/>
                     </reaction>
                   </comment>
            </entry>''',
            "CDM:3",
            {
                "E2": {
                    "evidence_type": "ECO:0000313",
                    "supporting_objects": None,
                    "publications": ["PMID:2222"]
                }
            },
            [
                {
                    "subject": "CDM:3",
                    "predicate": "catalyzes",
                    "object": "Rhea:12345",
                    "evidence_type": "ECO:0000313",
                    "supporting_objects": None,
                    "publications": ["PMID:2222"]
                }
            ]
        ),

        # Comment cofactor without evidence
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <comment type="cofactor">
                     <cofactor>
                       <dbReference type="ChEBI" id="CHEBI:15377"/>
                     </cofactor>
                   </comment>
            </entry>''',
            "CDM:4",
            {},
            [
                {
                    "subject": "CDM:4",
                    "predicate": "requires_cofactor",
                    "object": "ChEBI:CHEBI:15377"
                }
            ]
        ),

        # Several relevant relationship（with organism and dbReference）
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                   <organism>
                      <dbReference type="NCBI Taxonomy" id="562"/>
                   </organism>
                   <dbReference type="RefSeq" id="NP_414543"/>
            </entry>''',
            "CDM:5",
            {},
            [
                {'subject': "CDM:5", 'object': "NCBITaxon:562"},
                {'subject': "CDM:5", 'object': "RefSeq:NP_414543"}
            ]
        ),

        # if it is empty entry, return to []
        (
            '''<entry xmlns="https://uniprot.org/uniprot"></entry>''',
            "CDM:6",
            {},
            []
        ),
    ]
)
def test_parse_associations(xml_str, cdm_id, evidence_map, expected):
    entry = ET.fromstring(xml_str)
    result = parse_associations(entry, cdm_id, evidence_map)
    assert result == expected


from uniprot_to_cdm import parse_publications
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
@pytest.mark.parametrize(
    "xml_str, expected",
    [
        # Single PubMed 
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <reference>
                  <citation>
                    <dbReference type="PubMed" id="12345"/>
                  </citation>
                </reference>
            </entry>''',
            ["PMID:12345"]
        ),

        # 2. Multiple types include (PubMed, DOI, GenBank)
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <reference>
                  <citation>
                    <dbReference type="PubMed" id="55555"/>
                    <dbReference type="DOI" id="10.1000/j.jmb.2020.01.001"/>
                    <dbReference type="GenBank" id="AB123456"/>
                  </citation>
                </reference>
            </entry>''',
            [
                "PMID:55555",
                "DOI:10.1000/j.jmb.2020.01.001",
                "GenBank:AB123456"
            ]
        ),

        # Multiple references
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
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
            </entry>''',
            [
                "DOI:10.1000/jmb.123456",
                "PMID:98765"
            ]
        ),

        # dbReference: DDBJ and EMBL
        (
            '''<entry xmlns="https://uniprot.org/uniprot">
                <reference>
                  <citation>
                    <dbReference type="DDBJ" id="BA000001"/>
                    <dbReference type="EMBL" id="AB987654"/>
                  </citation>
                </reference>
            </entry>''',
            [
                "DDBJ:BA000001",
                "EMBL:AB987654"
            ]
        ),

        # no publication
        (
            '''<entry xmlns="https://uniprot.org/uniprot"></entry>''',
            []
        ),
    ]
)
def test_parse_publications(xml_str, expected):
    entry = ET.fromstring(xml_str)
    result = parse_publications(entry)
    assert result == expected



from uniprot_to_cdm import parse_uniprot_entry  

"""
    The function validates the correctness and integration of the parse_uniprot_entry function 
    by providing a minimal but representative UniProt XML entry.

    Ensure integration stability and correct results of your entire UniProt → CDM parsing process

"""

def test_parse_uniprot_entry_basic(monkeypatch):
    # Sets up a fixed UUID return value using monkeypatch
    monkeypatch.setattr("uniprot_to_cdm.generate_unique_id", lambda: "CDM:TESTUUID")

    # Defines a sample XML string with typical UniProt entry content
    # Accession, top-level name, recommended name, organism taxonomy, publication refer
    
    xml_str = '''
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
    '''

    # Parses the XML string into an Element object 
    entry = ET.fromstring(xml_str)

    # parse_uniprot_entry to extract all standardize CDM components from XML 
    record = parse_uniprot_entry(entry)
    
    # Assert the entity structures 
    assert record["entity"] == {
        "entity_id": "CDM:TESTUUID",
        "entity_type": "protein",
        "data_source": "UniProt archaea",
        "created": "2020-01-01",
        "updated": "2021-01-01",
        "version": "3"
    }

    # Only one identifiers 
    assert record["identifiers"] == [{
        "entity_id": "CDM:TESTUUID",
        "identifier": "UniProt:P12345",
        "source": "UniProt",
        "description": "UniProt accession"
    }]

    # Assert names strutures 
    assert {"entity_id": "CDM:TESTUUID", "name": "ProteinX", "description": "UniProt protein name", "source": "UniProt"} in record["names"]
    assert {"entity_id": "CDM:TESTUUID", "name": "ProteinX Full Name", "description": "UniProt recommended full name", "source": "UniProt"} in record["names"]
    # associations has organism
    assert {"subject": "CDM:TESTUUID", "object": "NCBITaxon:9606"} in record["associations"]
    # publications has PMID
    assert "PMID:99999" in record["publications"]

