"""Tests for the UniProtKB XML parser."""

import datetime
import re
from typing import Any

import pytest
from lxml.etree import XMLParser, fromstring

from cdm_data_loader_utils.parsers.uniprot.uniprot_kb import (
    CONTENT,
    DB,
    DESCRIPTION,
    ENTITY_ID,
    HTTPS_NS,
    KEY,
    NAME,
    NCBI_TAXON,
    XREF,
    dump_evidence_xml,
    parse_cross_references,
    parse_identifiers,
    parse_names,
    parse_organism,
    parse_protein_info,
    parse_references,
    parse_uniprot_entry,
)

TEST_ID = "uniprot:P01Q23"
FILE_PATH = "/path/to/file"


def entry_wrap(xml: str) -> str:
    """Wrap XML text in an `entry` element.

    :param xml: the XML to wrap
    :type xml: str
    :return: XML string with an entry around it.
    :rtype: str
    """
    return f"<entry xmlns='{HTTPS_NS}'>{xml}</entry>"


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        # Only <entry><name>
        (
            """
    <name>A4_HUMAN</name>""",
            [
                ("A4_HUMAN", "UniProt entry name"),
            ],
        ),
        # entry name + recommended full name
        (
            """
    <name>HumanA4</name>
    <protein>
      <recommendedName>
        <fullName evidence="75">Angiotensin-converting enzyme</fullName>
      </recommendedName>
    </protein>""",
            [
                ("HumanA4", "UniProt entry name"),
                ("Angiotensin-converting enzyme", "UniProt recommended full name"),
            ],
        ),
        # gene name
        (
            """
    <name>A4_HUMAN</name>
    <protein>
      <recommendedName />
    </protein>
    <gene>
      <name evidence="1 2" type="primary">atpB</name>
      <name evidence="3" type="synonym">ntpB</name>
      <name type="ordered locus">rrnAC3160</name>
    </gene>
""",
            [
                ("A4_HUMAN", "UniProt entry name"),
                ("atpB", "UniProt gene name, primary"),
                ("ntpB", "UniProt gene name, synonym"),
                ("rrnAC3160", "UniProt gene name, ordered locus"),
            ],
        ),
        # everything
        (
            """
    <name>A4_HUMAN</name>
    <name>HumanA4</name>
    <protein>
      <recommendedName>
        <fullName evidence="75">Angiotensin-converting enzyme</fullName>
        <shortName evidence="75">ACE</shortName>
        <ecNumber evidence="10 13 21 22 53 62">3.4.15.1</ecNumber>
      </recommendedName>
      <alternativeName>
        <fullName evidence="76">Dipeptidyl carboxypeptidase I</fullName>
      </alternativeName>
      <alternativeName>
        <fullName>DCP I</fullName>
      </alternativeName>
      <alternativeName>
        <fullName evidence="77">Kininase II</fullName>
        <shortName evidence="77">KNN II</shortName>
      </alternativeName>
    </protein>""",
            [
                ("A4_HUMAN", "UniProt entry name"),
                ("HumanA4", "UniProt entry name"),
                ("Angiotensin-converting enzyme", "UniProt recommended full name"),
                ("ACE", "UniProt recommended short name"),
                ("Dipeptidyl carboxypeptidase I", "UniProt alternative full name"),
                ("DCP I", "UniProt alternative full name"),
                ("Kininase II", "UniProt alternative full name"),
                ("KNN II", "UniProt alternative short name"),
            ],
        ),
        # submitted name
        (
            """
    <accession>P0C9F0</accession>
    <name>1001R_ASFK5</name>
    <protein>
      <submittedName>
        <fullName>Protein MGF 100-1R</fullName>
      </submittedName>
    </protein>""",
            [
                ("1001R_ASFK5", "UniProt entry name"),
                ("Protein MGF 100-1R", "UniProt submitted full name"),
            ],
        ),
        # domain and component names should be ignored
        (
            """
    <name>A4_HUMAN</name>
    <name>HumanA4</name>
    <protein>
      <cdAntigenName>CD143</cdAntigenName>
      <allergenName>Ulo b 1</allergenName>
      <allergenName>Alt a 1</allergenName>
      <innName evidence="13">Tozuleristide</innName>
      <innName>Glucarpidase</innName>
      <domain>
        <recommendedName>
          <fullName>Probable folylpolyglutamate synthase</fullName>
          <ecNumber>6.3.2.17</ecNumber>
        </recommendedName>
        <alternativeName>
          <fullName>Tetrahydrofolylpolyglutamate synthase</fullName>
          <shortName>Tetrahydrofolate synthase</shortName>
        </alternativeName>
      </domain>
      <component>
        <recommendedName>
          <fullName>C83</fullName>
        </recommendedName>
        <alternativeName>
          <fullName evidence="128">Alpha-secretase C-terminal fragment</fullName>
          <shortName evidence="128">Alpha-CTF</shortName>
        </alternativeName>
      </component>
    </protein>""",
            [
                ("A4_HUMAN", "UniProt entry name"),
                ("HumanA4", "UniProt entry name"),
                #   <cdAntigenName>CD143</cdAntigenName>
                #   <allergenName>Ulo b 1</allergenName>
                #   <allergenName>Alt a 1</allergenName>
                #   <innName evidence="13">Tozuleristide</innName>
                #   <innName>Glucarpidase</innName>
            ],
        ),
    ],
)
def test_parse_names(xml: str, expected: list[tuple[str, str]]) -> None:
    """Test the extraction of names from an entry."""
    xml_str = entry_wrap(xml)
    parsed_xrefs = parse_names(fromstring(xml_str))
    assert parsed_xrefs == [{NAME: e[0], DESCRIPTION: e[1]} for e in expected]


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        (
            # none of the useful fields present
            """""",
            {},
        ),
        (  # all fields present
            """
        <proteinExistence type="evidence at transcript level" />
        <sequence length="107" mass="12353" checksum="A0E89689D9BE1726" modified="2001-06-01" version="1">MVKAPRGYRNRTRRLLRKPVREKGSIPRLSTYLREYRVGDKVAIIINPSFPDWGMPHRRFHGLTGTVVGKRGEAYEVEVYLGRKRKTLFVPPVHLKPLSTAAERRGS</sequence>""",
            {
                "evidence_for_existence": "evidence at transcript level",
                "sequence": "MVKAPRGYRNRTRRLLRKPVREKGSIPRLSTYLREYRVGDKVAIIINPSFPDWGMPHRRFHGLTGTVVGKRGEAYEVEVYLGRKRKTLFVPPVHLKPLSTAAERRGS",
                "length": 107,
                "hash": "A0E89689D9BE1726",
            },
        ),
        (  # evidence for existence only
            """
        <proteinExistence type="evidence at transcript level" />
            """,
            {
                "evidence_for_existence": "evidence at transcript level",
            },
        ),
        (  # no sequence
            """
        <sequence length="107" mass="12353" checksum="A0E89689D9BE1726" modified="2001-06-01" version="1" />
            """,
            {
                "length": 107,
                "hash": "A0E89689D9BE1726",
            },
        ),
    ],
)
def test_parse_protein_info(xml: str, expected: dict[str, Any]) -> None:
    """Test extraction of core protein properties from the entry XML."""
    xml_str = entry_wrap(xml)
    result = parse_protein_info(fromstring(xml_str))
    assert result == expected


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        # No accession
        ("", []),
        # Single accession
        (
            """
    <accession>P12345</accession>""",
            [
                {
                    DB: "UniProt",
                    XREF: "P12345",
                    DESCRIPTION: "UniProt accession",
                }
            ],
        ),
        # Multiple accessions
        (
            """
    <accession>Q11111</accession>
    <accession>Q22222</accession>""",
            [
                {
                    DB: "UniProt",
                    XREF: "Q11111",
                    DESCRIPTION: "UniProt accession",
                },
                {
                    DB: "UniProt",
                    XREF: "Q22222",
                    DESCRIPTION: "UniProt accession",
                },
            ],
        ),
        # Empty accessions
        (
            """
    <accession>Q11111</accession>
    <accession></accession>
    <accession></accession>
    <accession>Q22222</accession>""",
            [
                {
                    DB: "UniProt",
                    XREF: "Q11111",
                    DESCRIPTION: "UniProt accession",
                },
                {
                    DB: "UniProt",
                    XREF: "Q22222",
                    DESCRIPTION: "UniProt accession",
                },
            ],
        ),
    ],
)
def test_parse_identifiers(xml: str, expected: list[dict[str, Any]]) -> None:
    """Test the parsing of accessions."""
    xml_str = entry_wrap(xml)
    result = parse_identifiers(fromstring(xml_str))
    assert result == expected


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        (  # simple evidence
            """
    <evidence type="ECO:0000255" key="4"/>
    <evidence type="ECO:0000305" key="7"/>""",
            [
                {KEY: "4", CONTENT: '<evidence type="ECO:0000255" key="4"/>'},
                {KEY: "7", CONTENT: '<evidence type="ECO:0000305" key="7"/>'},
            ],
        ),
        (  # evidence with source
            """
    <evidence type="ECO:0000250" key="1">
      <source>
        <dbReference type="UniProtKB" id="B9KDD4"/>
      </source>
    </evidence>
    <evidence type="ECO:0000250" key="2">
      <source>
        <dbReference type="UniProtKB" id="O29867"/>
      </source>
    </evidence>""",
            [
                {
                    KEY: "1",
                    CONTENT: '<evidence type="ECO:0000250" key="1"><source><dbReference type="UniProtKB" id="B9KDD4"/></source></evidence>',
                },
                {
                    KEY: "2",
                    CONTENT: '<evidence type="ECO:0000250" key="2"><source><dbReference type="UniProtKB" id="O29867"/></source></evidence>',
                },
            ],
        ),
        (  # evidence with ref
            """
    <evidence type="ECO:0000269" key="8">
      <source ref="3"/>
    </evidence>
    <evidence type="ECO:0000303" key="9">
      <source ref="3"/>
    </evidence>""",
            [
                {KEY: "8", CONTENT: '<evidence type="ECO:0000269" key="8"><source ref="3"/></evidence>'},
                {
                    KEY: "9",
                    CONTENT: '<evidence type="ECO:0000303" key="9"><source ref="3"/></evidence>',
                },
            ],
        ),
        (  # evidence imported from
            """
    <evidence type="string" key="10">
      <source ref="2" />
      <importedFrom>
        <dbReference type="string" id="string" />
      </importedFrom>
    </evidence>""",
            [
                {
                    KEY: "10",
                    CONTENT: '<evidence type="string" key="10"><source ref="2"/><importedFrom><dbReference type="string" id="string"/></importedFrom></evidence>',
                }
            ],
        ),
    ],
)
def test_parse_evidence(xml: str, expected: dict[str, str]) -> None:
    """Test the parsing of evidence data."""
    xml_str = entry_wrap(xml)
    parsed_evidence = dump_evidence_xml(fromstring(xml_str, parser=XMLParser(remove_blank_text=True)), TEST_ID)
    assert parsed_evidence == [{ENTITY_ID: TEST_ID, **e} for e in expected]


REF_TEXT = {
    1: """
    <reference evidence="26" key="1">
      <citation type="journal article" date="1991" name="Gene" volume="102" first="117" last="122">
        <title>Structure and organization of the gas vesicle gene cluster on the Halobacterium halobium plasmid pNRC100.</title>
        <authorList>
          <person name="Jones J.G."/>
          <person name="Young D.C."/>
          <person name="Dassarma S."/>
        </authorList>
        <dbReference type="PubMed" id="1864501"/>
        <dbReference type="PMCID" id="22222222"/>
        <dbReference type="DOI" id="10.1016/0378-1119(91)90549-q"/>
      </citation>
      <scope>NUCLEOTIDE SEQUENCE [GENOMIC DNA]</scope>
      <source>
        <strain>ATCC 700922 / JCM 11081 / NRC-1</strain>
        <plasmid>pNRC100</plasmid>
      </source>
    </reference>""",
    2: """
    <reference key="2">
      <citation type="patent" date="1999-03-11" number="WO9911821">
        <title>Method for screening restriction endonucleases.</title>
        <authorList>
          <person name="Noren C.J."/>
          <person name="Roberts R.J."/>
          <person name="Patti J."/>
          <person name="Byrd D.R."/>
          <person name="Morgan R.D."/>
        </authorList>
      </citation>
      <scope>CHARACTERIZATION</scope>
    </reference>""",
    3: """
    <reference key="3">
      <citation type="book" date="1980" name="Genetics and evolution of RNA polymerase, tRNA and ribosomes" first="585" last="599" publisher="University of Tokyo Press" city="Tokyo">
        <editorList>
          <person name="Osawa S."/>
          <person name="Ozeki H."/>
          <person name="Uchida H."/>
          <person name="Yura T."/>
        </editorList>
        <authorList>
          <person name="Yaguchi M."/>
          <person name="Matheson A.T."/>
          <person name="Visentin L.P."/>
          <person name="Zuker M."/>
        </authorList>
      </citation>
      <scope>PROTEIN SEQUENCE OF 1-76</scope>
    </reference>""",
    4: """
    <reference key="4">
      <citation type="unpublished observations" date="2001-05">
        <authorList>
          <person name="Medigue C."/>
          <person name="Bocs S."/>
        </authorList>
        <dbReference type="invalid_ref_type" id="whatever"/>
      </citation>
      <scope>IDENTIFICATION</scope>
    </reference>""",
    5: """
    <reference key="5">
      <citation type="submission" date="2005-03" db="EMBL/GenBank/DDBJ databases">
        <authorList>
          <person name="Totoki Y."/>
          <person name="Toyoda A."/>
          <person name="Takeda T."/>
          <person name="Sakaki Y."/>
          <person name="Tanaka A."/>
          <person name="Yokoyama S."/>
          <person name="Ohara O."/>
          <person name="Nagase T."/>
          <person name="Kikuno R.F."/>
        </authorList>
        <dbReference type="AGRICOLA" id="12345678"/>
      </citation>
      <scope>NUCLEOTIDE SEQUENCE [LARGE SCALE MRNA] OF 1-1239 (ISOFORM SOMATIC-1)</scope>
      <source>
        <tissue>Brain</tissue>
      </source>
    </reference>""",
    6: """
    <reference key="6">
        <citation type="journal article" date="2008" name="Appl. Environ. Microbiol." volume="74" first="682" last="692">
            <title>The genome sequence of the metal-mobilizing, extremely thermoacidophilic archaeon Metallosphaera sedula provides insights into bioleaching-associated metabolism.</title>
            <authorList>
                <person name="Auernik K.S."/>
                <person name="Maezato Y."/>
                <person name="Blum P.H."/>
                <person name="Kelly R.M."/>
            </authorList>
            <dbReference type="PubMed" id="18083856"/>
        </citation>
        <scope>NUCLEOTIDE SEQUENCE [LARGE SCALE GENOMIC DNA]</scope>
        <source>
            <strain>ATCC 51363 / DSM 5348 / JCM 9185 / NBRC 15509 / TH2</strain>
        </source>
    </reference>
    """,
}

PARSED_REFS = {
    1: {"publication_id": "DOI:10.1016/0378-1119(91)90549-q"},
    5: {"publication_id": "AGRICOLA:12345678"},
    6: {"publication_id": "PMID:18083856"},
}


@pytest.mark.parametrize(
    ("xml", "expected", "key"),
    [(REF_TEXT[key], [PARSED_REFS.get(key)] if PARSED_REFS.get(key) else [], str(key)) for key in REF_TEXT],
)
def test_parse_references(xml: str, expected: list[dict[str, str]] | None, key: str) -> None:
    """Test the parsing of reference data."""
    xml_str = entry_wrap(xml)
    parsed_refs = parse_references(fromstring(xml_str, parser=XMLParser(remove_blank_text=True)))
    assert parsed_refs["publication"] == expected
    no_ws_xml = re.sub(r">\s+<", "><", xml.strip())
    assert parsed_refs["all_xml"] == [{KEY: key, CONTENT: no_ws_xml}]


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        # simple db / xref
        ("""<dbReference type="PDBsum" id="1O86"/>""", [{DB: "PDBsum", XREF: "1O86"}]),
        ("""<dbReference type="AlphaFoldDB" id="P12821"/>""", [{DB: "AlphaFoldDB", XREF: "P12821"}]),
        ("""<dbReference type="EMDB" id="EMD-13797"/>""", [{DB: "EMDB", XREF: "EMD-13797"}]),
        ("""<dbReference type="SMR" id="P12821"/>""", [{DB: "SMR", XREF: "P12821"}]),
        ("""<dbReference type="CORUM" id="P12821"/>""", [{DB: "CORUM", XREF: "P12821"}]),
        ("""<dbReference type="MINT" id="P12821"/>""", [{DB: "MINT", XREF: "P12821"}]),
        (
            """<dbReference type="STRING" id="9606.ENSP00000290866"/>""",
            [{DB: "STRING", XREF: "9606.ENSP00000290866"}],
        ),
        ("""<dbReference type="BindingDB" id="P12821"/>""", [{DB: "BindingDB", XREF: "P12821"}]),
        ("""<dbReference type="ChEMBL" id="CHEMBL1808"/>""", [{DB: "ChEMBL", XREF: "CHEMBL1808"}]),
        ("""<dbReference type="DrugCentral" id="P12821"/>""", [{DB: "DrugCentral", XREF: "P12821"}]),
        ("""<dbReference type="GuidetoPHARMACOLOGY" id="1613"/>""", [{DB: "GuidetoPHARMACOLOGY", XREF: "1613"}]),
        ("""<dbReference type="MEROPS" id="M02.001"/>""", [{DB: "MEROPS", XREF: "M02.001"}]),
        ("""<dbReference type="DNASU" id="1636"/>""", [{DB: "DNASU", XREF: "1636"}]),
        ("""<dbReference type="GeneID" id="1636"/>""", [{DB: "GeneID", XREF: "1636"}]),
        ("""<dbReference type="KEGG" id="hsa:1636"/>""", [{DB: "KEGG", XREF: "hsa:1636"}]),
        # with evidence
        ("""<dbReference type="EC" id="3.4.15.1"/>""", [{DB: "EC", XREF: "3.4.15.1"}]),
        # ref containing a colon
        (
            """<dbReference type="FunFam" id="2.30.30.70:FF:000001">
      <property type="entry name" value="60S ribosomal protein L21"/>
      <property type="match status" value="1"/>
    </dbReference>""",
            [{DB: "FunFam", XREF: "2.30.30.70:FF:000001"}],
        ),
        (  # GO refs
            """<dbReference type="GO" id="GO:0005840">
      <property type="term" value="C:ribosome"/>
      <property type="evidence" value="ECO:0007669"/>
      <property type="project" value="UniProtKB-KW"/>
    </dbReference>
    <dbReference type="GO" id="GO:0003735">
      <property type="term" value="F:structural constituent of ribosome"/>
      <property type="evidence" value="ECO:0007669"/>
      <property type="project" value="InterPro"/>
    </dbReference>
    <dbReference type="GO" id="GO:0006412">
      <property type="term" value="P:translation"/>
      <property type="evidence" value="ECO:0007669"/>
      <property type="project" value="UniProtKB-UniRule"/>
    </dbReference>""",
            [
                {DB: "GO", XREF: "0005840"},
                {DB: "GO", XREF: "0003735"},
                {DB: "GO", XREF: "0006412"},
            ],
        ),
        # with molecule ID
        (
            """<dbReference type="CCDS" id="CCDS11637.1">
      <molecule id="P12821-1"/>
    </dbReference>""",
            [{DB: "CCDS", XREF: "CCDS11637.1", DESCRIPTION: "CCDS ID for UniProt:P12821-1"}],
        ),
        # property, not saved
        (
            """<dbReference type="BioGRID" id="108004">
      <property type="interactions" value="17"/>
    </dbReference>""",
            [{DB: "BioGRID", XREF: "108004"}],
        ),
        (
            """<dbReference type="FunCoup" id="P12821">
      <property type="interactions" value="323"/>
    </dbReference>""",
            [{DB: "FunCoup", XREF: "P12821"}],
        ),
        (
            """<dbReference type="IntAct" id="P12821">
      <property type="interactions" value="8"/>
    </dbReference>""",
            [{DB: "IntAct", XREF: "P12821"}],
        ),
        (
            """<dbReference type="DrugBank" id="DB14511">
      <property type="generic name" value="Acetate"/>
    </dbReference>""",
            [{DB: "DrugBank", XREF: "DB14511"}],
        ),
        (
            """<dbReference type="GlyConnect" id="1011">
      <property type="glycosylation" value="24 N-Linked glycans (7 sites)"/>
    </dbReference>""",
            [{DB: "GlyConnect", XREF: "1011"}],
        ),
        (
            """<dbReference type="Antibodypedia" id="31288">
      <property type="antibodies" value="911 antibodies from 43 providers"/>
    </dbReference>""",
            [{DB: "Antibodypedia", XREF: "31288"}],
        ),
        # with property
        (
            """<dbReference type="PIR" id="A31759">
      <property type="entry name" value="A31759"/>
    </dbReference>""",
            [{DB: "PIR", XREF: "A31759"}],
        ),
        # with properties, EMBL style
        (
            """<dbReference type="EMBL" id="J04144">
      <property type="protein sequence ID" value="AAA51684.1"/>
      <property type="molecule type" value="mRNA"/>
    </dbReference>""",
            [
                {DB: "genbank", XREF: "J04144", DESCRIPTION: "EMBL/GenBank mRNA ID"},
                {DB: "genbank", XREF: "AAA51684.1", DESCRIPTION: "EMBL/GenBank protein sequence ID"},
            ],
        ),
        # with properties, EMBL style mk II
        (
            """    <dbReference type="EMBL" id="AY261360">
      <property type="status" value="NOT_ANNOTATED_CDS"/>
      <property type="molecule type" value="Genomic_DNA"/>
    </dbReference>""",
            [{DB: "genbank", XREF: "AY261360", DESCRIPTION: "EMBL/GenBank Genomic_DNA ID"}],
        ),
        # multiple properties, but we don't care
        (
            """<dbReference type="PDB" id="1O86">
      <property type="method" value="X-ray"/>
      <property type="resolution" value="2.00 A"/>
      <property type="chains" value="A=642-1230"/>
    </dbReference>""",
            [{DB: "PDB", XREF: "1O86"}],
        ),
        # MANE-Select -- genbank and refseq entries
        (
            """<dbReference type="MANE-Select" id="ENST00000290866.10">
      <property type="protein sequence ID" value="ENSP00000290866.4"/>
      <property type="RefSeq nucleotide sequence ID" value="NM_000789.4"/>
      <property type="RefSeq protein sequence ID" value="NP_000780.1"/>
    </dbReference>""",
            [
                {
                    DB: "ensembl",
                    XREF: "ENST00000290866.10",
                    DESCRIPTION: "Ensembl transcript ID, via MANE-Select",
                },
                {
                    DB: "ensembl",
                    XREF: "ENSP00000290866.4",
                    DESCRIPTION: "Ensembl protein sequence ID, via MANE-Select",
                },
                {
                    DB: "refseq",
                    XREF: "NM_000789.4",
                    DESCRIPTION: "RefSeq nucleotide sequence ID, via MANE-Select",
                },
                {DB: "refseq", XREF: "NP_000780.1", DESCRIPTION: "RefSeq protein sequence ID, via MANE-Select"},
            ],
        ),
        # with molecule and property
        (
            """<dbReference type="RefSeq" id="NP_000780.1">
      <molecule id="P12821-1"/>
      <property type="nucleotide sequence ID" value="NM_000789.4"/>
    </dbReference>""",
            [
                {
                    DB: "refseq",
                    XREF: "NP_000780.1",
                    DESCRIPTION: "RefSeq protein sequence ID for UniProt:P12821-1",
                },
                {
                    DB: "refseq",
                    XREF: "NM_000789.4",
                    DESCRIPTION: "RefSeq nucleotide sequence ID for UniProt:P12821-1",
                },
            ],
        ),
        (
            """<dbReference type="Ensembl" id="ENST00000290863.10">
      <molecule id="P12821-3"/>
      <property type="protein sequence ID" value="ENSP00000290863.6"/>
      <property type="gene ID" value="ENSG00000159640.17"/>
    </dbReference>""",
            [
                {
                    DB: "ensembl",
                    XREF: "ENST00000290863.10",
                    DESCRIPTION: "Ensembl transcript ID for UniProt:P12821-3",
                },
                {
                    DB: "ensembl",
                    XREF: "ENSP00000290863.6",
                    DESCRIPTION: "Ensembl protein sequence ID for UniProt:P12821-3",
                },
                {DB: "ensembl", XREF: "ENSG00000159640.17", DESCRIPTION: "Ensembl gene ID for UniProt:P12821-3"},
            ],
        ),
        (
            """<dbReference type="Ensembl" id="ENST00000413513.7">
      <molecule id="P12821-4"/>
      <property type="protein sequence ID" value="ENSP00000392247.3"/>
      <property type="gene ID" value="ENSG00000159640.17"/>
    </dbReference>""",
            [
                {
                    DB: "ensembl",
                    XREF: "ENST00000413513.7",
                    DESCRIPTION: "Ensembl transcript ID for UniProt:P12821-4",
                },
                {
                    DB: "ensembl",
                    XREF: "ENSP00000392247.3",
                    DESCRIPTION: "Ensembl protein sequence ID for UniProt:P12821-4",
                },
                {DB: "ensembl", XREF: "ENSG00000159640.17", DESCRIPTION: "Ensembl gene ID for UniProt:P12821-4"},
            ],
        ),
        # keep molecule, discard property
        (
            """<dbReference type="UCSC" id="uc002jau.3">
      <molecule id="P12821-1"/>
      <property type="organism name" value="human"/>
    </dbReference>""",
            [{DB: "UCSC", XREF: "uc002jau.3", DESCRIPTION: "UCSC ID for UniProt:P12821-1"}],
        ),
        (  # invalid refs
            """<dbReference type="GO">
      <property type="term" value="C:ribosome"/>
      <property type="evidence" value="ECO:0007669"/>
      <property type="project" value="UniProtKB-KW"/>
    </dbReference>
    <dbReference id="GO:0003735">
      <property type="project" value="InterPro"/>
    </dbReference>
    <dbReference>
      <property type="project" value="UniProtKB-UniRule"/>
    </dbReference>
    <dbReference />""",
            [],
        ),
    ],
)
def test_parse_cross_references(xml: str, expected: list[dict[str, str]]) -> None:
    """Test the parsing of generic cross references."""
    xml_str = entry_wrap(xml)
    parsed_dbxrefs = parse_cross_references(fromstring(xml_str))
    assert parsed_dbxrefs == expected


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        (
            """
    <organism>
        <name type="scientific">Pyrococcus furiosus (strain ATCC 43587 / DSM 3638 / JCM 8422 / Vc1)</name>
        <dbReference type="NCBI Taxonomy" id="186497"/>
        <lineage>
            <taxon>Archaea</taxon>
            <taxon>Methanobacteriati</taxon>
            <taxon>Methanobacteriota</taxon>
            <taxon>Thermococci</taxon>
            <taxon>Thermococcales</taxon>
            <taxon>Thermococcaceae</taxon>
            <taxon>Pyrococcus</taxon>
        </lineage>
    </organism>
        """,
            [
                {
                    DB: NCBI_TAXON,
                    XREF: "186497",
                    DESCRIPTION: "UniProt taxon designation",
                    "relationship": "RO:0002162: in taxon",
                }
            ],
        ),
        (
            """
    <organism>
        <name type="scientific">African swine fever virus (isolate Pig/Kenya/KEN-50/1950)</name>
        <name type="common">ASFV</name>
        <dbReference type="NCBI Taxonomy" id="561445"/>
        <lineage>
            <taxon>Viruses</taxon>
            <taxon>Varidnaviria</taxon>
            <taxon>Bamfordvirae</taxon>
            <taxon>Nucleocytoviricota</taxon>
            <taxon>Pokkesviricetes</taxon>
            <taxon>Asfuvirales</taxon>
            <taxon>Asfarviridae</taxon>
            <taxon>Asfivirus</taxon>
            <taxon>African swine fever virus</taxon>
        </lineage>
    </organism>
    <organismHost>
        <name type="scientific">Ornithodoros</name>
        <name type="common">relapsing fever ticks</name>
        <dbReference type="NCBI Taxonomy" id="6937"/>
    </organismHost>
    <organismHost>
        <name type="scientific">Phacochoerus aethiopicus</name>
        <name type="common">Warthog</name>
        <dbReference type="NCBI Taxonomy" id="85517"/>
    </organismHost>
    <organismHost>
        <name type="scientific">Phacochoerus africanus</name>
        <name type="common">Warthog</name>
        <dbReference type="NCBI Taxonomy" id="41426"/>
    </organismHost>
    <organismHost>
        <name type="scientific">Potamochoerus larvatus</name>
        <name type="common">Bushpig</name>
        <dbReference type="NCBI Taxonomy" id="273792"/>
    </organismHost>
    <organismHost>
        <name type="scientific">Sus scrofa</name>
        <name type="common">Pig</name>
        <dbReference type="NCBI Taxonomy" id="9823"/>
    </organismHost>
    """,
            [
                {
                    DB: NCBI_TAXON,
                    XREF: "561445",
                    DESCRIPTION: "UniProt taxon designation",
                    "relationship": "RO:0002162: in taxon",
                },
            ],
        ),
        # no organism section -- would not happen
        ("""""", []),
        # no taxonomy ID -- would not happen
        (
            """
    <organism>
        <name type="scientific">African swine fever virus (isolate Pig/Kenya/KEN-50/1950)</name>
        <name type="common">ASFV</name>
        <lineage>
            <taxon>Viruses</taxon>
            <taxon>Varidnaviria</taxon>
            <taxon>Bamfordvirae</taxon>
            <taxon>Nucleocytoviricota</taxon>
            <taxon>Pokkesviricetes</taxon>
            <taxon>Asfuvirales</taxon>
            <taxon>Asfarviridae</taxon>
            <taxon>Asfivirus</taxon>
            <taxon>African swine fever virus</taxon>
        </lineage>
    </organism>
""",
            [],
        ),
    ],
)
def test_parse_organism(xml: str, expected: list[dict[str, str]]) -> None:
    """Test the parsing of organism taxon IDs."""
    xml_str = entry_wrap(xml)
    parsed_dbxrefs = parse_organism(fromstring(xml_str))
    assert parsed_dbxrefs == expected


@pytest.mark.parametrize(
    ("xml", "expected"),
    [
        ("", {"_parse_error": [{"xml": "<entry/>", "error": "No accession found", "source_file": FILE_PATH}]}),
        (
            "<dbReference/><name/>",
            {
                "_parse_error": [
                    {
                        "xml": "<entry><dbReference/><name/></entry>",
                        "error": "No accession found",
                        "source_file": FILE_PATH,
                    }
                ]
            },
        ),
        # TODO: put together a full example!
    ],
)
def test_parse_uniprot_entry(xml: str, expected: dict[str, list[dict[str, Any]]]) -> None:
    """Test the parsing of generic cross references."""
    xml_str = entry_wrap(xml)
    parsed_entry = parse_uniprot_entry(
        fromstring(xml_str, parser=XMLParser(remove_blank_text=True)),
        datetime.datetime.now(tz=datetime.UTC),
        FILE_PATH,
    )
    assert parsed_entry == expected
