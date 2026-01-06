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
   pytest tests/uniprot_refactor/test_uniprot_parsers.py

"""

import datetime
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from cdm_data_loader_utils.parsers.uniprot import (
    build_datasource_record,
    parse_associations,
    parse_cross_references,
    parse_evidence_map,
    parse_identifiers,
    parse_names,
    parse_protein_info,
    save_datasource_record,
)

NS_URI = "https://uniprot.org/uniprot"


@pytest.fixture(
    params=[
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz",
        "http://example.org/uniprot_test.xml.gz",
    ]
)
def xml_url(request):
    return request.param


def test_build_datasource_record(xml_url):
    record = build_datasource_record(xml_url)

    # ---- basic structure ----
    assert isinstance(record, dict)

    # ---- fixed fields ----
    assert record["name"] == "UniProt import"
    assert record["source"] == "UniProt"
    assert record["url"] == xml_url
    assert record["version"] == 115

    # ---- accessed field ----
    accessed = record.get("accessed")
    assert accessed is not None

    parsed = datetime.datetime.fromisoformat(accessed)
    assert parsed.tzinfo is not None
    assert parsed.tzinfo == datetime.UTC


def test_save_datasource_record(tmp_path: Path, xml_url):
    """
    save_datasource_record should:
    - create output directory if missing
    - write datasource.json
    - return the same content that is written to disk
    """
    output_dir = tmp_path / "output"

    # ---- call function ----
    result = save_datasource_record(xml_url, str(output_dir))

    # ---- return value sanity ----
    assert isinstance(result, dict)
    assert result["url"] == xml_url
    assert result["source"] == "UniProt"
    assert result["name"] == "UniProt import"
    assert "accessed" in result
    assert "version" in result

    # ---- file existence ----
    output_file = output_dir / "datasource.json"
    assert output_file.exists()
    assert output_file.is_file()

    # ---- file content correctness ----
    with open(output_file, encoding="utf-8") as f:
        on_disk = json.load(f)

    assert on_disk == result


def make_entry(names=None, protein_names=None):
    entry = ET.Element(f"{{{NS_URI}}}entry")

    # <entry><name>
    for n in names or []:
        e = ET.SubElement(entry, f"{{{NS_URI}}}name")
        e.text = n

    # <protein> block
    if protein_names:
        protein = ET.SubElement(entry, f"{{{NS_URI}}}protein")

        for tag, logical in [
            ("recommendedName", "recommended"),
            ("alternativeName", "alternative"),
        ]:
            if logical not in protein_names:
                continue

            block = ET.SubElement(protein, f"{{{NS_URI}}}{tag}")
            for xml_tag in ["fullName", "shortName"]:
                val = protein_names[logical].get(xml_tag.replace("Name", ""))
                if val:
                    e = ET.SubElement(block, f"{{{NS_URI}}}{xml_tag}")
                    e.text = val

    return entry


@pytest.mark.parametrize(
    "entry_kwargs, cdm_id, expected",
    [
        # Only <entry><name>
        (
            {"names": ["ProteinA"]},
            "cdm_1",
            {
                ("ProteinA", "UniProt entry name"),
            },
        ),
        # entry name + recommended full name
        (
            {
                "names": ["ProteinB"],
                "protein_names": {
                    "recommended": {"full": "Rec Full B", "short": None},
                },
            },
            "cdm_2",
            {
                ("ProteinB", "UniProt entry name"),
                ("Rec Full B", "UniProt recommended full name"),
            },
        ),
        # everything
        (
            {
                "names": ["ProteinC"],
                "protein_names": {
                    "recommended": {"full": "Rec Full C", "short": "Rec Short C"},
                    "alternative": {"full": "Alt Full C", "short": "Alt Short C"},
                },
            },
            "cdm_3",
            {
                ("ProteinC", "UniProt entry name"),
                ("Rec Full C", "UniProt recommended full name"),
                ("Rec Short C", "UniProt recommended short name"),
                ("Alt Full C", "UniProt alternative full name"),
                ("Alt Short C", "UniProt alternative short name"),
            },
        ),
    ],
)
def test_parse_names_parametrized(entry_kwargs, cdm_id, expected):
    entry = make_entry(**entry_kwargs)

    rows = parse_names(entry, cdm_id)

    # ---- row count ----
    assert len(rows) == len(expected)

    # ---- content ----
    observed = {(r["name"], r["description"]) for r in rows}
    assert observed == expected

    # ---- entity_id and source ----
    for r in rows:
        assert r["entity_id"] == cdm_id
        assert r["source"] == "UniProt"


@pytest.mark.parametrize(
    "build_entry, cdm_id, expected",
    [
        # --------------------------------------------------
        # Empty entry -> None
        # --------------------------------------------------
        (
            lambda: ET.Element(f"{{{NS_URI}}}entry"),
            "cdm_1",
            None,
        ),
        # --------------------------------------------------
        # Only EC numbers
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(
                        ET.SubElement(
                            ET.SubElement(entry, f"{{{NS_URI}}}protein"),
                            f"{{{NS_URI}}}recommendedName",
                        ),
                        f"{{{NS_URI}}}ecNumber",
                    ).__setattr__("text", "1.1.1.1"),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_2",
            {
                "ec_numbers": "1.1.1.1",
            },
        ),
        # --------------------------------------------------
        # Only sequence + entry modified
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    entry.set("modified", "2024-01-01"),
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}sequence",
                        {
                            "length": "100",
                            "mass": "12345",
                            "checksum": "ABC",
                            "version": "2",
                        },
                    ).__setattr__("text", "MKTIIALSY"),
                    entry,
                )[2]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_3",
            {
                "length": "100",
                "mass": "12345",
                "checksum": "ABC",
                "sequence_version": "2",
                "sequence": "MKTIIALSY",
                "entry_modified": "2024-01-01",
            },
        ),
        # --------------------------------------------------
        # Everything
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    entry.set("modified", "2024-02-02"),
                    # protein + EC
                    ET.SubElement(
                        ET.SubElement(
                            ET.SubElement(entry, f"{{{NS_URI}}}protein"),
                            f"{{{NS_URI}}}recommendedName",
                        ),
                        f"{{{NS_URI}}}ecNumber",
                    ).__setattr__("text", "3.5.4.4"),
                    # proteinExistence
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}proteinExistence",
                        {"type": "evidence at protein level"},
                    ),
                    # sequence
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}sequence",
                        {
                            "length": "250",
                            "mass": "99999",
                            "checksum": "XYZ",
                            "modified": "2023-12-01",
                            "version": "1",
                        },
                    ).__setattr__("text", "MADEUPSEQUENCE"),
                    entry,
                )[4]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_4",
            {
                "ec_numbers": "3.5.4.4",
                "protein_id": "cdm_4",
                "evidence_for_existence": "evidence at protein level",
                "length": "250",
                "mass": "99999",
                "checksum": "XYZ",
                "modified": "2023-12-01",
                "sequence_version": "1",
                "sequence": "MADEUPSEQUENCE",
                "entry_modified": "2024-02-02",
            },
        ),
    ],
)
def test_parse_protein_info(build_entry, cdm_id, expected):
    entry = build_entry()

    result = parse_protein_info(entry, cdm_id)

    if expected is None:
        assert result is None
    else:
        assert isinstance(result, dict)
        assert result == expected


@pytest.mark.parametrize(
    "build_xml, expected",
    [
        # --------------------------------------------------
        # No evidence elements
        # --------------------------------------------------
        (
            lambda: ET.Element(f"{{{NS_URI}}}entry"),
            {},
        ),
        # --------------------------------------------------
        # Evidence without key
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(entry, f"{{{NS_URI}}}evidence", {"type": "ECO:0000269"}),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            {},
        ),
        # --------------------------------------------------
        # Evidence with key, no source
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}evidence",
                        {"key": "1", "type": "ECO:0000313"},
                    ),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            {
                "1": {
                    "evidence_type": "ECO:0000313",
                }
            },
        ),
        # --------------------------------------------------
        # Evidence with PUBMED with other refs
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    lambda ev: (
                        ET.SubElement(
                            ET.SubElement(ev, f"{{{NS_URI}}}source"),
                            f"{{{NS_URI}}}dbReference",
                            {"type": "PubMed", "id": "12345"},
                        ),
                        ET.SubElement(
                            ET.SubElement(ev, f"{{{NS_URI}}}source"),
                            f"{{{NS_URI}}}dbReference",
                            {"type": "GO", "id": "GO:0008150"},
                        ),
                        entry,
                    )[2]
                )(
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}evidence",
                        {"key": "E2", "type": "ECO:0000269"},
                    )
                )
            )(ET.Element(f"{{{NS_URI}}}entry")),
            {
                "E2": {
                    "evidence_type": "ECO:0000269",
                    "publications": ["PMID:12345"],
                }
            },
        ),
    ],
)
def test_parse_evidence_map_parametrized(build_xml, expected):
    entry = build_xml()
    result = parse_evidence_map(entry)

    assert isinstance(result, dict)
    assert result == expected


@pytest.mark.parametrize(
    "build_xml, cdm_id, evidence_map, expected",
    [
        # --------------------------------------------------
        # Taxonomy association only
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(
                        ET.SubElement(entry, f"{{{NS_URI}}}organism"),
                        f"{{{NS_URI}}}dbReference",
                        {"type": "NCBI Taxonomy", "id": "1234"},
                    ),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_1",
            {},
            [
                {
                    "subject": "cdm_1",
                    "object": "NCBITaxon:1234",
                    "predicate": "in_taxon",
                }
            ],
        ),
        # --------------------------------------------------
        # Catalytic activity with evidence
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    lambda comment: (
                        lambda reaction: (
                            ET.SubElement(
                                reaction,
                                f"{{{NS_URI}}}dbReference",
                                {"type": "Rhea", "id": "RHEA:12345"},
                            ),
                            entry,
                        )[1]
                    )(
                        ET.SubElement(
                            comment,
                            f"{{{NS_URI}}}reaction",
                            {"evidence": "E1"},
                        )
                    )
                )(
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}comment",
                        {"type": "catalytic activity"},
                    )
                )
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_2",
            {
                "E1": {
                    "evidence_type": "ECO:0000269",
                    "publications": ["PMID:12345"],
                }
            },
            [
                {
                    "subject": "cdm_2",
                    "predicate": "catalyzes",
                    "object": "Rhea:RHEA:12345",
                    "evidence_type": "ECO:0000269",
                    "publications": ["PMID:12345"],
                }
            ],
        ),
        # --------------------------------------------------
        # Cofactor association
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    lambda comment: (
                        ET.SubElement(
                            ET.SubElement(
                                comment,
                                f"{{{NS_URI}}}cofactor",
                            ),
                            f"{{{NS_URI}}}dbReference",
                            {"type": "ChEBI", "id": "CHEBI:15377"},
                        ),
                        entry,
                    )[1]
                )(
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}comment",
                        {"type": "cofactor"},
                    )
                )
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_3",
            {},
            [
                {
                    "subject": "cdm_3",
                    "predicate": "requires_cofactor",
                    "object": "ChEBI:CHEBI:15377",
                }
            ],
        ),
    ],
)
def test_parse_associations_parametrized(build_xml, cdm_id, evidence_map, expected):
    entry = build_xml()

    result = parse_associations(entry, cdm_id, evidence_map)

    assert isinstance(result, list)
    assert result == expected


@pytest.mark.parametrize(
    "build_xml, cdm_id, expected",
    [
        # --------------------------------------------------
        # No dbReference
        # --------------------------------------------------
        (
            lambda: ET.Element(f"{{{NS_URI}}}entry"),
            "cdm_1",
            [],
        ),
        # --------------------------------------------------
        # dbReference with CURIE id
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}dbReference",
                        {"type": "GO", "id": "GO:0008150"},
                    ),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_2",
            [
                {
                    "entity_id": "cdm_2",
                    "xref_type": "GO",
                    "xref_value": "GO:0008150",
                    "xref": "GO:0008150",
                }
            ],
        ),
        # --------------------------------------------------
        # dbReference without CURIE (prefix)
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}dbReference",
                        {"type": "CDD", "id": "cd04253"},
                    ),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_3",
            [
                {
                    "entity_id": "cdm_3",
                    "xref_type": "CDD",
                    "xref_value": "cd04253",
                    "xref": "CDD:cd04253",
                }
            ],
        ),
        # --------------------------------------------------
        # Mixed dbReferences
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}dbReference",
                        {"type": "GO", "id": "GO:0003674"},
                    ),
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}dbReference",
                        {"type": "PDB", "id": "1ABC"},
                    ),
                    entry,
                )[2]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_4",
            [
                {
                    "entity_id": "cdm_4",
                    "xref_type": "GO",
                    "xref_value": "GO:0003674",
                    "xref": "GO:0003674",
                },
                {
                    "entity_id": "cdm_4",
                    "xref_type": "PDB",
                    "xref_value": "1ABC",
                    "xref": "PDB:1ABC",
                },
            ],
        ),
        # --------------------------------------------------
        # Missing type or id
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}dbReference",
                        {"type": "GO"},  # missing id
                    ),
                    ET.SubElement(
                        entry,
                        f"{{{NS_URI}}}dbReference",
                        {"id": "123"},  # missing type
                    ),
                    entry,
                )[2]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_5",
            [],
        ),
    ],
)
def test_parse_cross_references_parametrized(build_xml, cdm_id, expected):
    entry = build_xml()

    result = parse_cross_references(entry, cdm_id)

    assert isinstance(result, list)
    assert result == expected


@pytest.mark.parametrize(
    "build_xml, cdm_id, expected",
    [
        # --------------------------------------------------
        # No accession
        # --------------------------------------------------
        (
            lambda: ET.Element(f"{{{NS_URI}}}entry"),
            "cdm_1",
            [],
        ),
        # --------------------------------------------------
        # Single accession
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(entry, f"{{{NS_URI}}}accession").__setattr__("text", "P12345"),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_2",
            [
                {
                    "entity_id": "cdm_2",
                    "identifier": "UniProt:P12345",
                    "source": "UniProt",
                    "description": "UniProt accession",
                }
            ],
        ),
        # --------------------------------------------------
        # Multiple accessions
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(entry, f"{{{NS_URI}}}accession").__setattr__("text", "Q11111"),
                    ET.SubElement(entry, f"{{{NS_URI}}}accession").__setattr__("text", "Q22222"),
                    entry,
                )[2]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_3",
            [
                {
                    "entity_id": "cdm_3",
                    "identifier": "UniProt:Q11111",
                    "source": "UniProt",
                    "description": "UniProt accession",
                },
                {
                    "entity_id": "cdm_3",
                    "identifier": "UniProt:Q22222",
                    "source": "UniProt",
                    "description": "UniProt accession",
                },
            ],
        ),
        # --------------------------------------------------
        # parse_identifiers_generic already sets source/description → setdefault
        # --------------------------------------------------
        (
            lambda: (
                lambda entry: (
                    ET.SubElement(entry, f"{{{NS_URI}}}accession").__setattr__("text", "A0A000"),
                    entry,
                )[1]
            )(ET.Element(f"{{{NS_URI}}}entry")),
            "cdm_4",
            [
                {
                    "entity_id": "cdm_4",
                    "identifier": "UniProt:A0A000",
                    "source": "UniProt",  # remains
                    "description": "UniProt accession",  # remains
                }
            ],
        ),
    ],
)
def test_parse_identifiers_parametrized(build_xml, cdm_id, expected):
    entry = build_xml()

    result = parse_identifiers(entry, cdm_id)

    assert isinstance(result, list)
    assert result == expected
