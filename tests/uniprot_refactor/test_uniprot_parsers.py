import xml.etree.ElementTree as ET

from cdm_data_loader_utils.parsers.uniprot import parse_protein_info, parse_evidence_map

UNIPROT_NS = "https://uniprot.org/uniprot"
NSMAP = {"ns": UNIPROT_NS}


def make_entry(inner_xml: str):
    """
    Helper: wrap inner XML into a <uniprot><entry>...</entry></uniprot>
    and return the <entry> element.
    """
    xml = f"""
    <uniprot xmlns="{UNIPROT_NS}">
      <entry created="2020-01-01" modified="2020-01-02">
        {inner_xml}
      </entry>
    </uniprot>
    """
    root = ET.fromstring(xml)
    return root.find("ns:entry", NSMAP)


def test_parse_protein_info_full() -> None:
    entry = make_entry(
        """
        <protein>
          <recommendedName>
            <fullName>Example enzyme</fullName>
            <ecNumber>1.2.3.4</ecNumber>
          </recommendedName>
          <alternativeName>
            <ecNumber>5.6.7.8</ecNumber>
          </alternativeName>
        </protein>
        <proteinExistence type="evidence at protein level" />
        <sequence length="123" mass="45678" checksum="ABCDEF"
                  modified="2024-01-01" version="2">
          MPEPTIDE
        </sequence>
        """
    )

    info = parse_protein_info(entry, "CDM:123")

    # basic fields
    assert info["protein_id"] == "CDM:123"
    assert info["evidence_for_existence"] == "evidence at protein level"

    # EC numbers aggregated
    assert sorted(info["ec_numbers"]) == ["1.2.3.4", "5.6.7.8"]

    # sequence-related attributes (via get_attr / get_text)
    assert info["length"] == "123"
    assert info["mass"] == "45678"
    assert info["checksum"] == "ABCDEF"
    assert info["sequence_version"] == "2"
    assert info["sequence"].strip() == "MPEPTIDE"

    # entry modified/updated
    assert info["entry_modified"] == "2020-01-02"


def test_parse_protein_info_without_sequence() -> None:
    entry = make_entry(
        """
        <proteinExistence type="evidence at transcript level" />
        """
    )

    info = parse_protein_info(entry, "CDM:X")

    assert info["protein_id"] == "CDM:X"
    assert info["evidence_for_existence"] == "evidence at transcript level"
    # no <sequence> => no sequence-related keys
    assert "sequence" not in info
    assert "length" not in info


def test_parse_evidence_map_with_pubmed_and_other() -> None:
    entry = make_entry(
        """
        <evidence key="E1" type="ECO:0000255">
          <source>
            <dbReference type="PubMed" id="12345" />
            <dbReference type="SomeDB" id="XYZ" />
          </source>
        </evidence>
        """
    )

    ev_map = parse_evidence_map(entry)
    assert "E1" in ev_map
    ev = ev_map["E1"]

    assert ev["evidence_type"] == "ECO:0000255"
    # supporting_objects should contain non-publication cross-ref
    assert "SomeDB:XYZ" in (ev["supporting_objects"] or [])

    # we want PubMed normalized to PMID
    assert "PMID:12345" in (ev["publications"] or [])
