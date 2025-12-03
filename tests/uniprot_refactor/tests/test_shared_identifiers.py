import xml.etree.ElementTree as ET

from shared_identifiers import parse_identifiers_generic


def test_parse_identifiers_generic_basic():
    # <entry>
    #   <accession>P12345</accession>
    #   <accession>Q99999</accession>
    # </entry>
    ns = {"ns": "dummy"}
    entry = ET.Element("entry")

    a1 = ET.SubElement(entry, "accession")
    a1.text = "P12345"
    a2 = ET.SubElement(entry, "accession")
    a2.text = "Q99999"

    # Add namespace prefix to match xpath
    a1.tag = "{dummy}accession"
    a2.tag = "{dummy}accession"

    rows = parse_identifiers_generic(
        entry=entry,
        xpath="ns:accession",
        prefix="UniProt",
        ns=ns,
    )

    assert len(rows) == 2
    assert rows[0]["identifier"] == "UniProt:P12345"
    assert rows[1]["identifier"] == "UniProt:Q99999"
    assert rows[0]["source"] == "UniProt"
    assert rows[0]["description"] == "UniProt accession"

    