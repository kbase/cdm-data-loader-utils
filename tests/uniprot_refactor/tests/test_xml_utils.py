import xml.etree.ElementTree as ET

from xml_utils import (
    get_text,
    get_attr,
    parse_db_references,
    clean_dict,
)


def test_get_text_and_get_attr_basic():
    elem = ET.Element("tag", attrib={"id": "123"})
    elem.text = "  hello  "

    assert get_text(elem) == "hello"
    assert get_text(None) is None
    assert get_attr(elem, "id") == "123"
    assert get_attr(elem, "missing") is None


def test_parse_db_references_pub_and_others():
    # <source>
    #   <dbReference type="PubMed" id="12345"/>
    #   <dbReference type="DOI" id="10.1000/xyz"/>
    #   <dbReference type="PDB" id="1ABC"/>
    # </source>
    
    ns = {"ns": "dummy"}
    source = ET.Element("source")
    db1 = ET.SubElement(source, "dbReference", attrib={"type": "PubMed", "id": "12345"})
    db2 = ET.SubElement(source, "dbReference", attrib={"type": "DOI", "id": "10.1000/xyz"})
    db3 = ET.SubElement(source, "dbReference", attrib={"type": "PDB", "id": "1ABC"})

    # 手动加上 namespace 前缀，跟 parse_db_references 里的 "ns:dbReference" 对应
    db1.tag = "{dummy}dbReference"
    db2.tag = "{dummy}dbReference"
    db3.tag = "{dummy}dbReference"

    pubs, others = parse_db_references(source, ns)

    assert "PUBMED:12345" in pubs
    assert "DOI:10.1000/xyz" in pubs
    assert "PDB:1ABC" in others


def test_clean_dict_removes_nones_and_empty():
    d = {
        "a": 1,
        "b": None,
        "c": [],
        "d": {},
        "e": "ok",
    }
    cleaned = clean_dict(d)
    assert cleaned == {"a": 1, "e": "ok"}


