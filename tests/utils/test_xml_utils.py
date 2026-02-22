import xml.etree.ElementTree as ET

from cdm_data_loader_utils.utils.xml_utils import (
    clean_dict,
    get_attr,
    get_text,
    parse_db_references,
)


def test_get_text_and_get_attr_basic() -> None:
    elem = ET.Element("tag", attrib={"id": "123"})
    elem.text = "  hello  "

    assert get_text(elem) == "hello"
    assert get_text(None) is None
    assert get_attr(elem, "id") == "123"
    assert get_attr(elem, "missing") is None


def test_parse_db_references_pub_and_others() -> None:
    ns = {"ns": "dummy"}
    source = ET.Element("source")
    db1 = ET.SubElement(source, "dbReference", attrib={"type": "PubMed", "id": "12345"})
    db2 = ET.SubElement(source, "dbReference", attrib={"type": "DOI", "id": "10.1000/xyz"})
    db3 = ET.SubElement(source, "dbReference", attrib={"type": "PDB", "id": "1ABC"})

    db1.tag = "{dummy}dbReference"
    db2.tag = "{dummy}dbReference"
    db3.tag = "{dummy}dbReference"

    pubs, others = parse_db_references(source, ns)

    assert "PUBMED:12345" in pubs
    assert "DOI:10.1000/xyz" in pubs
    assert "PDB:1ABC" in others


def test_clean_dict_removes_nones_and_empty() -> None:
    """Test that clean_dict removes None and empty values."""
    d = {
        "a": 1,
        "b": None,
        "c": [],
        "d": {},
        "e": "ok",
    }
    cleaned = clean_dict(d)
    assert cleaned == {"a": 1, "e": "ok"}


# # TODO: add gzipped file
# @pytest.mark.parametrize("file", ["example.xml", "uniref_100.xml"])
# def test_parse_head_matter(test_data_dir: Path, file: str) -> None:
#     """Test the extraction of metadata from the head of an xml file."""
#     example_uniref_file = test_data_dir / "uniprot" / "uniref" / file
#     data_src = parse_head_matter(example_uniref_file)
#     assert data_src == {
#         "version": "2025_04",
#         "date_published": "2025-10-08",
#         "name": "UniRef50",
#     }


class FakeSparkDF:
    """A fake DataFrame returned by spark.read.format().load().select()."""

    def __init__(self, rows):
        self._rows = rows

    def collect(self):
        return self._rows


class FakeSparkReader:
    """Mock spark.read.format('delta').load().select() chain."""

    def __init__(self, rows=None, fail=False):
        self._rows = rows
        self._fail = fail

    def format(self, fmt):
        assert fmt == "delta"
        return self

    def load(self, path):
        if self._fail:
            raise Exception("Table does not exist")
        return self

    def select(self, *cols):
        return FakeSparkDF(self._rows)
