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
    PYTHONPATH=. pytest test_uniprot.py

"""

import json
import uuid
from pathlib import Path
import pytest
from unittest.mock import Mock
import requests

from cdm_data_loader_utils.parsers.uniprot import (
    download_file,
    save_datasource_record,
    parse_identifiers,
    prepare_local_xml,
    stable_cdm_id_from_accession,
    CDM_UUID_NAMESPACE,
    load_existing_identifiers,
    parse_names,
    parse_evidence_map,
    parse_associations,
    parse_publications,
)


@pytest.mark.parametrize(
    "xml_url",
    [
        "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz",
        "http://example.org/uniprot_test.xml.gz",
    ],
)
def test_save_datasource_record(tmp_path: Path, xml_url: str):
    """
    save_datasource_record should:
    - create datasource.json
    - return the same content as written to disk
    """

    output_dir = tmp_path / "output"
    result = save_datasource_record(xml_url, str(output_dir))

    # return value
    assert isinstance(result, dict)
    assert result["url"] == xml_url
    assert result["source"] == "UniProt"

    # file existence
    output_file = output_dir / "datasource.json"
    assert output_file.exists()

    # file content correctness
    with open(output_file, encoding="utf-8") as f:
        on_disk = json.load(f)

    assert on_disk == result


# ---------------------------------------------------------
# Helpers
# ---------------------------------------------------------
class FakeResponse:
    def __init__(self, content=b"test-data", status_code=200):
        self.content = content
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size=8192):
        yield self.content

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        pass


# ---------------------------------------------------------
# Tests
# ---------------------------------------------------------
def test_download_file_skip_if_exists(tmp_path, monkeypatch):
    output = tmp_path / "file.bin"
    output.write_bytes(b"existing")

    mock_get = Mock()
    monkeypatch.setattr(requests, "get", mock_get)

    download_file(
        url="http://example.com/file",
        output_path=str(output),
        overwrite=False,
    )

    # requests.get should not be called
    mock_get.assert_not_called()
    assert output.read_bytes() == b"existing"


def test_download_file_success(tmp_path, monkeypatch):
    """Successful download writes file"""
    output = tmp_path / "file.bin"

    monkeypatch.setattr(
        requests,
        "get",
        lambda *args, **kwargs: FakeResponse(b"hello"),
    )

    download_file(
        url="http://example.com/file",
        output_path=str(output),
        overwrite=True,
    )

    assert output.exists()
    assert output.read_bytes() == b"hello"


@pytest.mark.parametrize("status_code", [404, 500])
def test_download_file_http_error(tmp_path, monkeypatch, status_code):
    """HTTP error => exception + file removed"""
    output = tmp_path / "file.bin"

    monkeypatch.setattr(
        requests,
        "get",
        lambda *args, **kwargs: FakeResponse(
            content=b"bad",
            status_code=status_code,
        ),
    )

    with pytest.raises(requests.HTTPError):
        download_file(
            url="http://example.com/file",
            output_path=str(output),
            overwrite=True,
        )

    assert not output.exists()


def test_download_file_partial_cleanup(tmp_path, monkeypatch):
    """Exception during iter_content => partial file cleaned"""

    class BrokenResponse(FakeResponse):
        def iter_content(self, chunk_size=8192):
            yield b"partial"
            raise RuntimeError("connection lost")

    output = tmp_path / "file.bin"

    monkeypatch.setattr(
        requests,
        "get",
        lambda *args, **kwargs: BrokenResponse(),
    )

    with pytest.raises(RuntimeError):
        download_file(
            url="http://example.com/file",
            output_path=str(output),
            overwrite=True,
        )

    assert not output.exists()


def test_prepare_local_xml_download_called(monkeypatch, tmp_path):
    """
    Verify:
    1. output directory is created
    2. download_file is called with correct arguments
    3. returned path is correct
    """

    xml_url = "https://example.org/uniprot_sprot.xml.gz"
    output_dir = tmp_path / "output"

    called = {}

    def fake_download(url, path, *args, **kwargs):
        # record calls instead of real download
        called["url"] = url
        called["path"] = path

        # simulate downloaded file
        Path(path).write_bytes(b"fake gzip content")

    # monkeypatch download_file inside uniprot module
    monkeypatch.setattr("uniprot.download_file", fake_download)

    local_path = prepare_local_xml(xml_url, str(output_dir))
    expected_path = output_dir / "uniprot_sprot.xml.gz"

    assert output_dir.exists()
    assert Path(local_path) == expected_path
    assert called["url"] == xml_url
    assert Path(called["path"]) == expected_path
    assert expected_path.exists()


def test_stable_cdm_id_is_deterministic():
    accession = "P12345"

    id1 = stable_cdm_id_from_accession(accession)
    id2 = stable_cdm_id_from_accession(accession)

    assert id1 == id2


def test_stable_cdm_id_diff_accessions():
    id1 = stable_cdm_id_from_accession("P12345")
    id2 = stable_cdm_id_from_accession("Q99999")

    assert id1 != id2


def test_stable_cdm_id_prefix():
    accession = "P12345"
    prefix = "cdm_test_"

    result = stable_cdm_id_from_accession(accession, prefix=prefix)

    assert result.startswith(prefix)


def test_stable_cdm_id_uuid_v5_correctness():
    accession = "P12345"

    expected_uuid = uuid.uuid5(CDM_UUID_NAMESPACE, accession)
    result = stable_cdm_id_from_accession(accession)

    assert result == f"cdm_prot_{expected_uuid}"


@pytest.mark.parametrize(
    "accession",
    [
        "P12345",
        "Q9XYZ1",
        "A0A1234567",
    ],
)
def test_stable_cdm_id_parametrized(accession):
    id1 = stable_cdm_id_from_accession(accession)
    id2 = stable_cdm_id_from_accession(accession)

    assert id1 == id2


class FakeRow(dict):
    """Simple stand-in for pyspark.sql.Row"""

    def __getattr__(self, item):
        return self[item]


class FakeDataFrame:
    def __init__(self, rows):
        self._rows = rows

    def select(self, *cols):
        return self

    def collect(self):
        return self._rows


class FakeReader:
    def __init__(self, rows):
        self.rows = rows

    def format(self, _):
        return self

    def load(self, _):
        return FakeDataFrame(self.rows)


class FakeSpark:
    def __init__(self, rows):
        self.read = FakeReader(rows)


@pytest.mark.parametrize(
    "rows, expected",
    [
        # ---------------- valid UniProt ----------------
        (
            [
                FakeRow(identifier="UniProt:P12345", entity_id="cdm1"),
                FakeRow(identifier="UniProt:Q99999", entity_id="cdm2"),
            ],
            {
                "P12345": "cdm1",
                "Q99999": "cdm2",
            },
        ),
        # ---------------- mixed prefixes ----------------
        (
            [
                FakeRow(identifier="UniProt:P12345", entity_id="cdm1"),
                FakeRow(identifier="RefSeq:WP_000001", entity_id="cdmX"),
            ],
            {
                "P12345": "cdm1",
            },
        ),
        # ---------------- invalid formats ----------------
        (
            [
                FakeRow(identifier="UniProtP12345", entity_id="cdm1"),  # no colon
                FakeRow(identifier=None, entity_id="cdm2"),
            ],
            {},
        ),
        # ---------------- case-insensitive prefix ----------------
        (
            [
                FakeRow(identifier="uniprot:A0A123", entity_id="cdmA"),
            ],
            {
                "A0A123": "cdmA",
            },
        ),
    ],
)
def test_load_existing_identifiers_parametrized(tmp_path, rows, expected):
    output_dir = tmp_path
    namespace = "uniprot_db"

    # create fake identifiers dir
    id_dir = output_dir / namespace / "identifiers"
    id_dir.mkdir(parents=True)

    fake_spark = FakeSpark(rows)

    result = load_existing_identifiers(
        spark=fake_spark,
        output_dir=str(output_dir),
        namespace=namespace,
    )

    assert result == expected


@pytest.mark.parametrize(
    "generic_output, cdm_id, expected",
    [
        # case 1: single identifier
        (
            [
                {
                    "identifier": "P12345",
                    "identifier_type": "UniProt",
                }
            ],
            "cdm_prot_test_1",
            [
                {
                    "identifier": "P12345",
                    "identifier_type": "UniProt",
                    "entity_id": "cdm_prot_test_1",
                }
            ],
        ),
        # case 2: multiple identifiers
        (
            [
                {"identifier": "P1", "identifier_type": "UniProt"},
                {"identifier": "P2", "identifier_type": "UniProt"},
            ],
            "cdm_prot_test_2",
            [
                {
                    "identifier": "P1",
                    "identifier_type": "UniProt",
                    "entity_id": "cdm_prot_test_2",
                },
                {
                    "identifier": "P2",
                    "identifier_type": "UniProt",
                    "entity_id": "cdm_prot_test_2",
                },
            ],
        ),
        # case 3: empty list
        (
            [],
            "cdm_prot_test_3",
            [],
        ),
    ],
)
def test_parse_identifiers_adds_entity_id(monkeypatch, generic_output, cdm_id, expected):
    """
    parse_identifiers should:
    - call parse_identifiers_generic
    - add entity_id to each returned row
    """

    def mock_parse_identifiers_generic(*args, **kwargs):
        return generic_output

    monkeypatch.setattr(
        "uniprot.parse_identifiers_generic",
        mock_parse_identifiers_generic,
    )

    dummy_entry = object()
    result = parse_identifiers(dummy_entry, cdm_id)

    assert result == expected


@pytest.mark.parametrize(
    "top_level_names, protein_blocks, expected",
    [
        # -------------------------
        # case 1: only top-level
        # -------------------------
        (
            ["Protein A"],
            None,
            [
                {
                    "entity_id": "cdm_prot_1",
                    "name": "Protein A",
                    "description": "UniProt protein name",
                }
            ],
        ),
        # -------------------------
        # case 2: only protein names
        # -------------------------
        (
            [],
            [
                ("recommended", "full", "Protein Rec Full"),
                ("recommended", "short", "Protein Rec Short"),
                ("alternative", "full", "Protein Alt Full"),
            ],
            [
                {
                    "entity_id": "cdm_prot_1",
                    "name": "Protein Rec Full",
                    "description": "UniProt recommended full name",
                },
                {
                    "entity_id": "cdm_prot_1",
                    "name": "Protein Rec Short",
                    "description": "UniProt recommended short name",
                },
                {
                    "entity_id": "cdm_prot_1",
                    "name": "Protein Alt Full",
                    "description": "UniProt alternative full name",
                },
            ],
        ),
        # -------------------------
        # case 3: top-level + protein
        # -------------------------
        (
            ["Protein A"],
            [
                ("recommended", "full", "Protein Rec Full"),
            ],
            [
                {
                    "entity_id": "cdm_prot_1",
                    "name": "Protein A",
                    "description": "UniProt protein name",
                },
                {
                    "entity_id": "cdm_prot_1",
                    "name": "Protein Rec Full",
                    "description": "UniProt recommended full name",
                },
            ],
        ),
        # -------------------------
        # case 4: nothing
        # -------------------------
        (
            [],
            None,
            [],
        ),
    ],
)
def test_parse_names(
    monkeypatch,
    top_level_names,
    protein_blocks,
    expected,
):
    cdm_id = "cdm_prot_1"

    # -------------------------------------------------
    # mock: find_all_text
    # -------------------------------------------------
    monkeypatch.setattr(
        "uniprot.find_all_text",
        lambda entry, xpath, ns: top_level_names,
    )

    # -------------------------------------------------
    # Dummy element returned by get_text
    # -------------------------------------------------
    class DummyTextElem:
        def __init__(self, text):
            self.text = text

    # -------------------------------------------------
    # Name block: handles fullName / shortName
    # -------------------------------------------------
    class DummyNameBlock:
        def __init__(self, logical, values):
            self.logical = logical
            self.values = values

        def find(self, xpath, ns):
            if "fullName" in xpath:
                return DummyTextElem(self.values.get("full"))
            if "shortName" in xpath:
                return DummyTextElem(self.values.get("short"))
            return None

    # -------------------------------------------------
    # Protein element
    # -------------------------------------------------
    class DummyProtein:
        def findall(self, xpath, ns):
            if protein_blocks is None:
                return []

            blocks = {}
            for logical, length, text in protein_blocks:
                blocks.setdefault(logical, {})[length] = text

            if "recommendedName" in xpath and "recommended" in blocks:
                return [DummyNameBlock("recommended", blocks["recommended"])]
            if "alternativeName" in xpath and "alternative" in blocks:
                return [DummyNameBlock("alternative", blocks["alternative"])]

            return []

    # -------------------------------------------------
    # Entry
    # -------------------------------------------------
    class DummyEntry:
        def find(self, xpath, ns):
            if "protein" in xpath and protein_blocks is not None:
                return DummyProtein()
            return None

    entry = DummyEntry()

    # -------------------------------------------------
    # mock: get_text
    # -------------------------------------------------
    monkeypatch.setattr(
        "uniprot.get_text",
        lambda elem: elem.text if elem else None,
    )

    # -------------------------------------------------
    # mock: _make_name_record
    # -------------------------------------------------
    monkeypatch.setattr(
        "uniprot._make_name_record",
        lambda entity_id, name, desc: {
            "entity_id": entity_id,
            "name": name,
            "description": desc,
        },
    )

    result = parse_names(entry, cdm_id)

    assert result == expected


@pytest.mark.parametrize(
    "evidence_blocks, expected",
    [
        # -------------------------------------------------
        # case 1: no evidence
        # -------------------------------------------------
        (
            [],
            {},
        ),
        # -------------------------------------------------
        # case 2: evidence without key
        # -------------------------------------------------
        (
            [
                {
                    "key": None,
                    "type": "experimental",
                    "pubs": ["PUBMED:123"],
                    "others": [],
                }
            ],
            {},
        ),
        # -------------------------------------------------
        # case 3: evidence with key, no source
        # -------------------------------------------------
        (
            [
                {
                    "key": "E1",
                    "type": "experimental",
                    "pubs": None,
                    "others": None,
                }
            ],
            {
                "E1": {
                    "evidence_type": "experimental",
                }
            },
        ),
        # -------------------------------------------------
        # case 4: evidence with source + PUBMED normalize
        # -------------------------------------------------
        (
            [
                {
                    "key": "E2",
                    "type": "computational",
                    "pubs": ["PUBMED:12345", "DOI:10.1000/xyz"],
                    "others": ["CHEBI:1234"],
                }
            ],
            {
                "E2": {
                    "evidence_type": "computational",
                    "publications": ["PMID:12345", "DOI:10.1000/xyz"],
                    "supporting_objects": ["CHEBI:1234"],
                }
            },
        ),
    ],
)
def test_parse_evidence_map(monkeypatch, evidence_blocks, expected):
    """
    parse_evidence_map should:
    - skip evidence without key
    - normalize PUBMED -> PMID
    - preserve supporting objects
    - clean empty fields
    """

    # =====================================================
    # Dummy XML elements
    # =====================================================
    class DummyEvidence:
        def __init__(self, block):
            self.block = block

        def find(self, xpath, ns):
            if "source" in xpath and self.block.get("pubs") is not None:
                return object()
            return None

    class DummyEntry:
        def findall(self, xpath, ns):
            if "evidence" in xpath:
                return [DummyEvidence(b) for b in evidence_blocks]
            return []

    entry = DummyEntry()

    # =====================================================
    # mock: get_attr
    # =====================================================
    def mock_get_attr(elem, key):
        return elem.block.get(key)

    monkeypatch.setattr(
        "uniprot.get_attr",
        mock_get_attr,
    )

    # =====================================================
    # mock: parse_db_references
    # =====================================================
    def mock_parse_db_references(source, ns):
        # find which evidence this source belongs to
        for b in evidence_blocks:
            if b.get("pubs") is not None:
                return b.get("pubs", []), b.get("others", [])
        return [], []

    monkeypatch.setattr(
        "uniprot.parse_db_references",
        mock_parse_db_references,
    )

    # =====================================================
    # mock: clean_dict
    # =====================================================
    monkeypatch.setattr(
        "uniprot.clean_dict",
        lambda d: {k: v for k, v in d.items() if v is not None},
    )

    result = parse_evidence_map(entry)

    assert result == expected


@pytest.mark.parametrize(
    "organism_taxid, dbrefs, comments, expected",
    [
        # -------------------------------------------------
        # case 1: only taxonomy
        # -------------------------------------------------
        (
            "9606",
            [],
            [],
            [
                {
                    "entity_id": "cdm_prot_1",
                    "target": "NCBITaxon:9606",
                }
            ],
        ),
        # -------------------------------------------------
        # case 2: generic dbReferences
        # -------------------------------------------------
        (
            None,
            [
                {"type": "UniRef100", "id": "P12345", "evidence": "E1"},
                {"type": "PDB", "id": "1ABC", "evidence": None},
            ],
            [],
            [
                {
                    "entity_id": "cdm_prot_1",
                    "target": "UniRef100:P12345",
                    "evidence_key": "E1",
                },
                {
                    "entity_id": "cdm_prot_1",
                    "target": "PDB:1ABC",
                },
            ],
        ),
        # -------------------------------------------------
        # case 3: catalytic activity + cofactor
        # -------------------------------------------------
        (
            None,
            [],
            [
                {"type": "catalytic activity", "items": ["R1", "R2"]},
                {"type": "cofactor", "items": ["C1"]},
            ],
            [
                {"assoc": "reaction:R1"},
                {"assoc": "reaction:R2"},
                {"assoc": "cofactor:C1"},
            ],
        ),
        # -------------------------------------------------
        # case 4: everything together
        # -------------------------------------------------
        (
            "562",
            [
                {"type": "UniRef90", "id": "Q9XYZ", "evidence": "E2"},
            ],
            [
                {"type": "catalytic activity", "items": ["R1"]},
            ],
            [
                {
                    "entity_id": "cdm_prot_1",
                    "target": "NCBITaxon:562",
                },
                {
                    "entity_id": "cdm_prot_1",
                    "target": "UniRef90:Q9XYZ",
                    "evidence_key": "E2",
                },
                {"assoc": "reaction:R1"},
            ],
        ),
    ],
)
def test_parse_associations(
    monkeypatch,
    organism_taxid,
    dbrefs,
    comments,
    expected,
):
    cdm_id = "cdm_prot_1"
    evidence_map = {"E1": {}, "E2": {}}

    class DummyElem:
        def __init__(self, attrs=None):
            self.attrs = attrs or {}

        def get(self, key):
            return self.attrs.get(key)

    class DummyOrganism:
        def find(self, xpath, ns):
            if "NCBI Taxonomy" in xpath and organism_taxid:
                return DummyElem({"id": organism_taxid})
            return None

    class DummyDBRef(DummyElem):
        pass

    class DummyComment(DummyElem):
        def __init__(self, comment_type, items):
            super().__init__({"type": comment_type})
            self.items = items

        def findall(self, xpath, ns):
            return self.items

    class DummyEntry:
        def find(self, xpath, ns):
            if "organism" in xpath and organism_taxid:
                return DummyOrganism()
            return None

        def findall(self, xpath, ns):
            if xpath == "ns:dbReference":
                return [
                    DummyDBRef({
                        "type": d["type"],
                        "id": d["id"],
                        "evidence": d.get("evidence"),
                    })
                    for d in dbrefs
                ]
            if xpath == "ns:comment":
                return [DummyComment(c["type"], c["items"]) for c in comments]
            return []

    entry = DummyEntry()

    def mock_make_association(entity_id, target, evidence_key=None, evidence_map=None):
        out = {
            "entity_id": entity_id,
            "target": target,
        }
        if evidence_key is not None:
            out["evidence_key"] = evidence_key
        return out

    monkeypatch.setattr("uniprot._make_association", mock_make_association)
    monkeypatch.setattr(
        "uniprot.parse_reaction_association",
        lambda reaction, cdm_id, evidence_map: [{"assoc": f"reaction:{reaction}"}],
    )
    monkeypatch.setattr(
        "uniprot.parse_cofactor_association",
        lambda cofactor, cdm_id: [{"assoc": f"cofactor:{cofactor}"}],
    )

    result = parse_associations(entry, cdm_id, evidence_map)

    assert result == expected


@pytest.mark.parametrize(
    "references, expected",
    [
        # -------------------------------------------------
        # case 1: no reference
        # -------------------------------------------------
        (
            [],
            [],
        ),
        # -------------------------------------------------
        # case 2: reference without citation (ignored)
        # -------------------------------------------------
        (
            [
                {"pubs": ["PUBMED:12345"]},
            ],
            [],
        ),
        # -------------------------------------------------
        # case 3: single PUBMED
        # -------------------------------------------------
        (
            [
                {"citation": True, "pubs": ["PUBMED:12345"]},
            ],
            ["PMID:12345"],
        ),
        # -------------------------------------------------
        # case 4: PUBMED + DOI + duplicates
        # -------------------------------------------------
        (
            [
                {
                    "citation": True,
                    "pubs": ["PUBMED:12345", "DOI:10.1000/xyz"],
                },
                {
                    "citation": True,
                    "pubs": ["PUBMED:12345"],
                },
            ],
            ["PMID:12345", "DOI:10.1000/xyz"],
        ),
        # -------------------------------------------------
        # case 5: unsupported prefixes ignored
        # -------------------------------------------------
        (
            [
                {
                    "citation": True,
                    "pubs": ["PMC:999", "ISBN:123", "DOI:10.1/abc"],
                },
            ],
            ["DOI:10.1/abc"],
        ),
    ],
)
def test_parse_publications(monkeypatch, references, expected):
    """
    parse_publications should:
    - normalize PUBMED -> PMID
    - normalize DOI
    - ignore unsupported prefixes
    - deduplicate while preserving order
    """

    # =====================================================
    # Dummy elements
    # =====================================================
    class DummyReference:
        def __init__(self, block):
            self.block = block

        def find(self, xpath, ns):
            if "citation" in xpath and self.block.get("citation"):
                return object()
            return None

    class DummyEntry:
        def findall(self, xpath, ns):
            if "reference" in xpath:
                return [DummyReference(b) for b in references]
            return []

    entry = DummyEntry()

    # =====================================================
    # mock: parse_db_references
    # =====================================================
    def mock_parse_db_references(citation, ns):
        for b in references:
            if b.get("citation"):
                return b.get("pubs", []), []
        return [], []

    monkeypatch.setattr(
        "uniprot.parse_db_references",
        mock_parse_db_references,
    )

    result = parse_publications(entry)

    assert result == expected
