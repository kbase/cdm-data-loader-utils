import pytest
from refseq_importer.core.extractors import (
    _coalesce,
    _deep_find_str,
    _deep_collect_regex,
    extract_created_date,
    extract_assembly_name,
    extract_organism_name,
    extract_taxid,
    extract_biosample_ids,
    extract_bioproject_ids,
    extract_assembly_accessions,
    PAT_BIOSAMPLE
)


# ---------------------------------------------
# _coalesce
# ---------------------------------------------
@pytest.mark.parametrize(
    "vals, expected",
    [
        (["", "   ", "abc"], "abc"),
        (["  x ", None, ""], "x"),
        ([None, "", "   "], None),
        (["A", "B"], "A"),
    ],
)
def test_coalesce(vals, expected):
    assert _coalesce(*vals) == expected


# ---------------------------------------------
# _deep_find_str
# ---------------------------------------------
def test_deep_find_str():
    obj = {
        "level1": {
            "target": "VALUE",
            "other": [{"target": "SECOND"}],
        }
    }
    res = _deep_find_str(obj, {"target"})
    assert res == "VALUE"


# ---------------------------------------------
# _deep_collect_regex
# ---------------------------------------------
def test_deep_collect_regex():
    obj = {
        "a": "SAMN123",
        "b": ["xxx SAMN999 yyy", {"k": "SAMN555"}],
    }
    result = _deep_collect_regex(obj, PAT_BIOSAMPLE)
    assert result == ["SAMN123", "SAMN555", "SAMN999"]


# ---------------------------------------------
# extract_created_date
# ---------------------------------------------
def test_extract_created_date_refseq():
    rep = {
        "assembly_info": {
            "sourceDatabase": "SOURCE_DATABASE_REFSEQ",
            "releaseDate": "2020-01-01",
            "assemblyDate": "2019-01-01",
            "submissionDate": "2018-01-01",
        }
    }
    assert extract_created_date(rep) == "2020-01-01"


def test_extract_created_date_genbank():
    rep = {
        "assembly_info": {
            "sourceDatabase": "SOURCE_DATABASE_GENBANK",
            "submissionDate": "2018-01-01",
        }
    }
    assert extract_created_date(rep, allow_genbank_date=True) == "2018-01-01"


# ---------------------------------------------
# extract_assembly_name
# ---------------------------------------------
@pytest.mark.parametrize(
    "rep, expected",
    [
        ({"assemblyInfo": {"assemblyName": "ASM1"}}, "ASM1"),
        ({"assembly": {"assemblyName": "ASM2"}}, "ASM2"),
        ({"assembly": {"display_name": "ASM3"}}, "ASM3"),
        ({"assembly": {"displayName": "ASM4"}}, "ASM4"),
        ({"assembly": {"nested": {"display_name": "ASM5"}}}, "ASM5"),
    ],
)
def test_extract_assembly_name(rep, expected):
    assert extract_assembly_name(rep) == expected


# ---------------------------------------------
# extract_organism_name
# ---------------------------------------------
@pytest.mark.parametrize(
    "rep, expected",
    [
        ({"organism": {"scientificName": "E. coli"}}, "E. coli"),
        ({"organism": {"name": "Bacteria X"}}, "Bacteria X"),
        ({"assembly": {"organism": {"organismName": "ABC"}}}, "ABC"),
        ({"nested": {"organismName": "NNN"}}, "NNN"),
    ],
)
def test_extract_organism_name(rep, expected):
    assert extract_organism_name(rep) == expected


# ---------------------------------------------
# extract_taxid
# ---------------------------------------------
@pytest.mark.parametrize(
    "rep, expected",
    [
        ({"organism": {"taxId": 123}}, "123"),
        ({"organism": {"taxid": "456"}}, "456"),
        ({"organism": {"taxID": "789"}}, "789"),
        ({"nested": {"deep": {"taxid": 999}}}, "999"),
    ],
)
def test_extract_taxid(rep, expected):
    assert extract_taxid(rep) == expected


# ---------------------------------------------
# extract_biosample_ids
# ---------------------------------------------
@pytest.mark.parametrize(
    "rep, expected",
    [
        ({"biosample": ["BS1", "BS2"]}, ["BS1", "BS2"]),
        ({"assemblyInfo": {"biosample": {"accession": "BS3"}}}, ["BS3"]),
        ({"biosample": [{"biosampleAccession": "BS4"}]}, ["BS4"]),
        ({"text": "random SAMN111 text"}, ["SAMN111"]),
    ],
)
def test_extract_biosample_ids(rep, expected):
    assert extract_biosample_ids(rep) == expected


# ---------------------------------------------
# extract_bioproject_ids
# ---------------------------------------------
@pytest.mark.parametrize(
    "rep, expected",
    [
        ({"bioproject": ["BP1"]}, ["BP1"]),
        ({"assemblyInfo": {"bioproject": {"accession": "BP2"}}}, ["BP2"]),
        ({"bioproject": [{"bioprojectAccession": "BP3"}]}, ["BP3"]),
        ({"text": "abc PRJNA999 xyz"}, ["PRJNA999"]),
    ],
)
def test_extract_bioproject_ids(rep, expected):
    assert extract_bioproject_ids(rep) == expected


# ---------------------------------------------
# extract_assembly_accessions
# ---------------------------------------------
@pytest.mark.parametrize(
    "rep, expected_gcf, expected_gca",
    [
        (
            {"assembly": {"assembly_accession": ["GCF_0001.1"], "insdc_assembly_accession": ["GCA_0002.1"]}},
            ["GCF_0001.1"],
            ["GCA_0002.1"],
        ),
        (
            {"accession": "GCF_1111.1"},
            ["GCF_1111.1"],
            [],
        ),
        (
            {"assembly_info": {"paired_assembly": {"accession": "GCA_2222.1"}}},
            [],
            ["GCA_2222.1"],
        ),
    ],
)
def test_extract_assembly_accessions(rep, expected_gcf, expected_gca):
    gcf, gca = extract_assembly_accessions(rep)
    assert gcf == expected_gcf
    assert gca == expected_gca

