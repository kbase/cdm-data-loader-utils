import pytest
from pyspark.sql import SparkSession

from cdm_data_loader_utils.parsers.refseq_importer.core.cdm_builders import (
    build_cdm_contig_collection,
    build_cdm_entity,
    build_cdm_identifier_rows,
    build_cdm_name_rows,
    build_entity_id,
)

### pytest refseq_importer/tests/test_cdm_builders.py ###


# -------------------------------------------------------------
# Spark fixture (session scope — runs once for whole test suite)
# -------------------------------------------------------------
@pytest.fixture(scope="session")
def spark():
    spark = SparkSession.builder.master("local[1]").appName("refseq_importer_tests").getOrCreate()
    yield spark
    spark.stop()


# =============================================================
#                TEST build_entity_id
# =============================================================
@pytest.mark.parametrize("input_key", ["abc", "   hello  ", "", "123", "GCF_0001"])
def test_build_entity_id_prefix(input_key):
    eid = build_entity_id(input_key)
    assert eid.startswith("CDM:")
    assert len(eid) > 10  # UUID v5 non-empty
    assert isinstance(eid, str)


# =============================================================
#                TEST build_cdm_entity
# =============================================================
def test_build_cdm_entity_basic(spark):
    df, eid = build_cdm_entity(spark, key_for_uuid="ABC123", created_date="2020-01-01")

    row = df.collect()[0]
    assert row.entity_id == eid
    assert row.entity_type == "contig_collection"
    assert row.data_source == "RefSeq"
    assert row.created == "2020-01-01"


# =============================================================
#           TEST build_cdm_contig_collection
# =============================================================
@pytest.mark.parametrize("taxid, expected", [("1234", "NCBITaxon:1234"), (None, None), ("999", "NCBITaxon:999")])
def test_build_cdm_contig_collection_param(spark, taxid, expected):
    df = build_cdm_contig_collection(spark, entity_id="CDM:xyz", taxid=taxid)
    row = df.collect()[0]
    assert row.collection_id == "CDM:xyz"
    assert row.ncbi_taxon_id == expected


# =============================================================
#               TEST build_cdm_name_rows
# =============================================================
def test_build_cdm_name_rows(spark):
    rep = {"organism": {"name": "Escherichia coli"}, "assembly": {"display_name": "GCF_test_assembly"}}

    df = build_cdm_name_rows(spark, "CDM:abc", rep)
    rows = df.collect()

    names = {r.name for r in rows}
    assert "Escherichia coli" in names
    assert "GCF_test_assembly" in names


# =============================================================
#        TEST build_cdm_identifier_rows (parametrize!)
# =============================================================
@pytest.mark.parametrize(
    "rep, request_taxid, expected_identifiers",
    [
        # Case 1 – full fields
        (
            {"biosample": ["BS1"], "bioproject": ["BP1"], "taxid": "123"},
            "123",
            {"Biosample:BS1", "BioProject:BP1", "NCBITaxon:123"},
        ),
        # Case 2 – only taxid
        (
            {"biosample": [], "bioproject": [], "taxid": "999"},
            "999",
            {"NCBITaxon:999"},
        ),
        # Case 3 – fallback taxid used
        (
            {"biosample": ["X"], "bioproject": [], "taxid": None},
            "888",
            {"Biosample:X", "NCBITaxon:888"},
        ),
        # Case 4 – GCF/GCA accessions
        (
            {
                "biosample": ["BS"],
                "bioproject": [],
                "taxid": "555",
                "assembly": {"assembly_accession": ["GCF_0001"], "insdc_assembly_accession": ["GCA_0002"]},
            },
            "555",
            {"Biosample:BS", "NCBITaxon:555", "ncbi.assembly:GCF_0001", "insdc.gca:GCA_0002"},
        ),
    ],
)
def test_build_cdm_identifier_rows_param(rep, request_taxid, expected_identifiers):
    # Convert mock representation into what extract_assembly_accessions expects
    fake_rep = {
        "biosample": rep.get("biosample", []),
        "bioproject": rep.get("bioproject", []),
        "taxid": rep.get("taxid"),
        "assembly": rep.get("assembly", {}),
    }

    rows = build_cdm_identifier_rows("CDM:123", fake_rep, request_taxid)
    identifiers = {r["identifier"] for r in rows}

    assert identifiers == expected_identifiers
