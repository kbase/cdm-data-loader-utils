import uuid
import click
import pytest
import requests
import pandas as pd
from datetime import date, datetime
from unittest.mock import patch, MagicMock
from refseq_api import fetch_reports_by_taxon  
from refseq_api import _coalesce, _deep_find_str, _deep_collect_regex, PAT_BIOSAMPLE, PAT_BIOPROJECT, PAT_GCF, PAT_GCA
from refseq_api import extract_created_date
from refseq_api import extract_assembly_name
from refseq_api import extract_organism_name
from refseq_api import extract_taxid
from refseq_api import extract_biosample_ids
from refseq_api import extract_bioproject_ids
from refseq_api import extract_assembly_accessions
from refseq_api import build_cdm_datasource
from refseq_api import build_entity_id, CDM_NAMESPACE
from refseq_api import build_cdm_entity
from refseq_api import build_cdm_contig_collection
from refseq_api import build_cdm_name_rows
from refseq_api import build_cdm_identifier_rows
from refseq_api import process_report, process_taxon, finalize_tables, write_and_preview, main
from refseq_api import parse_taxid_args


@pytest.mark.parametrize(
    "taxid_arg, taxid_file_content, expected",
    [
        ("224325", None, ["224325"]),
        ("224325,2741724", None, ["224325", "2741724"]),
        ("TaxID:224325, abc2741724", None, ["224325", "2741724"]),
        ("224325,224325,2741724", None, ["224325", "2741724"]),
        ("", None, [])]
)


def test_parse_taxid_args_inline(taxid_arg, taxid_file_content, expected, tmp_path):
    taxid_file = None
    if taxid_file_content:
        taxid_file = tmp_path / "taxids.txt"
        taxid_file.write_text(taxid_file_content)
    assert parse_taxid_args(taxid_arg, taxid_file) == expected


def test_parse_taxid_file_not_found():
    with pytest.raises(click.BadParameter):
        parse_taxid_args(None, "nonexistent.txt")


def make_response(reports, next_token=None):
    payload = {"reports": reports}
    if next_token:
        payload["next_page_token"] = next_token
    mock_resp = MagicMock()
    mock_resp.json.return_value = payload
    mock_resp.raise_for_status.return_value = None
    return mock_resp


@pytest.mark.parametrize(
    "mock_reports, side_effect, expected_accessions",
    [
        # RefSeq only: should keep only RefSeq record
        (
            [
                {"accession": "GCA_123", "assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_GENBANK"}},
                {"accession": "GCF_456", "assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_REFSEQ"}},
            ],
            None,
            ["GCF_456"]
        ),
        # Empty reports
        ([], None, []),
        # Pagination
        (
            [
                {"accession": "GCF_1", "assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_REFSEQ"}}
            ],
            None,  # side_effect controlled separately
            ["GCF_1", "GCF_2"]
        ),
    ]
)

@patch("refseq_api.requests.get")

def test_fetch_reports(mock_get, mock_reports, side_effect, expected_accessions):
    if expected_accessions == ["GCF_1", "GCF_2"]:
        first_page = make_response(mock_reports, next_token="token123")
        second_page = make_response(
            [{"accession": "GCF_2", "assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_REFSEQ"}}]
        )
        mock_get.side_effect = [first_page, second_page]
    elif side_effect:
        mock_get.side_effect = side_effect
    else:
        mock_get.return_value = make_response(mock_reports)

    results = list(fetch_reports_by_taxon("1234"))
    accessions = [r["accession"] for r in results]
    assert accessions == expected_accessions


# Network error 
@patch("refseq_api.requests.get")

def test_network_error(mock_get):
    mock_get.side_effect = requests.RequestException("Network down")
    results = list(fetch_reports_by_taxon("1234"))
    assert results == []


# -------------------- _coalesce --------------------
@pytest.mark.parametrize("inputs, expected", [
    (["", " ", None, "abc"], "abc"),   # should pick first non-empty
    ([None, "   ", "xyz", "zzz"], "xyz"),
    ([None, "", "   "], None),         # all empty
])

def test_coalesce(inputs, expected):
    assert _coalesce(*inputs) == expected


# -------------------- _deep_find_str --------------------
@pytest.mark.parametrize("obj, keys, expected", [
    ({"assemblyDate": "2000-12-01"}, {"assemblyDate"}, "2000-12-01"),
    ({"nested": {"releaseDate": "2010-01-01"}}, {"releaseDate"}, "2010-01-01"),
    ({"list": [{"submissionDate": "2020-05-05"}]}, {"submissionDate"}, "2020-05-05"),
    ({"noDate": "xxx"}, {"releaseDate"}, None),
])

def test_deep_find_str(obj, keys, expected):
    assert _deep_find_str(obj, keys) == expected


# -------------------- _deep_collect_regex --------------------
@pytest.mark.parametrize("obj, pattern, expected", [
    ("Biosample SAMN12345 here", PAT_BIOSAMPLE, ["SAMN12345"]),
    ({"a": "Project PRJNA99999"}, PAT_BIOPROJECT, ["PRJNA99999"]),
    (["Genome GCF_000123456.1"], PAT_GCF, ["GCF_000123456.1"]),
    ({"list": ["Some GCA_000987654.2", "Other GCA_000987654.2"]}, PAT_GCA, ["GCA_000987654.2"]),  # dedup
    ("No matches here", PAT_BIOSAMPLE, []),
])

def test_deep_collect_regex(obj, pattern, expected):
    assert _deep_collect_regex(obj, pattern) == expected


@pytest.mark.parametrize(
    "rep, allow_genbank, expected",
    [
        # --- RefSeq record: release_date present ---
        (
            {"assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_REFSEQ", "releaseDate": "2020-01-01"}},
            False,
            "2020-01-01",
        ),
        # --- RefSeq record: no release_date, use assembly_date ---
        (
            {"assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_REFSEQ", "assemblyDate": "2021-05-05"}},
            False,
            "2021-05-05",
        ),
        # --- RefSeq record: no release/assembly, use submission_date ---
        (
            {"assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_REFSEQ", "submissionDate": "2022-03-03"}},
            False,
            "2022-03-03",
        ),
        # --- GenBank record: allow_genbank_date = False, expect None ---
        (
            {"assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_GENBANK", "submissionDate": "2019-09-09"}},
            False,
            None,
        ),
        # --- GenBank record: allow_genbank_date = True, expect submission_date ---
        (
            {"assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_GENBANK", "submissionDate": "2019-09-09"}},
            True,
            "2019-09-09",
        ),
        # --- No dates at all ---
        (
            {"assemblyInfo": {"sourceDatabase": "SOURCE_DATABASE_REFSEQ"}},
            False,
            None,
        ),
    ]
)

def test_extract_created_date(rep, allow_genbank, expected):
    result = extract_created_date(rep, allow_genbank_date=allow_genbank, debug=True)
    assert result == expected


@pytest.mark.parametrize(
    "rep, expected",
    [
        # Top-level assemblyName
        ({"assemblyName": "ASM1234v1"}, "ASM1234v1"),
        # Inside assemblyInfo
        ({"assemblyInfo": {"assemblyName": "ASM5678v1"}}, "ASM5678v1"),
        # Inside assembly
        ({"assembly": {"assemblyName": "ASM9999v1"}}, "ASM9999v1"),
        # Deep nested (fallback to _deep_find_str)
        ({"meta": {"nested": {"assembly_name": "ASM_DEEPv1"}}}, "ASM_DEEPv1"),
        # No name available
        ({"assemblyInfo": {}, "assembly": {}}, None),
    ]
)

def test_extract_assembly_name(rep, expected):
    assert extract_assembly_name(rep) == expected


@pytest.mark.parametrize(
    "rep, expected",
    [
        # Top-level organismName
        ({"organism": {"organismName": "Haloferax volcanii"}}, "Haloferax volcanii"),
        # Top-level scientificName
        ({"organism": {"scientificName": "Methanococcus maripaludis"}}, "Methanococcus maripaludis"),
        # Top-level taxName
        ({"organism": {"taxName": "Archaeoglobus fulgidus"}}, "Archaeoglobus fulgidus"),
        # assemblyInfo.organism.organismName
        ({"assemblyInfo": {"organism": {"organismName": "Sulfolobus islandicus"}}}, "Sulfolobus islandicus"),
        # assembly.organism.organismName
        ({"assembly": {"organism": {"organismName": "Thermococcus kodakarensis"}}}, "Thermococcus kodakarensis"),
        # Deep nested key (should trigger _deep_find_str)
        ({"meta": {"data": {"deep": {"scientificName": "Pyrococcus furiosus"}}}}, "Pyrococcus furiosus"),
        # No name available
        ({"assemblyInfo": {}, "assembly": {}, "organism": {}}, None)]
)

def test_extract_organism_name(rep, expected):
    assert extract_organism_name(rep) == expected


@pytest.mark.parametrize(
    "rep, expected",
    [
        # Top-level taxId as int
        ({"organism": {"taxId": 12345}}, "12345"),

        # Top-level taxid as string
        ({"organism": {"taxid": "67890"}}, "67890"),

        # Top-level taxID as float
        ({"organism": {"taxID": 11111.0}}, "11111"),

        # Nested dict with taxid
        ({"meta": {"organism": {"taxid": "22222"}}}, "22222"),

        # Nested deeply inside list
        ({"organisms": [{"nested": {"tax_id": "33333"}}]}, "33333"),

        # Invalid taxid (non-digit string)
        ({"organism": {"taxId": "abc"}}, None),

        # Missing entirely
        ({"organism": {}}, None),
    ]
)

def test_extract_taxid(rep, expected):
    assert extract_taxid(rep) == expected


@pytest.mark.parametrize(
    "rep, expected",
    [
        # biosample in assemblyInfo dict
        (
            {"assemblyInfo": {"biosample": {"accession": "SAMN12345"}}},
            ["SAMN12345"]
        ),

        # biosample in assembly dict
        (
            {"assembly": {"biosample": {"biosampleAccession": "SAMN67890"}}},
            ["SAMN67890"]
        ),

        # biosample at top level
        (
            {"biosample": {"accession": "SAMN11111"}},
            ["SAMN11111"]
        ),

        # biosample as list of dicts
        (
            {"biosample": [{"accession": "SAMN22222"}, {"biosampleAccession": "SAMN33333"}]},
            ["SAMN22222", "SAMN33333"]
        ),

        # fallback to regex
        (
            {"note": "Sample accession SAMN44444 found in text"},
            ["SAMN44444"]
        ),

        # no biosample at all
        (
            {"assemblyInfo": {}, "assembly": {}},
            []
        ),
    ]
)

def test_extract_biosample_ids(rep, expected):
    result = extract_biosample_ids(rep)
    assert result == expected


@pytest.mark.parametrize(
    "rep, expected",
    [
        # bioproject in assemblyInfo dict
        (
            {"assemblyInfo": {"bioproject": {"accession": "PRJNA12345"}}},
            ["PRJNA12345"]
        ),

        # bioproject in assembly dict
        (
            {"assembly": {"bioproject": {"bioprojectAccession": "PRJNA67890"}}},
            ["PRJNA67890"]
        ),

        # bioproject at top level
        (
            {"bioproject": {"accession": "PRJNA11111"}},
            ["PRJNA11111"]
        ),

        # bioproject as list of dicts
        (
            {"bioproject": [{"accession": "PRJNA22222"}, {"bioprojectAccession": "PRJNA33333"}]},
            ["PRJNA22222", "PRJNA33333"]
        ),

        # fallback to regex
        (
            {"note": "Bioproject accession PRJNA44444 found in free text"},
            ["PRJNA44444"]
        ),

        # no bioproject at all
        (
            {"assemblyInfo": {}, "assembly": {}},
            []
        ),
    ]
)

def test_extract_bioproject_ids(rep, expected):
    result = extract_bioproject_ids(rep)
    assert result == expected


@pytest.mark.parametrize(
    "rep, expected_gcf, expected_gca",
    [
        # Top-level accession is GCF
        ({"accession": "GCF_000123456.1"}, ["GCF_000123456.1"], []),

        # Top-level accession is GCA
        ({"accession": "GCA_000654321.1"}, [], ["GCA_000654321.1"]),

        # Paired assembly is GCF
        (
            {"assemblyInfo": {"paired_assembly": {"accession": "GCF_999999999.1"}}},
            ["GCF_999999999.1"],
            []
        ),

        # Paired assembly is GCA
        (
            {"assembly_info": {"paired_assembly": {"accession": "GCA_888888888.1"}}},
            [],
            ["GCA_888888888.1"]
        ),

        # Both top-level and paired (mixed GCF + GCA)
        (
            {
                "accession": "GCF_000123456.1",
                "assemblyInfo": {"paired_assembly": {"accession": "GCA_000654321.1"}},
            },
            ["GCF_000123456.1"],
            ["GCA_000654321.1"],
        ),

        # Invalid accession (not GCF/GCA) → ignore
        ({"accession": "XYZ_123"}, [], []),

        # Empty dict
        ({}, [], []),
    ]
)

def test_extract_assembly_accessions(rep, expected_gcf, expected_gca):
    gcf, gca = extract_assembly_accessions(rep)
    assert gcf == expected_gcf
    assert gca == expected_gca


@pytest.mark.parametrize(
    "expected_name, expected_source, expected_url, expected_version",
    [
        (
            "RefSeq",
            "NCBI RefSeq",
            "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/",
            "231",
        )
    ]
)

def test_build_cdm_datasource(expected_name, expected_source, expected_url, expected_version):
    df = build_cdm_datasource()

    # Should return exactly one row
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1

    record = df.iloc[0].to_dict()

    # Fixed fields
    assert record["name"] == expected_name
    assert record["source"] == expected_source
    assert record["url"] == expected_url
    assert record["version"] == expected_version

    # Accessed date should equal today's date
    assert record["accessed"] == date.today().isoformat()


@pytest.mark.parametrize(
    "key, expected_uuid",
    [
        # Normal case: deterministic UUID
        (
            "GCF_000008665.1",
            f"CDM:{uuid.uuid5(CDM_NAMESPACE, 'GCF_000008665.1')}",
        ),

        # Leading/trailing whitespace should not affect result
        (
            "  GCF_000008665.1  ",
            f"CDM:{uuid.uuid5(CDM_NAMESPACE, 'GCF_000008665.1')}",
        ),

        # Empty string should still return a valid UUID
        (
            "",
            f"CDM:{uuid.uuid5(CDM_NAMESPACE, '')}",
        ),

        # Another different key must yield a different UUID
        (
            "ASM12345v1",
            f"CDM:{uuid.uuid5(CDM_NAMESPACE, 'ASM12345v1')}",
        ),
    ]
)

def test_build_entity_id(key, expected_uuid):
    result = build_entity_id(key)
    assert result == expected_uuid
    assert result.startswith("CDM:")
    assert len(result) > 10


@pytest.mark.parametrize(
    "key, created_date, entity_type, data_source",
    [
        ("GCF_000008665.1", "2000-12-01", "contig_collection", "RefSeq"),
        ("ASM12345v1", None, "genome", "CustomSource"),
    ],
)
def test_build_cdm_entity(key, created_date, entity_type, data_source):
    # Run function
    df, entity_id = build_cdm_entity(
        key_for_uuid=key,
        created_date=created_date,
        entity_type=entity_type,
        data_source=data_source)

    # --- General checks ---
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 1  # exactly one row
    assert entity_id == build_entity_id(key)  # UUID must match

    # --- Row content checks ---
    row = df.iloc[0]
    assert row["entity_id"] == entity_id
    assert row["entity_type"] == entity_type
    assert row["data_source"] == data_source

    # If created_date was provided, must match; else today
    expected_created = created_date or date.today().isoformat()
    assert row["created"] == expected_created

    # updated must be a valid ISO8601 timestamp (seconds precision)
    try:
        datetime.fromisoformat(row["updated"])
    except ValueError:
        pytest.fail(f"Invalid updated timestamp: {row['updated']}")


@pytest.mark.parametrize(
    "entity_id, taxid, collection_type, expected",
    [
        (
            "CDM:1234", "224325", "isolate",
            {"collection_id": "CDM:1234", "contig_collection_type": "isolate", "ncbi_taxon_id": "NCBITaxon:224325", "gtdb_taxon_id": None}
        ),
        (
            "CDM:5678", None, "metagenome",
            {"collection_id": "CDM:5678", "contig_collection_type": "metagenome", "ncbi_taxon_id": None, "gtdb_taxon_id": None}
        ),
    ]
)

def test_build_cdm_contig_collection(entity_id, taxid, collection_type, expected):
    df = build_cdm_contig_collection(entity_id, taxid, collection_type)
    assert len(df) == 1
    row = df.iloc[0].to_dict()
    assert row == expected


@pytest.mark.parametrize(
    "entity_id, rep, expected_rows",
    [
        (
            "CDM:1234",
            {"organism": {"organismName": "Haloarcula marismortui"},
             "assemblyInfo": {"assemblyName": "ASM1234v1"}},
            [
                {
                    "entity_id": "CDM:1234",
                    "name": "Haloarcula marismortui",
                    "description": "RefSeq organism name",
                    "source": "RefSeq"
                },
                {
                    "entity_id": "CDM:1234",
                    "name": "ASM1234v1",
                    "description": "RefSeq assembly name",
                    "source": "RefSeq"
                }
            ]
        ),
        (
            "CDM:5678",
            {"organism": {"scientificName": "Methanocaldococcus jannaschii"}},
            [
                {
                    "entity_id": "CDM:5678",
                    "name": "Methanocaldococcus jannaschii",
                    "description": "RefSeq organism name",
                    "source": "RefSeq"
                }
            ]
        ),
        (
            "CDM:9012",
            {"assemblyInfo": {"assemblyName": "ASM9012v1"}},
            [
                {
                    "entity_id": "CDM:9012",
                    "name": "ASM9012v1",
                    "description": "RefSeq assembly name",
                    "source": "RefSeq"
                }
            ]
        ),
        (
            "CDM:0000",
            {},
            []
        ),
    ]
)

def test_build_cdm_name_rows(entity_id, rep, expected_rows):
    df = build_cdm_name_rows(entity_id, rep)
    result = df.to_dict(orient="records")
    assert result == expected_rows


ENTITY_ID = "CDM:12345"

@pytest.mark.parametrize(
    "rep, request_taxid, expected_identifiers",
    [
        # ---- BioSample ----
        (
            {"assemblyInfo": {"biosample": {"accession": "SAMN123"}}},
            None,
            ["Biosample:SAMN123"],
        ),

        # ---- BioProject ----
        (
            {"assemblyInfo": {"bioproject": {"accession": "PRJNA456"}}},
            None,
            ["BioProject:PRJNA456"],
        ),

        # ---- Taxon from rep ----
        (
            {"organism": {"taxId": "789"}},
            None,
            ["NCBITaxon:789"],
        ),

        # ---- Taxon from request_taxid fallback ----
        (
            {},
            "2468",
            ["NCBITaxon:2468"],
        ),

        # ---- Assembly Accessions (GCF + GCA) ----
        (
            {"accession": "GCF_000001.1"},
            None,
            ["ncbi.assembly:GCF_000001.1"],
        ),
        (
            {"accession": "GCA_000002.1"},
            None,
            ["insdc.gca:GCA_000002.1"],
        ),

        # ---- Deduplication (same identifier twice) ----
        (
            {"assemblyInfo": {"biosample": {"accession": "SAMN999"}}, 
             "assembly": {"biosample": {"biosampleAccession": "SAMN999"}}},
            None,
            ["Biosample:SAMN999"],
        ),
    ]
)

def test_build_cdm_identifier_rows(rep, request_taxid, expected_identifiers):
    rows = build_cdm_identifier_rows(ENTITY_ID, rep, request_taxid)
    identifiers = [r["identifier"] for r in rows]
    assert identifiers == expected_identifiers
    for r in rows:  
        assert r["entity_id"] == ENTITY_ID
        assert r["source"] == "RefSeq"


# ---- dummy data ----
dummy_rep = {
    "accession": "GCF_000001",
    "assemblyInfo": {"assemblyName": "ASM1", "sourceDatabase": "SOURCE_DATABASE_REFSEQ"},
    "organism": {"organismName": "Testus organism", "taxId": 1234},
    "biosample": {"accession": "SAMN0001"},
    "bioproject": {"accession": "PRJNA0001"}
    }


def test_process_report_new_entity():
    seen = set()
    e, c, n, i = process_report(dummy_rep, "1234", seen, debug=False, allow_genbank_date=False)
    # entity dataframe should have 1 row
    assert len(e) == 1
    assert not e[0].empty
    # collections should include taxid
    assert any("NCBITaxon:1234" in str(row) for row in c[0].to_dict("records"))
    # names contain organism and assembly
    assert any("Testus organism" in row["name"] for row in n)
    # identifiers contain accession
    assert any("SAMN0001" in row["identifier"] for row in i)


def test_process_report_duplicate():
    seen = {"GCF_000001"}  # pre-fill key
    e, c, n, i = process_report(dummy_rep, "1234", seen, debug=False, allow_genbank_date=False)
    # should skip
    assert e == [] and c == [] and n == [] and i == []


@patch("refseq_api.fetch_reports_by_taxon")
def test_process_taxon_unique(mock_fetch):
    rep1 = dict(dummy_rep)
    rep2 = dict(dummy_rep); rep2["accession"] = "GCF_000002"
    mock_fetch.return_value = [rep1, rep2]

    seen = set()
    e, c, n, i = process_taxon("1234", api_key=None, debug=False, allow_genbank_date=False, unique_per_taxon=True, seen=seen)
    assert len(e) == 1  # only one kept


def test_finalize_tables_dedup():
    df1 = pd.DataFrame([{"entity_id": "1"}])
    df2 = pd.DataFrame([{"entity_id": "1"}])  # duplicate
    pdf_entity, _, _, _ = finalize_tables([df1, df2], [], [], [])
    assert len(pdf_entity) == 1


@patch("refseq_api.write_delta")
@patch("refseq_api.preview_or_skip")

def test_write_and_preview(mock_preview, mock_write):
    spark = MagicMock()
    pdf = pd.DataFrame([{"entity_id": "1"}])
    write_and_preview(spark, "db", "overwrite", pdf, pdf, pdf, pdf)
    assert mock_write.call_count == 4
    assert mock_preview.call_count == 5


@patch("refseq_api.build_spark")
@patch("refseq_api.write_delta")
@patch("refseq_api.fetch_reports_by_taxon")

def test_main_orchestration(mock_fetch, mock_write, mock_build):
    mock_build.return_value = MagicMock()
    mock_fetch.return_value = [dummy_rep]
    main("1234", None, "db", "overwrite", debug=False)
    assert mock_write.call_count >= 1

