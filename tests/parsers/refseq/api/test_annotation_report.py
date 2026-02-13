### uv run pytest tests/parsers/refseq/api/test_annotation_report.py

import json
from pathlib import Path
from typing import Any
import pytest
from pyspark.sql import SparkSession

# from pyspark.testing import assertDataFrameEqual, assertSchemaEqual

from tests.helpers import assertDataFrameEqual
from pyspark.testing import assertSchemaEqual

from cdm_data_loader_utils.model.kbase_cdm_schema import CDM_SCHEMA
from cdm_data_loader_utils.parsers.refseq.api.annotation_report import (
    apply_prefix,
    load_contig_collection_x_feature,
    load_contig_collection_x_protein,
    load_contig_x_contig_collection,
    load_contigs,
    load_feature_records,
    load_feature_x_protein,
    load_identifiers,
    load_names,
    load_protein,
    load_contig_collections,
    load_contig_x_protein,
    load_contig_x_feature,
    parse_annotation_data,
    to_int,
)
from tests.conftest import TEST_NS

CDM_SCHEMA_LC = {k.lower(): v for k, v in CDM_SCHEMA.items()}


@pytest.mark.parametrize(
    ("input_id", "expected"),
    [
        ("123", "ncbigene:123"),
        ("YP_009725307.1", "refseq:YP_009725307.1"),
        ("NC_008187.1", "refseq:NC_008187.1"),
        ("GCF_000001405.39", "insdc.gcf:GCF_000001405.39"),
        ("random", "random"),
    ],
)
def test_apply_prefix(input_id: str, expected: str) -> None:
    assert apply_prefix(input_id) == expected


@pytest.mark.parametrize(("val", "expected"), [("123", 123), ("abc", None), ("", None)])
def test_to_int(val: str, expected: int | None) -> None:
    assert to_int(val) == expected


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "assembly_accession": "GCF_000869125",
                            "genomic_regions": [{"gene_range": {"accession_version": "NC_008187.1"}}],
                            "gene_id": "1234",
                            "proteins": [{"accession_version": "YP_654573.1"}],
                        }
                    }
                ]
            },
            [
                (
                    "insdc.gcf:GCF_000869125",
                    "insdc.gcf:GCF_000869125",
                    "RefSeq genome ID",
                    "RefSeq",
                    None,
                ),
                (
                    "refseq:NC_008187.1",
                    "refseq:NC_008187.1",
                    "RefSeq assembly ID",
                    "RefSeq",
                    None,
                ),
                (
                    "ncbigene:1234",
                    "ncbigene:1234",
                    "NCBI gene ID",
                    "RefSeq",
                    None,
                ),
                (
                    "refseq:YP_654573.1",
                    "refseq:YP_654573.1",
                    "RefSeq protein ID",
                    "RefSeq",
                    None,
                ),
            ],
        ),
        (
            {"reports": [{"annotation": {"gene_id": "5678"}}]},
            [
                (
                    "ncbigene:5678",
                    "ncbigene:5678",
                    "NCBI gene ID",
                    "RefSeq",
                    None,
                )
            ],
        ),
        ({"reports": [{"annotation": {"name": "no gene id here"}}]}, []),
    ],
)
def test_load_identifiers(input_data: dict[str, Any], expected_output: list[tuple]) -> None:
    result = load_identifiers(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: all name fields present
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "1234",
                            "symbol": "abc",
                            "name": "ABC protein",
                            "locus_tag": "LTG_1234",
                        }
                    }
                ]
            },
            [
                ("ncbigene:1234", "abc", "RefSeq gene symbol", "RefSeq"),
                ("ncbigene:1234", "ABC protein", "RefSeq gene name", "RefSeq"),
                ("ncbigene:1234", "LTG_1234", "RefSeq locus tag", "RefSeq"),
            ],
        ),
        # Case 2: only gene_name present
        (
            {"reports": [{"annotation": {"gene_id": "5678", "name": "Hypothetical protein"}}]},
            [
                (
                    "ncbigene:5678",
                    "Hypothetical protein",
                    "RefSeq gene name",
                    "RefSeq",
                )
            ],
        ),
        # Case 3: no gene_id
        (
            {"reports": [{"annotation": {"name": "Unnamed", "symbol": "XYZ"}}]},
            [],
        ),
        # Case 4: only locus_tag present
        (
            {"reports": [{"annotation": {"gene_id": "8888", "locus_tag": "LTG_8888"}}]},
            [("ncbigene:8888", "LTG_8888", "RefSeq locus tag", "RefSeq")],
        ),
        # Case 5: multiple reports
        (
            {
                "reports": [
                    {"annotation": {"gene_id": "1001", "symbol": "DEF"}},
                    {"annotation": {"gene_id": "1002", "name": "DEF protein"}},
                ]
            },
            [
                ("ncbigene:1001", "DEF", "RefSeq gene symbol", "RefSeq"),
                ("ncbigene:1002", "DEF protein", "RefSeq gene name", "RefSeq"),
            ],
        ),
    ],
)
def test_load_names(input_data: dict[str, Any], expected_output: list[tuple]) -> None:
    result = load_names(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: basic valid input with plus strand
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "1234",
                            "genomic_regions": [
                                {
                                    "gene_range": {
                                        "range": [
                                            {
                                                "begin": "100",
                                                "end": "200",
                                                "orientation": "plus",
                                            }
                                        ]
                                    }
                                }
                            ],
                        }
                    }
                ]
            },
            [
                (
                    "ncbigene:1234",
                    None,
                    None,
                    None,
                    200,
                    None,
                    100,
                    "positive",
                    "ncbigene",
                    None,
                    "unknown",
                )
            ],
        ),
        # Case 2: multiple ranges, different strands
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "5678",
                            "genomic_regions": [
                                {
                                    "gene_range": {
                                        "range": [
                                            {
                                                "begin": "300",
                                                "end": "500",
                                                "orientation": "minus",
                                            },
                                            {
                                                "begin": "600",
                                                "end": "800",
                                                "orientation": "plus",
                                            },
                                        ]
                                    }
                                }
                            ],
                        }
                    }
                ]
            },
            [
                (
                    "ncbigene:5678",
                    None,
                    None,
                    None,
                    500,
                    None,
                    300,
                    "negative",
                    "ncbigene",
                    None,
                    "unknown",
                ),
                (
                    "ncbigene:5678",
                    None,
                    None,
                    None,
                    800,
                    None,
                    600,
                    "positive",
                    "ncbigene",
                    None,
                    "unknown",
                ),
            ],
        ),
        # Case 3: missing orientation
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "9999",
                            "genomic_regions": [{"gene_range": {"range": [{"begin": "1", "end": "2"}]}}],
                        }
                    }
                ]
            },
            [
                (
                    "ncbigene:9999",
                    None,
                    None,
                    None,
                    2,
                    None,
                    1,
                    "unknown",
                    "ncbigene",
                    None,
                    "unknown",
                )
            ],
        ),
        # Case 4: no gene_id
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [
                                {
                                    "gene_range": {
                                        "range": [
                                            {
                                                "begin": "100",
                                                "end": "200",
                                                "orientation": "plus",
                                            }
                                        ]
                                    }
                                }
                            ]
                        }
                    }
                ]
            },
            [],
        ),
        # Case 5: non-integer start/end
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "1111",
                            "genomic_regions": [
                                {
                                    "gene_range": {
                                        "range": [
                                            {
                                                "begin": "abc",
                                                "end": "xyz",
                                                "orientation": "plus",
                                            }
                                        ]
                                    }
                                }
                            ],
                        }
                    }
                ]
            },
            [
                (
                    "ncbigene:1111",
                    None,
                    None,
                    None,
                    None,
                    None,
                    None,
                    "positive",
                    "ncbigene",
                    None,
                    "unknown",
                )
            ],
        ),
    ],
)
def test_load_feature_records(input_data: dict[str, Any], expected_output: list[tuple]):
    result = load_feature_records(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: valid mapping with assembly accession
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "12345",
                            "annotations": [{"assembly_accession": "GCF_000001.1"}],
                        }
                    }
                ]
            },
            [("insdc.gcf:GCF_000001.1", "ncbigene:12345")],
        ),
        # Case 2: missing gene_id
        (
            {"reports": [{"annotation": {"annotations": [{"assembly_accession": "GCF_000002.1"}]}}]},
            [],
        ),
        # Case 3: missing annotations field
        (
            {"reports": [{"annotation": {"gene_id": "67890"}}]},
            [],
        ),
        # Case 4: annotations list is empty
        (
            {"reports": [{"annotation": {"gene_id": "99999", "annotations": []}}]},
            [],
        ),
        # Case 5: missing assembly_accession in first annotations element
        (
            {"reports": [{"annotation": {"gene_id": "13579", "annotations": [{}]}}]},
            [],
        ),
    ],
)
def test_load_contig_collection_x_feature(input_data: dict[str, Any], expected_output: list[tuple[str, str]]) -> None:
    result = load_contig_collection_x_feature(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: Valid report with multiple proteins
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "annotations": [{"assembly_accession": "GCF_000001"}],
                            "proteins": [
                                {"accession_version": "XP_123"},
                                {"accession_version": "XP_456"},
                            ],
                        }
                    }
                ]
            },
            [
                ("insdc.gcf:GCF_000001", "refseq:XP_123"),
                ("insdc.gcf:GCF_000001", "refseq:XP_456"),
            ],
        ),
        # Case 2: No proteins
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "annotations": [{"assembly_accession": "GCF_000002"}],
                            "proteins": [],
                        }
                    }
                ]
            },
            [],
        ),
        # Case 3: No annotations
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "proteins": [{"accession_version": "XP_789"}],
                            # 'annotations' key missing
                        }
                    }
                ]
            },
            [],
        ),
        # Case 4: Missing assembly_accession in annotations[0]
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "annotations": [{}],
                            "proteins": [{"accession_version": "XP_789"}],
                        }
                    }
                ]
            },
            [],
        ),
        # Case 5: Some proteins missing accession_version
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "annotations": [{"assembly_accession": "GCF_000003"}],
                            "proteins": [
                                {"accession_version": "XP_111"},
                                {},
                                {"accession_version": "XP_222"},
                            ],
                        }
                    }
                ]
            },
            [
                ("insdc.gcf:GCF_000003", "refseq:XP_111"),
                ("insdc.gcf:GCF_000003", "refseq:XP_222"),
            ],
        ),
    ],
)
def test_load_contig_collection_x_protein(input_data: dict[str, Any], expected_output: list[tuple]) -> None:
    result = load_contig_collection_x_protein(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: valid gene with multiple proteins
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "4156311",
                            "proteins": [
                                {"accession_version": "XP_001"},
                                {"accession_version": "XP_002"},
                            ],
                        }
                    }
                ]
            },
            [
                ("ncbigene:4156311", "refseq:XP_001"),
                ("ncbigene:4156311", "refseq:XP_002"),
            ],
        ),
        # Case 2: no gene_id
        (
            {"reports": [{"annotation": {"proteins": [{"accession_version": "XP_999"}]}}]},
            [],
        ),
        # Case 3: gene with no proteins
        (
            {"reports": [{"annotation": {"gene_id": "4156312"}}]},
            [],
        ),
        # Case 4: some proteins missing accession_version
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "4156313",
                            "proteins": [
                                {"accession_version": "XP_777"},
                                {},
                                {"accession_version": "XP_888"},
                            ],
                        }
                    }
                ]
            },
            [
                ("ncbigene:4156313", "refseq:XP_777"),
                ("ncbigene:4156313", "refseq:XP_888"),
            ],
        ),
        # Case 5: empty report list
        ({"reports": []}, []),
    ],
)
def test_load_feature_x_protein(input_data: dict[str, Any], expected_output: list[tuple[str, str]]) -> None:
    result = load_feature_x_protein(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: Valid contig and assembly
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [{"gene_range": {"accession_version": "NC_000001.11"}}],
                            "annotations": [{"assembly_accession": "GCF_000001.1"}],
                        }
                    }
                ]
            },
            [("refseq:NC_000001.11", "insdc.gcf:GCF_000001.1")],
        ),
        # Case 2: Missing genomic_regions
        (
            {"reports": [{"annotation": {"annotations": [{"assembly_accession": "GCF_000002.1"}]}}]},
            [],
        ),
        # Case 3: Missing annotations
        (
            {"reports": [{"annotation": {"genomic_regions": [{"gene_range": {"accession_version": "NC_000003.11"}}]}}]},
            [],
        ),
        # Case 4: Missing accession_version in region
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [{"gene_range": {}}],
                            "annotations": [{"assembly_accession": "GCF_000004.1"}],
                        }
                    }
                ]
            },
            [],
        ),
        # Case 5: Missing assembly_accession in annotations
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [{"gene_range": {"accession_version": "NC_000005.11"}}],
                            "annotations": [{}],
                        }
                    }
                ]
            },
            [],
        ),
        # Case 6: Multiple reports, one valid
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [{"gene_range": {"accession_version": "NC_000006.11"}}],
                            "annotations": [{"assembly_accession": "GCF_000006.1"}],
                        }
                    },
                    {
                        "annotation": {
                            "genomic_regions": [{"gene_range": {"accession_version": "NC_000007.11"}}],
                            "annotations": [{}],
                        }
                    },
                ]
            },
            [("refseq:NC_000006.11", "insdc.gcf:GCF_000006.1")],
        ),
    ],
)
def test_load_contig_x_contig_collection(input_data: dict[str, Any], expected_output: list[tuple[str, str]]) -> None:
    result = load_contig_x_contig_collection(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: Valid contig with accession_version
        (
            {"reports": [{"annotation": {"genomic_regions": [{"gene_range": {"accession_version": "NC_000001.11"}}]}}]},
            [("refseq:NC_000001.11", None, None, None)],
        ),
        # Case 2: Multiple contigs, different accession_versions
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000001.11"}},
                                {"gene_range": {"accession_version": "NC_000002.12"}},
                            ]
                        }
                    }
                ]
            },
            [
                ("refseq:NC_000001.11", None, None, None),
                ("refseq:NC_000002.12", None, None, None),
            ],
        ),
        # Case 3: Duplicate accession versions
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000003.13"}},
                                {"gene_range": {"accession_version": "NC_000003.13"}},
                            ]
                        }
                    }
                ]
            },
            [("refseq:NC_000003.13", None, None, None)],
        ),
        # Case 4: Missing accession_version
        (
            {"reports": [{"annotation": {"genomic_regions": [{"gene_range": {}}]}}]},
            [],
        ),
        # Case 5: Empty reports
        (
            {"reports": []},
            [],
        ),
    ],
)
def test_load_contigs(input_data: dict[str, Any], expected_output: list[tuple]) -> None:
    result = load_contigs(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        # Case 1: Valid input with one contig
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "12345",
                            "genomic_regions": [{"gene_range": {"accession_version": "NC_000001.11"}}],
                        }
                    }
                ]
            },
            [("refseq:NC_000001.11", "ncbigene:12345")],
        ),
        # Case 2: Multiple contigs for the same gene
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "67890",
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000002.11"}},
                                {"gene_range": {"accession_version": "NC_000003.11"}},
                            ],
                        }
                    }
                ]
            },
            [
                ("refseq:NC_000002.11", "ncbigene:67890"),
                ("refseq:NC_000003.11", "ncbigene:67890"),
            ],
        ),
        # Case 3: Missing accession_version in gene_range
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "99999",
                            "genomic_regions": [{"gene_range": {}}],
                        }
                    }
                ]
            },
            [],
        ),
        # Case 4: Missing gene_id in annotation
        (
            {"reports": [{"annotation": {"genomic_regions": [{"gene_range": {"accession_version": "NC_000004.11"}}]}}]},
            [],
        ),
        # Case 5: accession_version starting with GCF_ (tests apply_prefix logic)
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "13579",
                            "genomic_regions": [{"gene_range": {"accession_version": "GCF_000005.1"}}],
                        }
                    }
                ]
            },
            [("insdc.gcf:GCF_000005.1", "ncbigene:13579")],
        ),
    ],
)
def test_load_contig_x_feature(input_data: dict[str, Any], expected_output: list[tuple[str, str]]) -> None:
    result = load_contig_x_feature(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "789",
                            "annotations": [
                                {
                                    "genomic_regions": [{"gene_range": {"accession_version": "NC_000003.11"}}],
                                    "proteins": [],
                                }
                            ],
                        }
                    }
                ]
            },
            [],
        ),
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "999",
                            "annotations": [{"proteins": [{"accession_version": "YP_999999.1"}]}],
                        }
                    }
                ]
            },
            [],
        ),
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "555",
                            "annotations": [
                                {
                                    "genomic_regions": [{"gene_range": {"accession_version": "NC_000004.11"}}],
                                    "proteins": [{}],
                                }
                            ],
                        }
                    }
                ]
            },
            [],
        ),
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "666",
                            "genomic_regions": [{"gene_range": {"accession_version": "NC_000001.11"}}],
                            "proteins": [{"accession_version": "YP_000001.1"}],
                        },
                    }
                ]
            },
            [
                (
                    apply_prefix("NC_000001.11"),
                    apply_prefix("YP_000001.1"),
                )
            ],
        ),
    ],
)
def test_load_contig_x_protein(
    input_data: dict[str, Any],
    expected_output: list[tuple[str, str]],
) -> None:
    result = load_contig_x_protein(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    ("input_data", "expected_ids"),
    [
        # Case 1: Basic valid accession
        (
            {"reports": [{"annotation": {"annotations": [{"assembly_accession": "GCF_000001.1"}]}}]},
            ["insdc.gcf:GCF_000001.1"],
        ),
        # Case 2: Multiple annotations with same accession (deduplicated)
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "annotations": [
                                {"assembly_accession": "GCF_000002.1"},
                                {"assembly_accession": "GCF_000002.1"},
                            ]
                        }
                    }
                ]
            },
            ["insdc.gcf:GCF_000002.1"],
        ),
        # Case 3: Mixed valid and missing accession
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "annotations": [
                                {"assembly_accession": "GCF_000003.1"},
                                {},
                            ]
                        }
                    }
                ]
            },
            ["insdc.gcf:GCF_000003.1"],
        ),
        # Case 4: No `annotations` field
        (
            {"reports": [{"annotation": {}}]},
            [],
        ),
        # Case 5: Empty `reports` list
        (
            {"reports": []},
            [],
        ),
        # Case 6: No `assembly_accession`
        (
            {"reports": [{"annotation": {"annotations": [{}]}}]},
            [],
        ),
    ],
)
def test_load_contig_collections(input_data: dict[str, Any], expected_ids: list[str]) -> None:
    result = load_contig_collections(input_data)

    # Extract contig_collection_id from dict
    result_ids = sorted(r["contig_collection_id"] for r in result)
    assert result_ids == sorted(expected_ids)


@pytest.mark.parametrize(
    ("input_data", "expected_output"),
    [
        (
            {"reports": [{"annotation": {"proteins": [{"name": "Nucleocapsid", "length": 419}]}}]},
            [],
        ),
        (
            {"reports": [{"annotation": {"proteins": []}}]},
            [],
        ),
        (
            {"reports": [{"annotation": {}}]},
            [],
        ),
    ],
)
def test_load_protein(
    input_data: dict[str, Any], expected_output: list[tuple[str, None, str | None, None, int | None, None]]
) -> None:
    result = load_protein(input_data)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.requires_spark
def test_parse_annotation_data(spark: SparkSession, test_data_dir: Path) -> None:
    """
    Test that parse_annotation_data produces expected tables with correct schemas and non-empty output.
    """

    spark.sql(f"CREATE DATABASE IF NOT EXISTS {TEST_NS}")

    # Load and parse test JSON
    input_path = test_data_dir / "refseq" / "annotation_report.json"
    dataset = json.load(input_path.open())

    # Run parser
    parse_annotation_data(spark, [dataset], TEST_NS, use_metastore=True)

    # Load expected schema definitions
    expected_path = test_data_dir / "refseq" / "annotation_report.parsed.json"
    expected_tables = json.load(expected_path.open())

    actual_tables = {table.name: spark.table(f"{TEST_NS}.{table.name}") for table in spark.catalog.listTables(TEST_NS)}

    # Check each produced table
    for table_name, actual_df in actual_tables.items():
        assert table_name in expected_tables, f"Unexpected table produced: {table_name}"

        expected_schema = CDM_SCHEMA_LC[table_name]

        # Schema must match exactly
        expected_df = spark.createDataFrame(
            expected_tables[table_name],
            schema=expected_schema,
        )
        assertSchemaEqual(expected_df.schema, actual_df.schema)

        # Must have at least one row
        assert actual_df.count() > 0

    # Ensure no expected table is missing
    assert set(expected_tables) == set(actual_tables)
