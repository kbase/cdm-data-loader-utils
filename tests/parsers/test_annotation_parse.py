import pytest
import json

from annotation_parse import (
    load_identifiers,
    load_names,
    load_feature_records,
    load_contig_collection_x_feature,
    load_contig_collection_x_protein,
    load_feature_x_protein,
    load_contig_x_contig_collection,
    load_contigs,
)


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "1234",
                            "name": "hypothetical protein",
                            "relationship": "RefSeq gene symbol",
                        }
                    }
                ]
            },
            [
                (
                    "ncbigene:1234",
                    "1234",
                    "hypothetical protein",
                    "RefSeq",
                    "RefSeq gene symbol",
                )
            ],
        ),
        (
            {"reports": [{"annotation": {"gene_id": "5678", "name": "some protein"}}]},
            [("ncbigene:5678", "5678", "some protein", "RefSeq", None)],
        ),
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "name": "no gene id here",
                            "relationship": "RefSeq locus tag",
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
                            "gene_id": "1001",
                            "name": "abc",
                            "relationship": "RefSeq gene symbol",
                        }
                    },
                    {"annotation": {"gene_id": "1002", "name": "xyz"}},
                ]
            },
            [
                ("ncbigene:1001", "1001", "abc", "RefSeq", "RefSeq gene symbol"),
                ("ncbigene:1002", "1002", "xyz", "RefSeq", None),
            ],
        ),
    ],
)
def test_load_identifiers(tmp_path, input_data, expected_output):
    input_file = tmp_path / "test.json"
    input_file.write_text(json.dumps(input_data))

    result = load_identifiers(input_file)
    assert result == expected_output


@pytest.mark.parametrize(
    "input_data, expected_output",
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
            {
                "reports": [
                    {"annotation": {"gene_id": "5678", "name": "Hypothetical protein"}}
                ]
            },
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
def test_load_names(tmp_path, input_data, expected_output):
    input_file = tmp_path / "test.json"
    input_file.write_text(json.dumps(input_data))

    result = load_names(input_file)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    "input_data, expected_output",
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
                    "RefSeq",
                    None,
                    "gene",
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
                    "RefSeq",
                    None,
                    "gene",
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
                    "RefSeq",
                    None,
                    "gene",
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
                            "genomic_regions": [
                                {"gene_range": {"range": [{"begin": "1", "end": "2"}]}}
                            ],
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
                    "RefSeq",
                    None,
                    "gene",
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
                    "RefSeq",
                    None,
                    "gene",
                )
            ],
        ),
    ],
)
def test_load_feature_records(tmp_path, input_data, expected_output):
    input_file = tmp_path / "features.json"
    input_file.write_text(json.dumps(input_data))

    result = load_feature_records(input_file)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        # Case 1: valid mapping
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "12345",
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000001.11"}}
                            ],
                        }
                    }
                ]
            },
            [("refseq:NC_000001.11", "ncbigene:12345")],
        ),
        # Case 2: no gene_id
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000002.11"}}
                            ]
                        }
                    }
                ]
            },
            [],
        ),
        # Case 3: no genomic_regions
        (
            {"reports": [{"annotation": {"gene_id": "67890"}}]},
            [],
        ),
        # Case 4: empty genomic_regions list
        (
            {"reports": [{"annotation": {"gene_id": "99999", "genomic_regions": []}}]},
            [],
        ),
        # Case 5: missing accession_version
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "gene_id": "13579",
                            "genomic_regions": [{"gene_range": {}}],
                        }
                    }
                ]
            },
            [],
        ),
    ],
)
def test_load_contig_collection_x_feature(tmp_path, input_data, expected_output):
    input_file = tmp_path / "contig_feature.json"
    input_file.write_text(json.dumps(input_data))

    result = load_contig_collection_x_feature(input_file)
    assert result == expected_output


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        # Case 1: Valid report with multiple proteins
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "proteins": [
                                {"accession_version": "XP_123"},
                                {"accession_version": "XP_456"},
                            ],
                            "annotations": [{"assembly_accession": "GCF_000001"}],
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
                            "proteins": [],
                            "annotations": [{"assembly_accession": "GCF_000002"}],
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
                    {"annotation": {"proteins": [{"accession_version": "XP_789"}]}}
                ]
            },
            [],
        ),
        # Case 4: Missing assembly_accession
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "proteins": [{"accession_version": "XP_789"}],
                            "annotations": [{}],
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
                            "proteins": [
                                {"accession_version": "XP_111"},
                                {},
                                {"accession_version": "XP_222"},
                            ],
                            "annotations": [{"assembly_accession": "GCF_000003"}],
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
def test_load_contig_collection_x_protein(tmp_path, input_data, expected_output):
    input_file = tmp_path / "protein_links.json"
    input_file.write_text(json.dumps(input_data))

    result = load_contig_collection_x_protein(input_file)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    "input_data, expected_output",
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
            {
                "reports": [
                    {"annotation": {"proteins": [{"accession_version": "XP_999"}]}}
                ]
            },
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
def test_load_feature_x_protein(tmp_path, input_data, expected_output):
    input_file = tmp_path / "feature_protein.json"
    input_file.write_text(json.dumps(input_data))

    result = load_feature_x_protein(input_file)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        # Case 1: Valid contig and assembly
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000001.11"}}
                            ],
                            "annotations": [{"assembly_accession": "GCF_000001.1"}],
                        }
                    }
                ]
            },
            [("refseq:NC_000001.11", "insdc.gcf:GCF_000001.1")],
        ),
        # Case 2: Missing genomic_regions
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "annotations": [{"assembly_accession": "GCF_000002.1"}]
                        }
                    }
                ]
            },
            [],
        ),
        # Case 3: Missing annotations
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000003.11"}}
                            ]
                        }
                    }
                ]
            },
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
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000005.11"}}
                            ],
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
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000006.11"}}
                            ],
                            "annotations": [{"assembly_accession": "GCF_000006.1"}],
                        }
                    },
                    {
                        "annotation": {
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000007.11"}}
                            ],
                            "annotations": [{}],
                        }
                    },
                ]
            },
            [("refseq:NC_000006.11", "insdc.gcf:GCF_000006.1")],
        ),
    ],
)
def test_load_contig_x_contig_collection(tmp_path, input_data, expected_output):
    input_file = tmp_path / "contig_collection.json"
    input_file.write_text(json.dumps(input_data))

    result = load_contig_x_contig_collection(input_file)
    assert sorted(result) == sorted(expected_output)


@pytest.mark.parametrize(
    "input_data, expected_output",
    [
        # Case 1: Valid contig with accession_version
        (
            {
                "reports": [
                    {
                        "annotation": {
                            "genomic_regions": [
                                {"gene_range": {"accession_version": "NC_000001.11"}}
                            ]
                        }
                    }
                ]
            },
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
def test_load_contigs(tmp_path, input_data, expected_output):
    input_file = tmp_path / "contig.json"
    input_file.write_text(json.dumps(input_data, indent=2))

    result = load_contigs(input_file)
    assert sorted(result) == sorted(expected_output)
