"""Expected PyArrow schema for pubmed_book_article (xml_tag='PubmedBookArticle')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field(
            "BookDocument",
            pa.struct(
                [
                    pa.field(
                        "PMID",
                        pa.struct(
                            [
                                pa.field("_Version", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                                pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=False,
                    ),
                    pa.field(
                        "ArticleIdList",
                        pa.struct(
                            [
                                pa.field(
                                    "ArticleId",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_IdType",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "#text",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=False,
                    ),
                    pa.field(
                        "Book",
                        pa.struct(
                            [
                                pa.field(
                                    "Publisher",
                                    pa.struct(
                                        [
                                            pa.field(
                                                "PublisherName",
                                                pa.string(),
                                                nullable=False,
                                                metadata={"xsd_type": "xs:string"},
                                            ),
                                            pa.field(
                                                "PublisherLocation",
                                                pa.string(),
                                                nullable=True,
                                                metadata={"xsd_type": "xs:string"},
                                            ),
                                        ]
                                    ),
                                    nullable=False,
                                ),
                                pa.field(
                                    "BookTitle",
                                    pa.struct(
                                        [
                                            pa.field(
                                                "_book", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "_part", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "_sec", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                        ]
                                    ),
                                    nullable=False,
                                ),
                                pa.field(
                                    "PubDate",
                                    pa.struct(
                                        [
                                            pa.field(
                                                "Year", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Month", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Day", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Season", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "MedlineDate",
                                                pa.string(),
                                                nullable=True,
                                                metadata={"xsd_type": "xs:string"},
                                            ),
                                        ]
                                    ),
                                    nullable=False,
                                ),
                                pa.field(
                                    "BeginningDate",
                                    pa.struct(
                                        [
                                            pa.field(
                                                "Year", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Month", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Day", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Season", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                        ]
                                    ),
                                    nullable=True,
                                ),
                                pa.field(
                                    "EndingDate",
                                    pa.struct(
                                        [
                                            pa.field(
                                                "Year", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Month", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Day", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "Season", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                        ]
                                    ),
                                    nullable=True,
                                ),
                                pa.field(
                                    "AuthorList",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_CompleteYN",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_Type",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Author",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_ValidYN",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "_EqualContrib",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "LastName",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "ForeName",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "Initials",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "Suffix",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "CollectiveName",
                                                                            pa.struct(
                                                                                [
                                                                                    pa.field(
                                                                                        "_Investigators",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:IDREF"
                                                                                        },
                                                                                    ),
                                                                                    pa.field(
                                                                                        "#text",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                ]
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                        pa.field(
                                                                            "Identifier",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.struct(
                                                                                        [
                                                                                            pa.field(
                                                                                                "_Source",
                                                                                                pa.string(),
                                                                                                nullable=False,
                                                                                                metadata={
                                                                                                    "xsd_type": "xs:string"
                                                                                                },
                                                                                            ),
                                                                                            pa.field(
                                                                                                "#text",
                                                                                                pa.string(),
                                                                                                nullable=True,
                                                                                                metadata={
                                                                                                    "xsd_type": "xs:string"
                                                                                                },
                                                                                            ),
                                                                                        ]
                                                                                    ),
                                                                                    nullable=True,
                                                                                ),
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                        pa.field(
                                                                            "AffiliationInfo",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.struct(
                                                                                        [
                                                                                            pa.field(
                                                                                                "Affiliation",
                                                                                                pa.string(),
                                                                                                nullable=False,
                                                                                                metadata={
                                                                                                    "xsd_type": "xs:string"
                                                                                                },
                                                                                            ),
                                                                                            pa.field(
                                                                                                "Identifier",
                                                                                                pa.list_(
                                                                                                    pa.field(
                                                                                                        "item",
                                                                                                        pa.struct(
                                                                                                            [
                                                                                                                pa.field(
                                                                                                                    "_Source",
                                                                                                                    pa.string(),
                                                                                                                    nullable=False,
                                                                                                                    metadata={
                                                                                                                        "xsd_type": "xs:string"
                                                                                                                    },
                                                                                                                ),
                                                                                                                pa.field(
                                                                                                                    "#text",
                                                                                                                    pa.string(),
                                                                                                                    nullable=True,
                                                                                                                    metadata={
                                                                                                                        "xsd_type": "xs:string"
                                                                                                                    },
                                                                                                                ),
                                                                                                            ]
                                                                                                        ),
                                                                                                        nullable=True,
                                                                                                    ),
                                                                                                ),
                                                                                                nullable=True,
                                                                                            ),
                                                                                        ]
                                                                                    ),
                                                                                    nullable=True,
                                                                                ),
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ),
                                                        nullable=False,
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=True,
                                ),
                                pa.field(
                                    "InvestigatorList",
                                    pa.struct(
                                        [
                                            pa.field("_ID", pa.string(), nullable=True, metadata={"xsd_type": "xs:ID"}),
                                            pa.field(
                                                "Investigator",
                                                pa.list_(
                                                    pa.field(
                                                        "item",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_ValidYN",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "LastName",
                                                                    pa.string(),
                                                                    nullable=False,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "ForeName",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "Initials",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "Suffix",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "Identifier",
                                                                    pa.list_(
                                                                        pa.field(
                                                                            "item",
                                                                            pa.struct(
                                                                                [
                                                                                    pa.field(
                                                                                        "_Source",
                                                                                        pa.string(),
                                                                                        nullable=False,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                    pa.field(
                                                                                        "#text",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                ]
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                    ),
                                                                    nullable=True,
                                                                ),
                                                                pa.field(
                                                                    "AffiliationInfo",
                                                                    pa.list_(
                                                                        pa.field(
                                                                            "item",
                                                                            pa.struct(
                                                                                [
                                                                                    pa.field(
                                                                                        "Affiliation",
                                                                                        pa.string(),
                                                                                        nullable=False,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                    pa.field(
                                                                                        "Identifier",
                                                                                        pa.list_(
                                                                                            pa.field(
                                                                                                "item",
                                                                                                pa.struct(
                                                                                                    [
                                                                                                        pa.field(
                                                                                                            "_Source",
                                                                                                            pa.string(),
                                                                                                            nullable=False,
                                                                                                            metadata={
                                                                                                                "xsd_type": "xs:string"
                                                                                                            },
                                                                                                        ),
                                                                                                        pa.field(
                                                                                                            "#text",
                                                                                                            pa.string(),
                                                                                                            nullable=True,
                                                                                                            metadata={
                                                                                                                "xsd_type": "xs:string"
                                                                                                            },
                                                                                                        ),
                                                                                                    ]
                                                                                                ),
                                                                                                nullable=True,
                                                                                            ),
                                                                                        ),
                                                                                        nullable=True,
                                                                                    ),
                                                                                ]
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                    ),
                                                                    nullable=True,
                                                                ),
                                                            ]
                                                        ),
                                                        nullable=True,
                                                    ),
                                                ),
                                                nullable=False,
                                            ),
                                        ]
                                    ),
                                    nullable=True,
                                ),
                                pa.field("Volume", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("VolumeTitle", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("Edition", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field(
                                    "CollectionTitle",
                                    pa.struct(
                                        [
                                            pa.field(
                                                "_book", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "_part", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "_sec", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                            pa.field(
                                                "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                        ]
                                    ),
                                    nullable=True,
                                ),
                                pa.field(
                                    "Isbn",
                                    pa.list_(
                                        pa.field(
                                            "item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ),
                                    nullable=True,
                                ),
                                pa.field(
                                    "ELocationID",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_EIdType",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_ValidYN",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "#text",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=True,
                                ),
                                pa.field("Medium", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field(
                                    "ReportNumber", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                ),
                            ]
                        ),
                        nullable=False,
                    ),
                    pa.field(
                        "LocationLabel",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_Type", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ]
                                ),
                                nullable=True,
                            ),
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "ArticleTitle",
                        pa.struct(
                            [
                                pa.field("_book", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("_part", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("_sec", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field("VernacularTitle", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "Pagination",
                        pa.struct(
                            [
                                pa.field("StartPage", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("EndPage", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("MedlinePgn", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "Language",
                        pa.list_(
                            pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "AuthorList",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_CompleteYN",
                                            pa.string(),
                                            nullable=True,
                                            metadata={"xsd_type": "xs:string"},
                                        ),
                                        pa.field(
                                            "_Type", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "Author",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_ValidYN",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "_EqualContrib",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "LastName",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "ForeName",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "Initials",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "Suffix",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "CollectiveName",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_Investigators",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:IDREF"},
                                                                        ),
                                                                        pa.field(
                                                                            "#text",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "Identifier",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_Source",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "#text",
                                                                                    pa.string(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                            ]
                                                                        ),
                                                                        nullable=True,
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "AffiliationInfo",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "Affiliation",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "Identifier",
                                                                                    pa.list_(
                                                                                        pa.field(
                                                                                            "item",
                                                                                            pa.struct(
                                                                                                [
                                                                                                    pa.field(
                                                                                                        "_Source",
                                                                                                        pa.string(),
                                                                                                        nullable=False,
                                                                                                        metadata={
                                                                                                            "xsd_type": "xs:string"
                                                                                                        },
                                                                                                    ),
                                                                                                    pa.field(
                                                                                                        "#text",
                                                                                                        pa.string(),
                                                                                                        nullable=True,
                                                                                                        metadata={
                                                                                                            "xsd_type": "xs:string"
                                                                                                        },
                                                                                                    ),
                                                                                                ]
                                                                                            ),
                                                                                            nullable=True,
                                                                                        ),
                                                                                    ),
                                                                                    nullable=True,
                                                                                ),
                                                                            ]
                                                                        ),
                                                                        nullable=True,
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ]
                                                    ),
                                                    nullable=True,
                                                ),
                                            ),
                                            nullable=False,
                                        ),
                                    ]
                                ),
                                nullable=True,
                            ),
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "InvestigatorList",
                        pa.struct(
                            [
                                pa.field("_ID", pa.string(), nullable=True, metadata={"xsd_type": "xs:ID"}),
                                pa.field(
                                    "Investigator",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_ValidYN",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "LastName",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "ForeName",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Initials",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Suffix",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Identifier",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_Source",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "#text",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                    pa.field(
                                                        "AffiliationInfo",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "Affiliation",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "Identifier",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.struct(
                                                                                        [
                                                                                            pa.field(
                                                                                                "_Source",
                                                                                                pa.string(),
                                                                                                nullable=False,
                                                                                                metadata={
                                                                                                    "xsd_type": "xs:string"
                                                                                                },
                                                                                            ),
                                                                                            pa.field(
                                                                                                "#text",
                                                                                                pa.string(),
                                                                                                nullable=True,
                                                                                                metadata={
                                                                                                    "xsd_type": "xs:string"
                                                                                                },
                                                                                            ),
                                                                                        ]
                                                                                    ),
                                                                                    nullable=True,
                                                                                ),
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "PublicationType",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_UI", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ]
                                ),
                                nullable=True,
                            ),
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "Abstract",
                        pa.struct(
                            [
                                pa.field(
                                    "AbstractText",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_Label",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_NlmCategory",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "#text",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                                pa.field(
                                    "CopyrightInformation",
                                    pa.string(),
                                    nullable=True,
                                    metadata={"xsd_type": "xs:string"},
                                ),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "Sections",
                        pa.struct(
                            [
                                pa.field(
                                    "Section",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "LocationLabel",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_Type",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "#text",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                            ]
                                                        ),
                                                        nullable=True,
                                                    ),
                                                    pa.field(
                                                        "SectionTitle",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_book",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "_part",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "_sec",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                                pa.field(
                                                                    "#text",
                                                                    pa.string(),
                                                                    nullable=True,
                                                                    metadata={"xsd_type": "xs:string"},
                                                                ),
                                                            ]
                                                        ),
                                                        nullable=False,
                                                    ),
                                                    pa.field(
                                                        "Section",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "LocationLabel",
                                                                            pa.struct(
                                                                                [
                                                                                    pa.field(
                                                                                        "_Type",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                    pa.field(
                                                                                        "#text",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                ]
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                        pa.field(
                                                                            "SectionTitle",
                                                                            pa.struct(
                                                                                [
                                                                                    pa.field(
                                                                                        "_book",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                    pa.field(
                                                                                        "_part",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                    pa.field(
                                                                                        "_sec",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                    pa.field(
                                                                                        "#text",
                                                                                        pa.string(),
                                                                                        nullable=True,
                                                                                        metadata={
                                                                                            "xsd_type": "xs:string"
                                                                                        },
                                                                                    ),
                                                                                ]
                                                                            ),
                                                                            nullable=False,
                                                                        ),
                                                                        pa.field(
                                                                            "Section",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item", pa.string(), nullable=True
                                                                                ),
                                                                            ),
                                                                            nullable=True,
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "KeywordList",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_Owner", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "Keyword",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_MajorTopicYN",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "#text",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                        ]
                                                    ),
                                                    nullable=True,
                                                ),
                                            ),
                                            nullable=False,
                                        ),
                                    ]
                                ),
                                nullable=True,
                            ),
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "ContributionDate",
                        pa.struct(
                            [
                                pa.field("Year", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                                pa.field("Month", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("Day", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field("Season", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "DateRevised",
                        pa.struct(
                            [
                                pa.field("Year", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                                pa.field("Month", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                                pa.field("Day", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "GrantList",
                        pa.struct(
                            [
                                pa.field("_CompleteYN", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                pa.field(
                                    "Grant",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "GrantID",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Acronym",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Agency",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Country",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "ItemList",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_ListType", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "Item",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.string(),
                                                    nullable=True,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                            ),
                                            nullable=False,
                                        ),
                                    ]
                                ),
                                nullable=True,
                            ),
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "ReferenceList",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "Title", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "Reference",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "Citation",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "ArticleIdList",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "ArticleId",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.struct(
                                                                                        [
                                                                                            pa.field(
                                                                                                "_IdType",
                                                                                                pa.string(),
                                                                                                nullable=True,
                                                                                                metadata={
                                                                                                    "xsd_type": "xs:string"
                                                                                                },
                                                                                            ),
                                                                                            pa.field(
                                                                                                "#text",
                                                                                                pa.string(),
                                                                                                nullable=True,
                                                                                                metadata={
                                                                                                    "xsd_type": "xs:string"
                                                                                                },
                                                                                            ),
                                                                                        ]
                                                                                    ),
                                                                                    nullable=True,
                                                                                ),
                                                                            ),
                                                                            nullable=False,
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ]
                                                    ),
                                                    nullable=True,
                                                ),
                                            ),
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "ReferenceList",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "Title",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "Reference",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "Citation",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "ArticleIdList",
                                                                                    pa.struct(
                                                                                        [
                                                                                            pa.field(
                                                                                                "ArticleId",
                                                                                                pa.list_(
                                                                                                    pa.field(
                                                                                                        "item",
                                                                                                        pa.struct(
                                                                                                            [
                                                                                                                pa.field(
                                                                                                                    "_IdType",
                                                                                                                    pa.string(),
                                                                                                                    nullable=True,
                                                                                                                    metadata={
                                                                                                                        "xsd_type": "xs:string"
                                                                                                                    },
                                                                                                                ),
                                                                                                                pa.field(
                                                                                                                    "#text",
                                                                                                                    pa.string(),
                                                                                                                    nullable=True,
                                                                                                                    metadata={
                                                                                                                        "xsd_type": "xs:string"
                                                                                                                    },
                                                                                                                ),
                                                                                                            ]
                                                                                                        ),
                                                                                                        nullable=True,
                                                                                                    ),
                                                                                                ),
                                                                                                nullable=False,
                                                                                            ),
                                                                                        ]
                                                                                    ),
                                                                                    nullable=True,
                                                                                ),
                                                                            ]
                                                                        ),
                                                                        nullable=True,
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "ReferenceList",
                                                                pa.list_(
                                                                    pa.field("item", pa.string(), nullable=True),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ]
                                                    ),
                                                    nullable=True,
                                                ),
                                            ),
                                            nullable=True,
                                        ),
                                    ]
                                ),
                                nullable=True,
                            ),
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=False,
        ),
        pa.field(
            "PubmedBookData",
            pa.struct(
                [
                    pa.field(
                        "History",
                        pa.struct(
                            [
                                pa.field(
                                    "PubMedPubDate",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_PubStatus",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Year",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Month",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Day",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Hour",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Minute",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Second",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field("PublicationStatus", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "ArticleIdList",
                        pa.struct(
                            [
                                pa.field(
                                    "ArticleId",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_IdType",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "#text",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=False,
                    ),
                    pa.field(
                        "ObjectList",
                        pa.struct(
                            [
                                pa.field(
                                    "Object",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_Type",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "Param",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_Name",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "#text",
                                                                            pa.string(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=True,
        ),
    ]
)
