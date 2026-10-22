"""PyArrow schema for uniprot (xml_tag='entry')."""

import pyarrow as pa

SCHEMA = pa.schema(
    [
        pa.field("_dataset", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field("_created", pa.date32(), nullable=False, metadata={"xsd_type": "xs:date"}),
        pa.field("_modified", pa.date32(), nullable=False, metadata={"xsd_type": "xs:date"}),
        pa.field("_version", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}),
        pa.field(
            "accession",
            pa.list_(
                pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
            ),
            nullable=False,
        ),
        pa.field(
            "name",
            pa.list_(
                pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
            ),
            nullable=False,
        ),
        pa.field(
            "protein",
            pa.struct(
                [
                    pa.field(
                        "recommendedName",
                        pa.struct(
                            [
                                pa.field(
                                    "fullName",
                                    pa.struct(
                                        [
                                            pa.field(
                                                "_evidence",
                                                pa.list_(
                                                    pa.field(
                                                        "item",
                                                        pa.int32(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:int"},
                                                    ),
                                                ),
                                                nullable=True,
                                            ),
                                            pa.field(
                                                "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                            ),
                                        ]
                                    ),
                                    nullable=False,
                                ),
                                pa.field(
                                    "shortName",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                    "ecNumber",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                    pa.field(
                        "alternativeName",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "fullName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                            "shortName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "ecNumber",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                        nullable=True,
                    ),
                    pa.field(
                        "submittedName",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "fullName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                            "ecNumber",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                        nullable=True,
                    ),
                    pa.field(
                        "allergenName",
                        pa.struct(
                            [
                                pa.field(
                                    "_evidence",
                                    pa.list_(
                                        pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                    ),
                                    nullable=True,
                                ),
                                pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "biotechName",
                        pa.struct(
                            [
                                pa.field(
                                    "_evidence",
                                    pa.list_(
                                        pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                    ),
                                    nullable=True,
                                ),
                                pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "cdAntigenName",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_evidence",
                                            pa.list_(
                                                pa.field(
                                                    "item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}
                                                ),
                                            ),
                                            nullable=True,
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
                        "innName",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_evidence",
                                            pa.list_(
                                                pa.field(
                                                    "item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}
                                                ),
                                            ),
                                            nullable=True,
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
                        "domain",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "recommendedName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "fullName",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_evidence",
                                                                    pa.list_(
                                                                        pa.field(
                                                                            "item",
                                                                            pa.int32(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:int"},
                                                                        ),
                                                                    ),
                                                                    nullable=True,
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
                                                        "shortName",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                                        "ecNumber",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                        pa.field(
                                            "alternativeName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "fullName",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                                                "shortName",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_evidence",
                                                                                    pa.list_(
                                                                                        pa.field(
                                                                                            "item",
                                                                                            pa.int32(),
                                                                                            nullable=True,
                                                                                            metadata={
                                                                                                "xsd_type": "xs:int"
                                                                                            },
                                                                                        ),
                                                                                    ),
                                                                                    nullable=True,
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
                                                                "ecNumber",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_evidence",
                                                                                    pa.list_(
                                                                                        pa.field(
                                                                                            "item",
                                                                                            pa.int32(),
                                                                                            nullable=True,
                                                                                            metadata={
                                                                                                "xsd_type": "xs:int"
                                                                                            },
                                                                                        ),
                                                                                    ),
                                                                                    nullable=True,
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
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "submittedName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "fullName",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                                                "ecNumber",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_evidence",
                                                                                    pa.list_(
                                                                                        pa.field(
                                                                                            "item",
                                                                                            pa.int32(),
                                                                                            nullable=True,
                                                                                            metadata={
                                                                                                "xsd_type": "xs:int"
                                                                                            },
                                                                                        ),
                                                                                    ),
                                                                                    nullable=True,
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
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "allergenName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                            "biotechName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                            "cdAntigenName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "innName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                        nullable=True,
                        metadata={
                            "comment": 'Describes names of "domains". Equivalent to the flat file DE-line Includes: section.'
                        },
                    ),
                    pa.field(
                        "component",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "recommendedName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "fullName",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_evidence",
                                                                    pa.list_(
                                                                        pa.field(
                                                                            "item",
                                                                            pa.int32(),
                                                                            nullable=True,
                                                                            metadata={"xsd_type": "xs:int"},
                                                                        ),
                                                                    ),
                                                                    nullable=True,
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
                                                        "shortName",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                                        "ecNumber",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                        pa.field(
                                            "alternativeName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "fullName",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                                                "shortName",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_evidence",
                                                                                    pa.list_(
                                                                                        pa.field(
                                                                                            "item",
                                                                                            pa.int32(),
                                                                                            nullable=True,
                                                                                            metadata={
                                                                                                "xsd_type": "xs:int"
                                                                                            },
                                                                                        ),
                                                                                    ),
                                                                                    nullable=True,
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
                                                                "ecNumber",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_evidence",
                                                                                    pa.list_(
                                                                                        pa.field(
                                                                                            "item",
                                                                                            pa.int32(),
                                                                                            nullable=True,
                                                                                            metadata={
                                                                                                "xsd_type": "xs:int"
                                                                                            },
                                                                                        ),
                                                                                    ),
                                                                                    nullable=True,
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
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "submittedName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "fullName",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_evidence",
                                                                            pa.list_(
                                                                                pa.field(
                                                                                    "item",
                                                                                    pa.int32(),
                                                                                    nullable=True,
                                                                                    metadata={"xsd_type": "xs:int"},
                                                                                ),
                                                                            ),
                                                                            nullable=True,
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
                                                                "ecNumber",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_evidence",
                                                                                    pa.list_(
                                                                                        pa.field(
                                                                                            "item",
                                                                                            pa.int32(),
                                                                                            nullable=True,
                                                                                            metadata={
                                                                                                "xsd_type": "xs:int"
                                                                                            },
                                                                                        ),
                                                                                    ),
                                                                                    nullable=True,
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
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "allergenName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                            "biotechName",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                            "cdAntigenName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "innName",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                        nullable=True,
                        metadata={
                            "comment": "Describes names of processed products. Equivalent to the flat file DE-line Contains: section."
                        },
                    ),
                ]
            ),
            nullable=False,
            metadata={
                "comment": "Describes the names for the protein and parts thereof. Equivalent to the flat file DE-line."
            },
        ),
        pa.field(
            "gene",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "name",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_evidence",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.int32(),
                                                            nullable=True,
                                                            metadata={"xsd_type": "xs:int"},
                                                        ),
                                                    ),
                                                    nullable=True,
                                                ),
                                                pa.field(
                                                    "_type",
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
                                nullable=False,
                                metadata={
                                    "comment": "Describes different types of gene designations. Equivalent to the flat file GN-line."
                                },
                            ),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
            metadata={"comment": "Describes a gene. Equivalent to the flat file GN-line."},
        ),
        pa.field(
            "organism",
            pa.struct(
                [
                    pa.field(
                        "_evidence",
                        pa.list_(
                            pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                        ),
                        nullable=True,
                    ),
                    pa.field(
                        "name",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_type", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ]
                                ),
                                nullable=True,
                            ),
                        ),
                        nullable=False,
                        metadata={
                            "comment": "Describes the names of the source organism. Equivalent to the flat file OS-line."
                        },
                    ),
                    pa.field(
                        "dbReference",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_type",
                                            pa.string(),
                                            nullable=False,
                                            metadata={
                                                "comment": "Describes the name of the database.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_id",
                                            pa.string(),
                                            nullable=False,
                                            metadata={
                                                "comment": "Describes a unique database identifier.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_evidence",
                                            pa.list_(
                                                pa.field(
                                                    "item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}
                                                ),
                                            ),
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "molecule",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_id",
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
                                            metadata={"comment": "Describes a molecule by name or unique identifier."},
                                        ),
                                        pa.field(
                                            "property",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_type",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "_value",
                                                                pa.string(),
                                                                nullable=False,
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
                        metadata={
                            "comment": "Describes a cross-reference to the NCBI taxonomy database. Equivalent to the flat file OX-line."
                        },
                    ),
                    pa.field(
                        "lineage",
                        pa.struct(
                            [
                                pa.field(
                                    "taxon",
                                    pa.list_(
                                        pa.field(
                                            "item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ),
                                    nullable=False,
                                ),
                            ]
                        ),
                        nullable=True,
                        metadata={
                            "comment": "Describes the lineage of the source organism. Equivalent to the flat file OC-line."
                        },
                    ),
                ]
            ),
            nullable=False,
            metadata={"comment": "Describes the source organism."},
        ),
        pa.field(
            "organismHost",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "_evidence",
                                pa.list_(
                                    pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                ),
                                nullable=True,
                            ),
                            pa.field(
                                "name",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_type",
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
                                nullable=False,
                                metadata={
                                    "comment": "Describes the names of the source organism. Equivalent to the flat file OS-line."
                                },
                            ),
                            pa.field(
                                "dbReference",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_type",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={
                                                        "comment": "Describes the name of the database.",
                                                        "xsd_type": "xs:string",
                                                    },
                                                ),
                                                pa.field(
                                                    "_id",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={
                                                        "comment": "Describes a unique database identifier.",
                                                        "xsd_type": "xs:string",
                                                    },
                                                ),
                                                pa.field(
                                                    "_evidence",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.int32(),
                                                            nullable=True,
                                                            metadata={"xsd_type": "xs:int"},
                                                        ),
                                                    ),
                                                    nullable=True,
                                                ),
                                                pa.field(
                                                    "molecule",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_id",
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
                                                    metadata={
                                                        "comment": "Describes a molecule by name or unique identifier."
                                                    },
                                                ),
                                                pa.field(
                                                    "property",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.struct(
                                                                [
                                                                    pa.field(
                                                                        "_type",
                                                                        pa.string(),
                                                                        nullable=False,
                                                                        metadata={"xsd_type": "xs:string"},
                                                                    ),
                                                                    pa.field(
                                                                        "_value",
                                                                        pa.string(),
                                                                        nullable=False,
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
                                metadata={
                                    "comment": "Describes a cross-reference to the NCBI taxonomy database. Equivalent to the flat file OX-line."
                                },
                            ),
                            pa.field(
                                "lineage",
                                pa.struct(
                                    [
                                        pa.field(
                                            "taxon",
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
                                metadata={
                                    "comment": "Describes the lineage of the source organism. Equivalent to the flat file OC-line."
                                },
                            ),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
            metadata={"comment": "Describes the source organism."},
        ),
        pa.field(
            "geneLocation",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_type", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                            pa.field(
                                "_evidence",
                                pa.list_(
                                    pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                ),
                                nullable=True,
                            ),
                            pa.field(
                                "name",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_status",
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
                                metadata={"comment": "Indicates whether the name of a plasmid is known or unknown."},
                            ),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
            metadata={
                "comment": "Describes non-nuclear gene locations (organelles and plasmids). Equivalent to the flat file OG-line."
            },
        ),
        pa.field(
            "reference",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "_evidence",
                                pa.list_(
                                    pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                ),
                                nullable=True,
                            ),
                            pa.field("_key", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                            pa.field(
                                "citation",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_type",
                                            pa.string(),
                                            nullable=False,
                                            metadata={
                                                "comment": "Describes the type of a citation.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_date",
                                            pa.string(),
                                            nullable=True,
                                            metadata={"xsd_type": "xs:date | xs:gYearMonth | xs:gYear"},
                                        ),
                                        pa.field(
                                            "_name",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the name of an (online) journal or book.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_volume",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the volume of a journal or book.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_first",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the first page of an article.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_last",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the last page of an article.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_publisher",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the publisher of a book.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_city",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the city where a book was published.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_db",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the database name of submissions.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_number",
                                            pa.string(),
                                            nullable=True,
                                            metadata={"comment": "Describes a patent number.", "xsd_type": "xs:string"},
                                        ),
                                        pa.field(
                                            "_institute",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the institute where a thesis was made.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_country",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the country where a thesis was made.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "title",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the title of a citation. Equivalent to the flat file RT-line.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "editorList",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "consortium",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_name",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ),
                                                        nullable=True,
                                                        metadata={
                                                            "comment": "Describes the authors of a citation when these are represented by a consortium. Equivalent to the flat file RG-line."
                                                        },
                                                    ),
                                                    pa.field(
                                                        "person",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_name",
                                                                            pa.string(),
                                                                            nullable=False,
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
                                            metadata={
                                                "comment": "Describes the editors of a book (only used for books). Equivalent to part of the flat file RL-line of books."
                                            },
                                        ),
                                        pa.field(
                                            "authorList",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "consortium",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_name",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                    ]
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ),
                                                        nullable=True,
                                                        metadata={
                                                            "comment": "Describes the authors of a citation when these are represented by a consortium. Equivalent to the flat file RG-line."
                                                        },
                                                    ),
                                                    pa.field(
                                                        "person",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_name",
                                                                            pa.string(),
                                                                            nullable=False,
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
                                            metadata={
                                                "comment": "Describes the authors of a citation. Equivalent to the flat file RA-line."
                                            },
                                        ),
                                        pa.field(
                                            "locator",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes the location (URL) of an online journal article. No flat file equivalent.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "dbReference",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_type",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes the name of the database.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_id",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes a unique database identifier.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "molecule",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_id",
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
                                                                metadata={
                                                                    "comment": "Describes a molecule by name or unique identifier."
                                                                },
                                                            ),
                                                            pa.field(
                                                                "property",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_type",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "_value",
                                                                                    pa.string(),
                                                                                    nullable=False,
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
                                            nullable=True,
                                            metadata={
                                                "comment": "Describes cross-references to bibliography databases (MEDLINE, PubMed, AGRICOLA) or other online resources (DOI). Equivalent to the flat file RX-line."
                                            },
                                        ),
                                    ]
                                ),
                                nullable=False,
                                metadata={
                                    "comment": "Describes different types of citations. Equivalent to the flat file RX-, RG-, RA-, RT- and RL-lines."
                                },
                            ),
                            pa.field(
                                "scope",
                                pa.list_(
                                    pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                ),
                                nullable=False,
                                metadata={
                                    "comment": "Describes the scope of a citation. Equivalent to the flat file RP-line."
                                },
                            ),
                            pa.field(
                                "source",
                                pa.struct(
                                    [
                                        pa.field(
                                            "strain",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "plasmid",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "transposon",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "tissue",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                metadata={
                                    "comment": "Describes the source of the sequence according to the citation. Equivalent to the flat file RC-line."
                                },
                            ),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=False,
            metadata={
                "comment": "Describes a citation and a summary of its content. Equivalent to the flat file RN-, RP-, RC-, RX-, RG-, RA-, RT- and RL-lines."
            },
        ),
        pa.field(
            "comment",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "_type",
                                pa.string(),
                                nullable=True,
                                metadata={
                                    "comment": 'Describes the type of a general annotation. Equivalent to the flat file CC comment topics (except for "DATABASE" which is translated to "online information").',
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "_locationType",
                                pa.string(),
                                nullable=True,
                                metadata={
                                    "comment": 'Describes the type of sequence location in \'RNA editing\' annotations. Common values are "Not_applicable" and "Undetermined".',
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "_name",
                                pa.string(),
                                nullable=True,
                                metadata={
                                    "comment": "Describes an optional name for an 'online information'.",
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "_mass",
                                pa.float64(),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the molecular mass in 'mass spectrometry' annotations.",
                                    "xsd_type": "xs:float",
                                },
                            ),
                            pa.field(
                                "_error",
                                pa.string(),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the error of the mass measurement in 'mass spectrometry' annotations.",
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "_method",
                                pa.string(),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the experimental method in 'mass spectrometry' annotations.",
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "_evidence",
                                pa.list_(
                                    pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                ),
                                nullable=True,
                            ),
                            pa.field(
                                "molecule",
                                pa.struct(
                                    [
                                        pa.field("_id", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                        pa.field(
                                            "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={"comment": "Describes a molecule by name or unique identifier."},
                            ),
                            pa.field(
                                "absorption",
                                pa.struct(
                                    [
                                        pa.field(
                                            "max",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
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
                                            "text",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                            pa.field(
                                "kinetics",
                                pa.struct(
                                    [
                                        pa.field(
                                            "KM",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "Vmax",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                            "text",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                            pa.field(
                                "phDependence",
                                pa.struct(
                                    [
                                        pa.field(
                                            "text",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                            pa.field(
                                "redoxPotential",
                                pa.struct(
                                    [
                                        pa.field(
                                            "text",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                            pa.field(
                                "temperatureDependence",
                                pa.struct(
                                    [
                                        pa.field(
                                            "text",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                            pa.field(
                                "reaction",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_evidence",
                                            pa.list_(
                                                pa.field(
                                                    "item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}
                                                ),
                                            ),
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "text", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "dbReference",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_type",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes the name of the database.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_id",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes a unique database identifier.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "molecule",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_id",
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
                                                                metadata={
                                                                    "comment": "Describes a molecule by name or unique identifier."
                                                                },
                                                            ),
                                                            pa.field(
                                                                "property",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_type",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "_value",
                                                                                    pa.string(),
                                                                                    nullable=False,
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
                                            metadata={
                                                "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                            },
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={"comment": "Describes a chemical reaction."},
                            ),
                            pa.field(
                                "physiologicalReaction",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_direction",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                                pa.field(
                                                    "_evidence",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.int32(),
                                                            nullable=True,
                                                            metadata={"xsd_type": "xs:int"},
                                                        ),
                                                    ),
                                                    nullable=True,
                                                ),
                                                pa.field(
                                                    "dbReference",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_type",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes the name of the database.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_id",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes a unique database identifier.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "molecule",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_id",
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
                                                                metadata={
                                                                    "comment": "Describes a molecule by name or unique identifier."
                                                                },
                                                            ),
                                                            pa.field(
                                                                "property",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_type",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "_value",
                                                                                    pa.string(),
                                                                                    nullable=False,
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
                                                    nullable=False,
                                                    metadata={
                                                        "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                                    },
                                                ),
                                            ]
                                        ),
                                        nullable=True,
                                    ),
                                ),
                                nullable=True,
                                metadata={"comment": "Describes a physiological reaction."},
                            ),
                            pa.field(
                                "cofactor",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_evidence",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.int32(),
                                                            nullable=True,
                                                            metadata={"xsd_type": "xs:int"},
                                                        ),
                                                    ),
                                                    nullable=True,
                                                ),
                                                pa.field(
                                                    "name",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                                pa.field(
                                                    "dbReference",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_type",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes the name of the database.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_id",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes a unique database identifier.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "molecule",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_id",
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
                                                                metadata={
                                                                    "comment": "Describes a molecule by name or unique identifier."
                                                                },
                                                            ),
                                                            pa.field(
                                                                "property",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_type",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "_value",
                                                                                    pa.string(),
                                                                                    nullable=False,
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
                                                    nullable=False,
                                                    metadata={
                                                        "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                                    },
                                                ),
                                            ]
                                        ),
                                        nullable=True,
                                    ),
                                ),
                                nullable=True,
                                metadata={"comment": "Describes a cofactor."},
                            ),
                            pa.field(
                                "subcellularLocation",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "location",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.struct(
                                                                [
                                                                    pa.field(
                                                                        "_evidence",
                                                                        pa.list_(
                                                                            pa.field(
                                                                                "item",
                                                                                pa.int32(),
                                                                                nullable=True,
                                                                                metadata={"xsd_type": "xs:int"},
                                                                            ),
                                                                        ),
                                                                        nullable=True,
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
                                                    "topology",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.struct(
                                                                [
                                                                    pa.field(
                                                                        "_evidence",
                                                                        pa.list_(
                                                                            pa.field(
                                                                                "item",
                                                                                pa.int32(),
                                                                                nullable=True,
                                                                                metadata={"xsd_type": "xs:int"},
                                                                            ),
                                                                        ),
                                                                        nullable=True,
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
                                                    "orientation",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.struct(
                                                                [
                                                                    pa.field(
                                                                        "_evidence",
                                                                        pa.list_(
                                                                            pa.field(
                                                                                "item",
                                                                                pa.int32(),
                                                                                nullable=True,
                                                                                metadata={"xsd_type": "xs:int"},
                                                                            ),
                                                                        ),
                                                                        nullable=True,
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
                                nullable=True,
                                metadata={
                                    "comment": "Describes the subcellular location and optionally the topology and orientation of a molecule."
                                },
                            ),
                            pa.field(
                                "conflict",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_type", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "_ref",
                                            pa.string(),
                                            nullable=True,
                                            metadata={
                                                "comment": "Refers to the 'key' attribute of a 'reference' element.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "sequence",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_resource",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_id",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_version",
                                                        pa.int32(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:int"},
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={"comment": "Used in 'sequence caution' annotations."},
                            ),
                            pa.field(
                                "link",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_uri",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={"xsd_type": "xs:anyURI"},
                                                ),
                                            ]
                                        ),
                                        nullable=True,
                                    ),
                                ),
                                nullable=True,
                                metadata={"comment": "Used in 'online information' annotations."},
                            ),
                            pa.field(
                                "event",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_type",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                            ]
                                        ),
                                        nullable=True,
                                    ),
                                ),
                                nullable=True,
                                metadata={"comment": "Describes the type of events that cause alternative products."},
                            ),
                            pa.field(
                                "isoform",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "id",
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
                                                pa.field(
                                                    "name",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.struct(
                                                                [
                                                                    pa.field(
                                                                        "_evidence",
                                                                        pa.list_(
                                                                            pa.field(
                                                                                "item",
                                                                                pa.int32(),
                                                                                nullable=True,
                                                                                metadata={"xsd_type": "xs:int"},
                                                                            ),
                                                                        ),
                                                                        nullable=True,
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
                                                    "sequence",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_type",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "_ref",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                        ]
                                                    ),
                                                    nullable=False,
                                                ),
                                                pa.field(
                                                    "text",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.struct(
                                                                [
                                                                    pa.field(
                                                                        "_evidence",
                                                                        pa.list_(
                                                                            pa.field(
                                                                                "item",
                                                                                pa.int32(),
                                                                                nullable=True,
                                                                                metadata={"xsd_type": "xs:int"},
                                                                            ),
                                                                        ),
                                                                        nullable=True,
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
                                nullable=True,
                                metadata={"comment": "Describes isoforms in 'alternative products' annotations."},
                            ),
                            pa.field(
                                "interactant",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_intactId",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                                pa.field(
                                                    "id", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                                ),
                                                pa.field(
                                                    "label",
                                                    pa.string(),
                                                    nullable=True,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                                pa.field(
                                                    "dbReference",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_type",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes the name of the database.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_id",
                                                                pa.string(),
                                                                nullable=False,
                                                                metadata={
                                                                    "comment": "Describes a unique database identifier.",
                                                                    "xsd_type": "xs:string",
                                                                },
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                            pa.field(
                                                                "molecule",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_id",
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
                                                                metadata={
                                                                    "comment": "Describes a molecule by name or unique identifier."
                                                                },
                                                            ),
                                                            pa.field(
                                                                "property",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.struct(
                                                                            [
                                                                                pa.field(
                                                                                    "_type",
                                                                                    pa.string(),
                                                                                    nullable=False,
                                                                                    metadata={"xsd_type": "xs:string"},
                                                                                ),
                                                                                pa.field(
                                                                                    "_value",
                                                                                    pa.string(),
                                                                                    nullable=False,
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
                                                    metadata={
                                                        "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                                    },
                                                ),
                                            ]
                                        ),
                                        nullable=True,
                                    ),
                                ),
                                nullable=True,
                            ),
                            pa.field("organismsDiffer", pa.bool_(), nullable=True, metadata={"xsd_type": "xs:boolean"}),
                            pa.field("experiments", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                            pa.field(
                                "disease",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "name", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "acronym", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "description",
                                            pa.string(),
                                            nullable=False,
                                            metadata={"xsd_type": "xs:string"},
                                        ),
                                        pa.field(
                                            "dbReference",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_type",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes the name of the database.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_id",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes a unique database identifier.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                    pa.field(
                                                        "molecule",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_id",
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
                                                        metadata={
                                                            "comment": "Describes a molecule by name or unique identifier."
                                                        },
                                                    ),
                                                    pa.field(
                                                        "property",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_type",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "_value",
                                                                            pa.string(),
                                                                            nullable=False,
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
                                            nullable=False,
                                            metadata={
                                                "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                            },
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={"comment": "Used in 'disease' annotations."},
                            ),
                            pa.field(
                                "location",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_sequence",
                                                    pa.string(),
                                                    nullable=True,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                                pa.field(
                                                    "begin",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_position",
                                                                pa.uint64(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:unsignedLong"},
                                                            ),
                                                            pa.field(
                                                                "_status",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ]
                                                    ),
                                                    nullable=True,
                                                ),
                                                pa.field(
                                                    "end",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_position",
                                                                pa.uint64(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:unsignedLong"},
                                                            ),
                                                            pa.field(
                                                                "_status",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
                                                            ),
                                                        ]
                                                    ),
                                                    nullable=True,
                                                ),
                                                pa.field(
                                                    "position",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_position",
                                                                pa.uint64(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:unsignedLong"},
                                                            ),
                                                            pa.field(
                                                                "_status",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "_evidence",
                                                                pa.list_(
                                                                    pa.field(
                                                                        "item",
                                                                        pa.int32(),
                                                                        nullable=True,
                                                                        metadata={"xsd_type": "xs:int"},
                                                                    ),
                                                                ),
                                                                nullable=True,
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
                                metadata={"comment": "Used in 'mass spectrometry' and 'sequence caution' annotations."},
                            ),
                            pa.field(
                                "text",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_evidence",
                                                    pa.list_(
                                                        pa.field(
                                                            "item",
                                                            pa.int32(),
                                                            nullable=True,
                                                            metadata={"xsd_type": "xs:int"},
                                                        ),
                                                    ),
                                                    nullable=True,
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
            nullable=True,
            metadata={
                "comment": "Describes different types of general annotations. Equivalent to the flat file CC-line."
            },
        ),
        pa.field(
            "dbReference",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "_type",
                                pa.string(),
                                nullable=False,
                                metadata={"comment": "Describes the name of the database.", "xsd_type": "xs:string"},
                            ),
                            pa.field(
                                "_id",
                                pa.string(),
                                nullable=False,
                                metadata={
                                    "comment": "Describes a unique database identifier.",
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "_evidence",
                                pa.list_(
                                    pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                ),
                                nullable=True,
                            ),
                            pa.field(
                                "molecule",
                                pa.struct(
                                    [
                                        pa.field("_id", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                        pa.field(
                                            "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={"comment": "Describes a molecule by name or unique identifier."},
                            ),
                            pa.field(
                                "property",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_type",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                                pa.field(
                                                    "_value",
                                                    pa.string(),
                                                    nullable=False,
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
            nullable=True,
            metadata={"comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."},
        ),
        pa.field(
            "proteinExistence",
            pa.struct(
                [
                    pa.field("_type", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                ]
            ),
            nullable=False,
            metadata={
                "comment": "Describes the evidence for the protein's existence. Equivalent to the flat file PE-line."
            },
        ),
        pa.field(
            "keyword",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "_evidence",
                                pa.list_(
                                    pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                ),
                                nullable=True,
                            ),
                            pa.field("_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                            pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
        pa.field(
            "feature",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "_type",
                                pa.string(),
                                nullable=False,
                                metadata={
                                    "comment": "Describes the type of a sequence annotation. Equivalent to the flat file FT feature keys, but using full terms instead of acronyms.",
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field("_id", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("_description", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field(
                                "_evidence",
                                pa.list_(
                                    pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                                ),
                                nullable=True,
                            ),
                            pa.field("_ref", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field(
                                "original",
                                pa.string(),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the original sequence in annotations that describe natural or artifical sequence variations.",
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "variation",
                                pa.list_(
                                    pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                                ),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the variant sequence in annotations that describe natural or artifical sequence variations."
                                },
                            ),
                            pa.field(
                                "location",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_sequence", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "begin",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_position",
                                                        pa.uint64(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:unsignedLong"},
                                                    ),
                                                    pa.field(
                                                        "_status",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "end",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_position",
                                                        pa.uint64(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:unsignedLong"},
                                                    ),
                                                    pa.field(
                                                        "_status",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                        pa.field(
                                            "position",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_position",
                                                        pa.uint64(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:unsignedLong"},
                                                    ),
                                                    pa.field(
                                                        "_status",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                ]
                                            ),
                                            nullable=True,
                                        ),
                                    ]
                                ),
                                nullable=False,
                                metadata={"comment": "Describes the sequence coordinates of the annotation."},
                            ),
                            pa.field(
                                "ligand",
                                pa.struct(
                                    [
                                        pa.field(
                                            "name", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "dbReference",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_type",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes the name of the database.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_id",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes a unique database identifier.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                    pa.field(
                                                        "molecule",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_id",
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
                                                        metadata={
                                                            "comment": "Describes a molecule by name or unique identifier."
                                                        },
                                                    ),
                                                    pa.field(
                                                        "property",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_type",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "_value",
                                                                            pa.string(),
                                                                            nullable=False,
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
                                            metadata={
                                                "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                            },
                                        ),
                                        pa.field(
                                            "label", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "note", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the chemical entity that is bound in annotations that describe binding sites."
                                },
                            ),
                            pa.field(
                                "ligandPart",
                                pa.struct(
                                    [
                                        pa.field(
                                            "name", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "dbReference",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_type",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes the name of the database.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_id",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes a unique database identifier.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                    pa.field(
                                                        "molecule",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_id",
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
                                                        metadata={
                                                            "comment": "Describes a molecule by name or unique identifier."
                                                        },
                                                    ),
                                                    pa.field(
                                                        "property",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_type",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "_value",
                                                                            pa.string(),
                                                                            nullable=False,
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
                                            metadata={
                                                "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                            },
                                        ),
                                        pa.field(
                                            "label", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "note", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the specific part of a molecule that is bound in annotations that describe binding sites."
                                },
                            ),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
            metadata={
                "comment": "Describes different types of sequence annotations. Equivalent to the flat file FT-line."
            },
        ),
        pa.field(
            "evidence",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field(
                                "_type",
                                pa.string(),
                                nullable=False,
                                metadata={
                                    "comment": "Describes the type of an evidence using the Evidence Code Ontology (http://www.obofoundry.org/cgi-bin/detail.cgi?id=evidence_code).",
                                    "xsd_type": "xs:string",
                                },
                            ),
                            pa.field(
                                "_key",
                                pa.int64(),
                                nullable=False,
                                metadata={
                                    "comment": "A unique key to link annotations (via 'evidence' attributes) to evidences.",
                                    "xsd_type": "xs:integer",
                                },
                            ),
                            pa.field(
                                "source",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_ref", pa.int64(), nullable=True, metadata={"xsd_type": "xs:integer"}
                                        ),
                                        pa.field(
                                            "dbReference",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_type",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes the name of the database.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_id",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes a unique database identifier.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                    pa.field(
                                                        "molecule",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_id",
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
                                                        metadata={
                                                            "comment": "Describes a molecule by name or unique identifier."
                                                        },
                                                    ),
                                                    pa.field(
                                                        "property",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_type",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "_value",
                                                                            pa.string(),
                                                                            nullable=False,
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
                                            metadata={
                                                "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                            },
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the source of the data using a database cross-reference (or a 'ref' attribute when the source cannot be found in a public data source, such as PubMed, and is cited only within the UniProtKB entry)."
                                },
                            ),
                            pa.field(
                                "importedFrom",
                                pa.struct(
                                    [
                                        pa.field(
                                            "dbReference",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_type",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes the name of the database.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_id",
                                                        pa.string(),
                                                        nullable=False,
                                                        metadata={
                                                            "comment": "Describes a unique database identifier.",
                                                            "xsd_type": "xs:string",
                                                        },
                                                    ),
                                                    pa.field(
                                                        "_evidence",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.int32(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:int"},
                                                            ),
                                                        ),
                                                        nullable=True,
                                                    ),
                                                    pa.field(
                                                        "molecule",
                                                        pa.struct(
                                                            [
                                                                pa.field(
                                                                    "_id",
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
                                                        metadata={
                                                            "comment": "Describes a molecule by name or unique identifier."
                                                        },
                                                    ),
                                                    pa.field(
                                                        "property",
                                                        pa.list_(
                                                            pa.field(
                                                                "item",
                                                                pa.struct(
                                                                    [
                                                                        pa.field(
                                                                            "_type",
                                                                            pa.string(),
                                                                            nullable=False,
                                                                            metadata={"xsd_type": "xs:string"},
                                                                        ),
                                                                        pa.field(
                                                                            "_value",
                                                                            pa.string(),
                                                                            nullable=False,
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
                                            nullable=False,
                                            metadata={
                                                "comment": "Describes a database cross-reference. Equivalent to the flat file DR-line."
                                            },
                                        ),
                                    ]
                                ),
                                nullable=True,
                                metadata={
                                    "comment": "Describes the source of the evidence, when it is not assigned by UniProt, but imported from an external database."
                                },
                            ),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
            metadata={"comment": "Describes the evidence for an annotation. No flat file equivalent."},
        ),
        pa.field(
            "sequence",
            pa.struct(
                [
                    pa.field("_length", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}),
                    pa.field("_mass", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}),
                    pa.field("_checksum", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                    pa.field("_modified", pa.date32(), nullable=False, metadata={"xsd_type": "xs:date"}),
                    pa.field("_version", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}),
                    pa.field("_precursor", pa.bool_(), nullable=True, metadata={"xsd_type": "xs:boolean"}),
                    pa.field("_fragment", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                ]
            ),
            nullable=False,
        ),
    ]
)
