"""Expected PyArrow schema for uniref (xml_tag='entry')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field(
            "_id", pa.string(), nullable=False, metadata={"comment": "The entry identifier.", "xsd_type": "xs:ID"}
        ),
        pa.field(
            "_updated",
            pa.date32(),
            nullable=False,
            metadata={"comment": "The update date for the entry.", "xsd_type": "xs:date"},
        ),
        pa.field("name", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "property",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_type", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                            pa.field("_value", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
        pa.field(
            "representativeMember",
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
                                        "comment": "The name of the database this cross-reference is referring to.",
                                        "xsd_type": "xs:string",
                                    },
                                ),
                                pa.field(
                                    "_id",
                                    pa.string(),
                                    nullable=False,
                                    metadata={"comment": "The ID referred to.", "xsd_type": "xs:string"},
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
                        metadata={"comment": "Database cross-references"},
                    ),
                    pa.field(
                        "sequence",
                        pa.struct(
                            [
                                pa.field("_length", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}),
                                pa.field("_checksum", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                                pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=False,
            metadata={"comment": "UniRef50/90/100 cluster members"},
        ),
        pa.field(
            "member",
            pa.list_(
                pa.field(
                    "item",
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
                                                "comment": "The name of the database this cross-reference is referring to.",
                                                "xsd_type": "xs:string",
                                            },
                                        ),
                                        pa.field(
                                            "_id",
                                            pa.string(),
                                            nullable=False,
                                            metadata={"comment": "The ID referred to.", "xsd_type": "xs:string"},
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
                                metadata={"comment": "Database cross-references"},
                            ),
                            pa.field(
                                "sequence",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_length", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}
                                        ),
                                        pa.field(
                                            "_checksum", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
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
            metadata={"comment": "UniRef50/90/100 cluster members"},
        ),
    ]
)
