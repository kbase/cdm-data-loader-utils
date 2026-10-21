"""Expected PyArrow schema for 10_combinations (xml_tag='record')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field("plain", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "plainList",
            pa.list_(
                pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
            ),
            nullable=True,
        ),
        pa.field(
            "bounded",
            pa.list_(
                pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
            ),
            nullable=True,
        ),
        pa.field(
            "label",
            pa.struct(
                [
                    pa.field("_lang", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("_source", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "labelList",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_lang", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("_source", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
        pa.field("textNoAttr", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "measure",
            pa.struct(
                [
                    pa.field("_unit", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("#text", pa.float64(), nullable=True, metadata={"xsd_type": "xs:decimal"}),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "flag",
            pa.struct(
                [
                    pa.field("_on", pa.bool_(), nullable=True, metadata={"xsd_type": "xs:boolean"}),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "flagList",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_on", pa.bool_(), nullable=True, metadata={"xsd_type": "xs:boolean"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
        pa.field(
            "marker",
            pa.list_(
                pa.field("item", pa.string(), nullable=True),
            ),
            nullable=True,
        ),
        pa.field(
            "group",
            pa.struct(
                [
                    pa.field("name", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "alias",
                        pa.list_(
                            pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "box",
            pa.struct(
                [
                    pa.field(
                        "item",
                        pa.list_(
                            pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "entity",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_type", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("_name", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("name", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field(
                                "prop",
                                pa.list_(
                                    pa.field(
                                        "item",
                                        pa.struct(
                                            [
                                                pa.field(
                                                    "_key",
                                                    pa.string(),
                                                    nullable=False,
                                                    metadata={"xsd_type": "xs:string"},
                                                ),
                                                pa.field(
                                                    "_value",
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
            "para",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_style", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
        pa.field(
            "note",
            pa.struct(
                [
                    pa.field("_kind", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                ]
            ),
            nullable=True,
        ),
        pa.field("remark", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "optA",
            pa.list_(
                pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
            ),
            nullable=True,
        ),
        pa.field(
            "optB",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_lang", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("_source", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
        pa.field(
            "pairs",
            pa.struct(
                [
                    pa.field(
                        "k",
                        pa.list_(
                            pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ),
                        nullable=False,
                    ),
                    pa.field(
                        "v",
                        pa.list_(
                            pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ),
                        nullable=False,
                    ),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "book",
            pa.struct(
                [
                    pa.field("title", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "chapter",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "title",
                                            pa.list_(
                                                pa.field(
                                                    "item",
                                                    pa.struct(
                                                        [
                                                            pa.field(
                                                                "_lang",
                                                                pa.string(),
                                                                nullable=True,
                                                                metadata={"xsd_type": "xs:string"},
                                                            ),
                                                            pa.field(
                                                                "_source",
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
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "left",
            pa.struct(
                [
                    pa.field(
                        "node",
                        pa.struct(
                            [
                                pa.field("leaf", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                            ]
                        ),
                        nullable=False,
                    ),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "right",
            pa.struct(
                [
                    pa.field(
                        "node",
                        pa.struct(
                            [
                                pa.field(
                                    "leaf",
                                    pa.list_(
                                        pa.field(
                                            "item",
                                            pa.struct(
                                                [
                                                    pa.field(
                                                        "_lang",
                                                        pa.string(),
                                                        nullable=True,
                                                        metadata={"xsd_type": "xs:string"},
                                                    ),
                                                    pa.field(
                                                        "_source",
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
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "tree",
            pa.struct(
                [
                    pa.field("_label", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "branch",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "_label", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}
                                        ),
                                        pa.field(
                                            "branch",
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
        pa.field(
            "shared",
            pa.list_(
                pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
            ),
            nullable=True,
        ),
        pa.field(
            "derived",
            pa.struct(
                [
                    pa.field("_baseAttr", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("_extraAttr", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("baseField", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "extraField",
                        pa.list_(
                            pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=True,
        ),
    ]
)
