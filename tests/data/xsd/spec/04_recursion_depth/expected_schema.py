"""Expected PyArrow schema for 04_recursion_depth (xml_tag='record')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "tree",
            pa.struct(
                [
                    pa.field("_label", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("_weight", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
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
                                        pa.field("_weight", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
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
            "Section",
            pa.struct(
                [
                    pa.field("SectionTitle", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "Section",
                        pa.list_(
                            pa.field(
                                "item",
                                pa.struct(
                                    [
                                        pa.field(
                                            "SectionTitle",
                                            pa.string(),
                                            nullable=False,
                                            metadata={"xsd_type": "xs:string"},
                                        ),
                                        pa.field(
                                            "Section",
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
    ]
)

EXPECTED_SCHEMA_DEPTH_1 = pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "tree",
            pa.struct(
                [
                    pa.field("_label", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("_weight", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
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
        pa.field(
            "Section",
            pa.struct(
                [
                    pa.field("SectionTitle", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                    pa.field(
                        "Section",
                        pa.list_(
                            pa.field("item", pa.string(), nullable=True),
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=True,
        ),
    ]
)
