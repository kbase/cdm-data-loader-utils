"""Expected PyArrow schema for 05_groups_and_compositors (xml_tag='item')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("_createdBy", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field("_version", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
        pa.field(
            "pagination",
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
            "interaction",
            pa.struct(
                [
                    pa.field(
                        "interactant",
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
            "meta",
            pa.struct(
                [
                    pa.field("author", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
                    pa.field("year", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                ]
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
                            pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
                        ),
                        nullable=False,
                    ),
                ]
            ),
            nullable=True,
        ),
    ]
)
