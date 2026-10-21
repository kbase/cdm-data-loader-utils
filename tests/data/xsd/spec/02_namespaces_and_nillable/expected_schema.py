"""Expected PyArrow schema for 02_namespaces_and_nillable (xml_tag='entry')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field("title", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "label",
            pa.struct(
                [
                    pa.field("_lang", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "comment",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_type", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
    ]
)
