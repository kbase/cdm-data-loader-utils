"""Expected PyArrow schema for 03_mixed_content_raw_text (xml_tag='article')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("title", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "abstractText",
            pa.list_(
                pa.field(
                    "item",
                    pa.struct(
                        [
                            pa.field("_Label", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                            pa.field("#text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                        ]
                    ),
                    nullable=True,
                ),
            ),
            nullable=True,
        ),
    ]
)
