"""Expected PyArrow schema for 09_empty_marker_elements (xml_tag='book')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("title", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field("BookInProcess", pa.string(), nullable=True),
        pa.field(
            "flag",
            pa.list_(
                pa.field("item", pa.string(), nullable=True),
            ),
            nullable=True,
        ),
    ]
)
