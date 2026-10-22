"""Expected PyArrow schema for 08_nullability_and_empty_tags (xml_tag='entry')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field("_source", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field("name", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field("summary", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field("optionalNote", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field("geneId", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
        pa.field("proteinId", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "comment",
            pa.struct(
                [
                    pa.field("_type", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                    pa.field("text", pa.string(), nullable=True, metadata={"xsd_type": "xs:string"}),
                ]
            ),
            nullable=True,
        ),
    ]
)
