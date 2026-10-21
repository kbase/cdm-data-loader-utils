"""Expected PyArrow schema for 01_scalars_and_simple_types (xml_tag='record')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
        pa.field(
            "_id", pa.string(), nullable=False, metadata={"comment": "Unique record identifier.", "xsd_type": "xs:ID"}
        ),
        pa.field(
            "_status",
            pa.string(),
            nullable=False,
            metadata={"comment": "Record lifecycle status.", "xsd_type": "xs:string"},
        ),
        pa.field(
            "_evidence",
            pa.list_(
                pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
            ),
            nullable=True,
        ),
        pa.field("_version", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}),
        pa.field(
            "flag", pa.bool_(), nullable=False, metadata={"comment": "Primary boolean flag.", "xsd_type": "xs:boolean"}
        ),
        pa.field("smallInt", pa.int8(), nullable=True, metadata={"xsd_type": "xs:byte"}),
        pa.field("medInt", pa.int16(), nullable=True, metadata={"xsd_type": "xs:short"}),
        pa.field("stdInt", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
        pa.field("bigInt", pa.int64(), nullable=True, metadata={"xsd_type": "xs:integer"}),
        pa.field("uByte", pa.uint8(), nullable=True, metadata={"xsd_type": "xs:unsignedByte"}),
        pa.field("uShort", pa.uint16(), nullable=True, metadata={"xsd_type": "xs:unsignedShort"}),
        pa.field("uInt", pa.uint32(), nullable=True, metadata={"xsd_type": "xs:unsignedInt"}),
        pa.field("uLong", pa.uint64(), nullable=True, metadata={"xsd_type": "xs:unsignedLong"}),
        pa.field("f32", pa.float64(), nullable=True, metadata={"xsd_type": "xs:float"}),
        pa.field("f64", pa.float64(), nullable=True, metadata={"xsd_type": "xs:double"}),
        pa.field("dec", pa.float64(), nullable=True, metadata={"xsd_type": "xs:decimal"}),
        pa.field("created", pa.date32(), nullable=True, metadata={"xsd_type": "xs:date"}),
        pa.field("updatedAt", pa.string(), nullable=True, metadata={"xsd_type": "xs:dateTime"}),
        pa.field("published", pa.string(), nullable=True, metadata={"xsd_type": "xs:date | xs:gYearMonth | xs:gYear"}),
        pa.field("homogeneousUnion", pa.int64(), nullable=True, metadata={"xsd_type": "xs:long | xs:integer"}),
        pa.field(
            "scores",
            pa.list_(
                pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
            ),
            nullable=True,
        ),
        pa.field(
            "tokens",
            pa.list_(
                pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": "xs:NMTOKEN"}),
            ),
            nullable=True,
        ),
        pa.field("rawHex", pa.binary(), nullable=True, metadata={"xsd_type": "xs:hexBinary"}),
        pa.field("rawB64", pa.binary(), nullable=True, metadata={"xsd_type": "xs:base64Binary"}),
    ]
)
