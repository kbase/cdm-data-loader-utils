"""Expected PyArrow schema for 06_context_dependent_paths (xml_tag='root')."""

# ruff: noqa: INP001

import pyarrow as pa

EXPECTED_SCHEMA = pa.schema(
    [
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
            nullable=False,
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
            nullable=False,
        ),
    ]
)
