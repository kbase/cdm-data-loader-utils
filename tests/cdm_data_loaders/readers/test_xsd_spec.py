"""Specification tests for XSD-to-PyArrow schema generation and deterministic XML-to-Arrow conversion."""

import copy
import datetime as dt
import importlib.util
import io
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import xmlschema
from _pytest.mark.structures import ParameterSet
from lxml import etree

from cdm_data_loaders.readers.xsd import (
    convert_xml_element,
    generate_pyarrow_schema,
    load_schema,
)


def _list_field(name: str, item_type: pa.DataType, *, nullable: bool = True) -> pa.Field:
    """Build a PyArrow list field with standard nullable item field."""
    return pa.field(name, pa.list_(pa.field("item", item_type, nullable=True)), nullable=nullable)


def _load_expected_schema(schema_file: Path, attr_name: str = "EXPECTED_SCHEMA") -> pa.Schema:
    """Load a pa.Schema constant from an expected_schema.py fixture file."""
    spec = importlib.util.spec_from_file_location(schema_file.stem, schema_file)
    if spec is None or spec.loader is None:
        msg = f"Could not load schema module from {schema_file}"
        raise RuntimeError(msg)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    schema = getattr(module, attr_name)
    if not isinstance(schema, pa.Schema):
        msg = f"Expected {attr_name} in {schema_file} to be a pa.Schema, got {type(schema)}"
        raise TypeError(msg)
    return schema


def _coerce_json_rows_for_arrow(value: Any, dtype: pa.DataType) -> Any:  # noqa: ANN401
    """Coerce ISO date strings and binary strings from expected JSON into Python objects for pa.Table.from_pylist."""
    if value is None:
        return None
    if pa.types.is_date32(dtype) and isinstance(value, str):
        return dt.date.fromisoformat(value)
    if pa.types.is_binary(dtype) and isinstance(value, str):
        return value.encode("utf-8")
    if pa.types.is_list(dtype) and isinstance(value, list):
        return [_coerce_json_rows_for_arrow(item, dtype.value_type) for item in value]
    if pa.types.is_struct(dtype) and isinstance(value, dict):
        return {
            field.name: _coerce_json_rows_for_arrow(value[field.name], field.type)
            for field in dtype
            if field.name in value
        }
    return value


def _coerce_rows_to_schema(rows: list[dict[str, Any]], schema: pa.Schema) -> list[dict[str, Any]]:
    """Coerce JSON-loaded rows using field types in schema so pa.Table.from_pylist can materialize them."""
    struct_type = pa.struct(list(schema))
    return [_coerce_json_rows_for_arrow(row, struct_type) for row in rows]


def _json_default(obj: object) -> str:
    """Serialize Python date and bytes values from Arrow tables for comparison with expected JSON."""
    if isinstance(obj, dt.date):
        return obj.isoformat()
    if isinstance(obj, bytes):
        return obj.decode("utf-8")
    msg = f"Object of type {type(obj)} is not JSON serializable"
    raise TypeError(msg)


def _strip_none(value: object) -> object:
    """Recursively remove keys with None values from dicts to compare Arrow-materialized rows with sparse expected.json."""
    if isinstance(value, dict):
        return {k: _strip_none(v) for k, v in value.items() if v is not None}
    if isinstance(value, list):
        return [_strip_none(item) for item in value]
    return value


def _extract_matching_elements(xml_path: Path, xml_tag: str) -> list[etree._Element]:
    """Parse xml_path and return all elements matching xml_tag at the root or as a direct child of root."""
    local_tag = xml_tag.rsplit("}", maxsplit=1)[-1] if "}" in xml_tag else xml_tag
    tree = etree.parse(str(xml_path))
    root = tree.getroot()
    root_local = root.tag.rsplit("}", maxsplit=1)[-1] if "}" in root.tag else root.tag
    if root_local == local_tag:
        return [root]
    return [
        child
        for child in root
        if isinstance(child.tag, str)
        and (child.tag.rsplit("}", maxsplit=1)[-1] if "}" in child.tag else child.tag) == local_tag
    ]


SPEC_CASES: list[ParameterSet] = [
    pytest.param("01_scalars_and_simple_types", "record", id="01-scalars-and-simple-types"),
    pytest.param("02_namespaces_and_nillable", "entry", id="02-namespaces-and-nillable"),
    pytest.param("03_mixed_content_raw_text", "article", id="03-mixed-content-raw-text"),
    pytest.param("04_recursion_depth", "record", id="04-recursion-depth"),
    pytest.param("05_groups_and_compositors", "item", id="05-groups-and-compositors"),
    pytest.param("06_context_dependent_paths", "root", id="06-context-dependent-paths"),
    pytest.param("07_multi_root_shared_type", "entry", id="07-multi-root-shared-type"),
    pytest.param("08_nullability_and_empty_tags", "entry", id="08-nullability-and-empty-tags"),
    pytest.param("09_empty_marker_elements", "book", id="09-empty-marker-elements"),
    pytest.param("10_combinations", "record", id="10-combinations"),
]

REAL_WORLD_SCHEMAS: list[ParameterSet] = [
    pytest.param("uniref", "uniprot/uniref/uniref.xsd", "entry", id="uniref"),
    pytest.param("uniprot", "uniprot/uniprot_kb/uniprot.xsd", "entry", id="uniprot"),
    pytest.param("pubmed_article", "pubmed/pubmed.xsd", "PubmedArticle", id="pubmed-article"),
    pytest.param("pubmed_book_article", "pubmed/pubmed.xsd", "PubmedBookArticle", id="pubmed-book-article"),
]

REAL_WORLD_EXAMPLES: list[ParameterSet] = [
    pytest.param("uniref", "uniprot/uniref/uniref.xsd", "entry", "minimal", id="uniref-minimal"),
    pytest.param("uniref", "uniprot/uniref/uniref.xsd", "entry", "maximal", id="uniref-maximal"),
    pytest.param("uniprot", "uniprot/uniprot_kb/uniprot.xsd", "entry", "minimal", id="uniprot-minimal"),
    pytest.param("uniprot", "uniprot/uniprot_kb/uniprot.xsd", "entry", "maximal", id="uniprot-maximal"),
    pytest.param("pubmed_article", "pubmed/pubmed.xsd", "PubmedArticle", "minimal", id="pubmed-article-minimal"),
    pytest.param("pubmed_article", "pubmed/pubmed.xsd", "PubmedArticle", "maximal", id="pubmed-article-maximal"),
    pytest.param(
        "pubmed_book_article",
        "pubmed/pubmed.xsd",
        "PubmedBookArticle",
        "minimal",
        id="pubmed-book-article-minimal",
    ),
    pytest.param(
        "pubmed_book_article",
        "pubmed/pubmed.xsd",
        "PubmedBookArticle",
        "maximal",
        id="pubmed-book-article-maximal",
    ),
]


@pytest.mark.parametrize(("case_dir", "xml_tag"), SPEC_CASES)
def test_spec_fixtures_pass_xml_and_expected_schema_validate(
    test_data_dir: Path,
    case_dir: str,
    xml_tag: str,
) -> None:
    """Every specification data.xml validates against schema.xsd and expected.json materializes with expected_schema.py."""
    fixture_dir = test_data_dir / "xsd" / "spec" / case_dir
    schema = load_schema(fixture_dir / "schema.xsd")
    schema.validate(str(fixture_dir / "data.xml"))

    elements = _extract_matching_elements(fixture_dir / "data.xml", xml_tag)
    expected = json.loads((fixture_dir / "expected.json").read_text(encoding="utf-8"))
    assert len(elements) == len(expected)

    expected_pa_schema = _load_expected_schema(fixture_dir / "expected_schema.py")
    coerced_rows = _coerce_rows_to_schema(expected, expected_pa_schema)
    table = pa.Table.from_pylist(coerced_rows, schema=expected_pa_schema)
    table.validate(full=True)
    assert table.num_rows == len(expected)

    dumped_rows = json.loads(json.dumps(_strip_none(table.to_pylist()), default=_json_default))
    assert dumped_rows == expected

    buf = io.BytesIO()
    pq.write_table(table, buf)
    buf.seek(0)
    roundtrip_table = pq.read_table(buf)
    roundtrip_rows = json.loads(json.dumps(_strip_none(roundtrip_table.to_pylist()), default=_json_default))
    assert roundtrip_rows == expected


@pytest.mark.parametrize(
    ("example_dir", "xsd_rel_path", "xml_tag", "variant"),
    REAL_WORLD_EXAMPLES,
)
def test_real_world_examples_pass_xml_and_expected_schema_validate(
    test_data_dir: Path,
    example_dir: str,
    xsd_rel_path: str,
    xml_tag: str,
    variant: str,
) -> None:
    """Every minimal and maximal real-world XML validates against its XSD and materializes with expected_schema.py."""
    schema = load_schema(test_data_dir / xsd_rel_path)
    example_base = test_data_dir / "xsd" / "examples" / example_dir
    xml_path = example_base / f"{variant}.xml"
    expected_path = example_base / f"{variant}_expected.json"

    schema.validate(str(xml_path))
    elements = _extract_matching_elements(xml_path, xml_tag)
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    assert len(elements) == 1
    assert len(expected) == 1

    expected_pa_schema = _load_expected_schema(example_base / "expected_schema.py")
    coerced_rows = _coerce_rows_to_schema(expected, expected_pa_schema)
    table = pa.Table.from_pylist(coerced_rows, schema=expected_pa_schema)
    table.validate(full=True)
    assert table.num_rows == 1

    dumped_rows = json.loads(json.dumps(_strip_none(table.to_pylist()), default=_json_default))
    assert dumped_rows == expected

    buf = io.BytesIO()
    pq.write_table(table, buf)
    buf.seek(0)
    roundtrip_table = pq.read_table(buf)
    roundtrip_rows = json.loads(json.dumps(_strip_none(roundtrip_table.to_pylist()), default=_json_default))
    assert roundtrip_rows == expected


@pytest.mark.parametrize(
    "example_dir",
    [
        pytest.param("uniref", id="uniref-maximal"),
        pytest.param("uniprot", id="uniprot-maximal"),
        pytest.param("pubmed_article", id="pubmed-article-maximal"),
        pytest.param("pubmed_book_article", id="pubmed-book-article-maximal"),
    ],
)
def test_real_world_maximal_examples_pass_repeated_elements_have_distinct_content(
    test_data_dir: Path,
    example_dir: str,
) -> None:
    """Every repeated element in a maximal XML example contains distinct data."""
    xml_path = test_data_dir / "xsd" / "examples" / example_dir / "maximal.xml"
    root = etree.parse(str(xml_path)).getroot()

    by_tag: dict[str, list[str]] = defaultdict(list)
    for elem in root.iter():
        if elem is root or not isinstance(elem.tag, str):
            continue
        clone = copy.deepcopy(elem)
        for node in clone.iter():
            if node.text is not None:
                node.text = node.text.strip()
            if node.tail is not None:
                node.tail = node.tail.strip()
        local = elem.tag.rsplit("}", maxsplit=1)[-1] if "}" in elem.tag else elem.tag
        by_tag[local].append(etree.tostring(clone, encoding="unicode"))

    duplicates = {tag: sigs for tag, sigs in by_tag.items() if len(sigs) > 1 and len(set(sigs)) < len(sigs)}
    assert not duplicates


@pytest.mark.parametrize(("case_dir", "xml_tag"), SPEC_CASES)
def test_generate_pyarrow_schema_pass_spec_cases(
    test_data_dir: Path,
    case_dir: str,
    xml_tag: str,
) -> None:
    """Verify XSD specification schemas generate the exact expected PyArrow schema in expected_schema.py."""
    fixture_dir = test_data_dir / "xsd" / "spec" / case_dir
    schema = load_schema(fixture_dir / "schema.xsd")
    pa_schema = generate_pyarrow_schema(schema, xml_tag, max_recursion_depth=2)
    expected_schema = _load_expected_schema(fixture_dir / "expected_schema.py")
    assert pa_schema.equals(expected_schema, check_metadata=True)


def test_generate_pyarrow_schema_pass_clark_notation_tag(test_data_dir: Path) -> None:
    """Verify Clark-notation {namespace}local_name resolves identically to local_name."""
    fixture_dir = test_data_dir / "xsd" / "spec" / "02_namespaces_and_nillable"
    schema = load_schema(fixture_dir / "schema.xsd")
    pa_clark = generate_pyarrow_schema(schema, "{http://example.org/ns}entry")
    expected_schema = _load_expected_schema(fixture_dir / "expected_schema.py")
    assert pa_clark.equals(expected_schema, check_metadata=True)


def test_generate_pyarrow_schema_pass_04_recursion_depth_1(test_data_dir: Path) -> None:
    """Verify structural recursion with max_recursion_depth=1 matches EXPECTED_SCHEMA_DEPTH_1."""
    fixture_dir = test_data_dir / "xsd" / "spec" / "04_recursion_depth"
    schema = load_schema(fixture_dir / "schema.xsd")
    schema_d1 = generate_pyarrow_schema(schema, "record", max_recursion_depth=1)
    expected_d1 = _load_expected_schema(fixture_dir / "expected_schema.py", "EXPECTED_SCHEMA_DEPTH_1")
    assert schema_d1.equals(expected_d1, check_metadata=True)


def test_generate_pyarrow_schema_pass_mutual_and_self_referential_xsds(test_data_dir: Path) -> None:
    """Verify existing self_referential.xsd and mutual_recursion.xsd schemas terminate at max_recursion_depth."""
    self_ref_schema = load_schema(test_data_dir / "xsd" / "self_referential.xsd")
    pa_self = generate_pyarrow_schema(self_ref_schema, "node", max_recursion_depth=2)
    expected_self = pa.schema(
        [
            _list_field(
                "child",
                pa.struct(
                    [
                        _list_field("child", pa.string(), nullable=True),
                    ]
                ),
                nullable=True,
            ),
        ]
    )
    assert pa_self.equals(expected_self, check_metadata=True)

    mutual_schema = load_schema(test_data_dir / "xsd" / "mutual_recursion.xsd")
    pa_mutual = generate_pyarrow_schema(mutual_schema, "a", max_recursion_depth=1)
    expected_mutual = pa.schema(
        [
            _list_field(
                "b",
                pa.struct(
                    [
                        pa.field("a", pa.string(), nullable=True),
                    ]
                ),
                nullable=True,
            ),
        ]
    )
    assert pa_mutual.equals(expected_mutual, check_metadata=True)


@pytest.mark.parametrize(
    ("example_dir", "xsd_rel_path", "xml_tag"),
    REAL_WORLD_SCHEMAS,
)
def test_generate_pyarrow_schema_pass_real_world_schemas(
    test_data_dir: Path,
    example_dir: str,
    xsd_rel_path: str,
    xml_tag: str,
) -> None:
    """Verify PyArrow schema generation on UniRef, UniProtKB, and PubMed XSDs matches expected_schema.py."""
    schema = load_schema(test_data_dir / xsd_rel_path)
    pa_schema = generate_pyarrow_schema(schema, xml_tag, max_recursion_depth=2)
    expected_schema = _load_expected_schema(test_data_dir / "xsd" / "examples" / example_dir / "expected_schema.py")
    assert pa_schema.equals(expected_schema, check_metadata=True)


@pytest.mark.parametrize(
    ("error_xsd_filename", "xml_tag", "expected_match"),
    [
        pytest.param("ambiguous_xml_tag.xsd", "entry", "Ambiguous", id="E2-ambiguous-xml-tag"),
        pytest.param("scalar_top_level.xsd", "leaf", "complex type", id="E3-scalar-top-level"),
        pytest.param("empty_marker_top_level.xsd", "marker", "complex type", id="E3-empty-marker-top-level"),
        pytest.param("wildcard_any_element.xsd", "entry", "xs:any", id="E4-wildcard-any-element"),
        pytest.param("wildcard_any_attribute.xsd", "entry", "xs:anyAttribute", id="E5-wildcard-any-attribute"),
        pytest.param("abstract_element.xsd", "entry", "abstract", id="E6-abstract-element"),
        pytest.param("abstract_type.xsd", "entry", "abstract", id="E6-abstract-type"),
        pytest.param("substitution_group.xsd", "entry", "substitutionGroup", id="E7-substitution-group"),
        pytest.param("attr_child_name_collision.xsd", "entry", "collision", id="E8-attr-child-name-collision"),
    ],
)
def test_generate_pyarrow_schema_fail_unsupported_or_invalid_xsd(
    test_data_dir: Path,
    error_xsd_filename: str,
    xml_tag: str,
    expected_match: str,
) -> None:
    """Verify unsupported or ambiguous XSD features raise a descriptive ValueError."""
    schema = load_schema(test_data_dir / "xsd" / "spec" / "errors" / error_xsd_filename)
    with pytest.raises(ValueError, match=expected_match):
        generate_pyarrow_schema(schema, xml_tag)


@pytest.mark.parametrize(
    ("xsd_filename", "xml_tag", "expected_match"),
    [
        pytest.param("empty_schema.xsd", "root", "not found", id="E1-empty-schema"),
        pytest.param("no_children.xsd", "leaf", "complex type", id="E3-no-children-scalar-root"),
        pytest.param("cdata_variety.xsd", "root", "xs:anyAttribute", id="E5-cdata-variety-has-anyAttribute"),
        pytest.param("combinations.xsd", "record", "xs:anyAttribute", id="E5-legacy-combinations-has-anyAttribute"),
    ],
)
def test_generate_pyarrow_schema_fail_existing_invalid_schemas(
    test_data_dir: Path,
    xsd_filename: str,
    xml_tag: str,
    expected_match: str,
) -> None:
    """Verify existing XSDs with scalar roots, empty schemas, or xs:anyAttribute raise ValueError."""
    schema = load_schema(test_data_dir / "xsd" / xsd_filename)
    with pytest.raises(ValueError, match=expected_match):
        generate_pyarrow_schema(schema, xml_tag)


def test_generate_pyarrow_schema_fail_nonexistent_xml_tag(test_data_dir: Path) -> None:
    """Verify requesting an xml_tag that does not exist in the XSD raises ValueError."""
    schema = load_schema(test_data_dir / "uniprot" / "uniref" / "uniref.xsd")
    with pytest.raises(ValueError, match="not found"):
        generate_pyarrow_schema(schema, "nonExistentElement")


@pytest.mark.parametrize(
    "forbidden_xsd_type",
    [
        pytest.param("xs:date", id="xs-date"),
        pytest.param("xs:dateTime", id="xs-dateTime"),
        pytest.param("xs:hexBinary", id="xs-hexBinary"),
        pytest.param("xs:base64Binary", id="xs-base64Binary"),
    ],
)
def test_generate_pyarrow_schema_fail_non_json_serializable_type_under_recursion(
    forbidden_xsd_type: str,
) -> None:
    """Verify xs:date, xs:dateTime, xs:hexBinary, and xs:base64Binary under a recursive element raise ValueError."""
    xsd_with_child = f"""<?xml version="1.0" encoding="UTF-8"?>
    <xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
      <xs:complexType name="NodeType">
        <xs:sequence>
          <xs:element name="payload" type="{forbidden_xsd_type}" minOccurs="0"/>
          <xs:element name="node" type="NodeType" minOccurs="0" maxOccurs="unbounded"/>
        </xs:sequence>
      </xs:complexType>
      <xs:element name="root" type="NodeType"/>
    </xs:schema>
    """
    schema_child = xmlschema.XMLSchema(xsd_with_child)
    with pytest.raises(ValueError, match="recursive type"):
        generate_pyarrow_schema(schema_child, "root")

    xsd_with_attr = f"""<?xml version="1.0" encoding="UTF-8"?>
    <xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
      <xs:complexType name="NodeType">
        <xs:sequence>
          <xs:element name="node" type="NodeType" minOccurs="0" maxOccurs="unbounded"/>
        </xs:sequence>
        <xs:attribute name="attrPayload" type="{forbidden_xsd_type}"/>
      </xs:complexType>
      <xs:element name="root" type="NodeType"/>
    </xs:schema>
    """
    schema_attr = xmlschema.XMLSchema(xsd_with_attr)
    with pytest.raises(ValueError, match="recursive type"):
        generate_pyarrow_schema(schema_attr, "root")


@pytest.mark.parametrize(("case_dir", "xml_tag"), SPEC_CASES)
def test_convert_xml_element_and_arrow_table_pass_spec_fixtures(
    test_data_dir: Path,
    case_dir: str,
    xml_tag: str,
) -> None:
    """Verify XML elements convert to expected JSON and build a valid PyArrow Table and Parquet file."""
    fixture_dir = test_data_dir / "xsd" / "spec" / case_dir
    schema = load_schema(fixture_dir / "schema.xsd")
    pa_schema = generate_pyarrow_schema(schema, xml_tag, max_recursion_depth=2)

    elements = _extract_matching_elements(fixture_dir / "data.xml", xml_tag)
    rows = [convert_xml_element(elem, schema, xml_tag, max_recursion_depth=2) for elem in elements]

    expected = json.loads((fixture_dir / "expected.json").read_text(encoding="utf-8"))
    serialized_rows = json.loads(json.dumps(rows, default=_json_default))
    assert serialized_rows == expected

    table = pa.Table.from_pylist(rows, schema=pa_schema)
    table.validate(full=True)
    assert table.num_rows == len(expected)

    dumped_rows = json.loads(json.dumps(_strip_none(table.to_pylist()), default=_json_default))
    assert dumped_rows == expected

    buf = io.BytesIO()
    pq.write_table(table, buf)
    buf.seek(0)
    roundtrip_table = pq.read_table(buf)
    roundtrip_rows = json.loads(json.dumps(_strip_none(roundtrip_table.to_pylist()), default=_json_default))
    assert roundtrip_rows == expected


def test_convert_xml_element_and_arrow_table_pass_recursion_depth_1(test_data_dir: Path) -> None:
    """Verify max_recursion_depth=1 truncates recursive structures one level earlier than depth=2."""
    fixture_dir = test_data_dir / "xsd" / "spec" / "04_recursion_depth"
    schema = load_schema(fixture_dir / "schema.xsd")
    pa_schema = generate_pyarrow_schema(schema, "record", max_recursion_depth=1)

    elements = _extract_matching_elements(fixture_dir / "data.xml", "record")
    rows = [convert_xml_element(elem, schema, "record", max_recursion_depth=1) for elem in elements]

    expected = json.loads((fixture_dir / "expected_depth_1.json").read_text(encoding="utf-8"))
    serialized_rows = json.loads(json.dumps(rows, default=_json_default))
    assert serialized_rows == expected

    table = pa.Table.from_pylist(rows, schema=pa_schema)
    table.validate(full=True)
    dumped_rows = json.loads(json.dumps(_strip_none(table.to_pylist()), default=_json_default))
    assert dumped_rows == expected

    buf = io.BytesIO()
    pq.write_table(table, buf)
    buf.seek(0)
    roundtrip_table = pq.read_table(buf)
    roundtrip_rows = json.loads(json.dumps(_strip_none(roundtrip_table.to_pylist()), default=_json_default))
    assert roundtrip_rows == expected


@pytest.mark.parametrize(
    ("example_dir", "xsd_rel_path", "xml_tag", "variant"),
    REAL_WORLD_EXAMPLES,
)
def test_convert_xml_element_and_arrow_table_pass_real_world_examples(
    test_data_dir: Path,
    example_dir: str,
    xsd_rel_path: str,
    xml_tag: str,
    variant: str,
) -> None:
    """Verify minimal and maximal UniRef, UniProt, and PubMed entries convert to expected JSON and Arrow tables."""
    schema = load_schema(test_data_dir / xsd_rel_path)
    xml_path = test_data_dir / "xsd" / "examples" / example_dir / f"{variant}.xml"
    expected_path = test_data_dir / "xsd" / "examples" / example_dir / f"{variant}_expected.json"

    pa_schema = generate_pyarrow_schema(schema, xml_tag, max_recursion_depth=2)
    elements = _extract_matching_elements(xml_path, xml_tag)
    rows = [convert_xml_element(elem, schema, xml_tag, max_recursion_depth=2) for elem in elements]

    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    serialized_rows = json.loads(json.dumps(rows, default=_json_default))
    assert serialized_rows == expected

    table = pa.Table.from_pylist(rows, schema=pa_schema)
    table.validate(full=True)
    assert table.num_rows == 1

    dumped_rows = json.loads(json.dumps(_strip_none(table.to_pylist()), default=_json_default))
    assert dumped_rows == expected

    buf = io.BytesIO()
    pq.write_table(table, buf)
    buf.seek(0)
    roundtrip_table = pq.read_table(buf)
    assert roundtrip_table.num_rows == 1
    roundtrip_rows = json.loads(json.dumps(_strip_none(roundtrip_table.to_pylist()), default=_json_default))
    assert roundtrip_rows == expected
