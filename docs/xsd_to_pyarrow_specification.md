# XSD-to-PyArrow Specification

This document specifies how an XML Schema (`.xsd`) and target `xml_tag` determine:
1. The **PyArrow `pa.Schema`** for the table yielded by the pipeline.
2. The **deterministic Python / JSON row dictionaries** produced by schema-driven XML element conversion (`convert_xml_element`) before constructing `pa.Table.from_pylist(rows, schema=pa_schema)`.
3. The **error conditions** raised when an XSD or `xml_tag` uses unsupported or ambiguous constructs.

---

## 1. Summary of Conversion Rules

### 1.1 Top-Level `xml_tag` Resolution & Unwrapping
- **Resolution:** Search the XSD for `xml_tag` starting at depth 0 (global root elements) and depth 1 (direct children of global root elements), deduplicating matches by `id(element.type)`.
  - If **0** elements match (even after falling back to deeper elements), raise an error.
  - If **>1** distinct `XsdType` definitions match at the shallowest matching depth, raise an ambiguity error.
- **Unwrapping:** Each streamed `<{xml_tag}>` element becomes one row in the output `pa.Table`.
  - `xml_tag` must have a `complexType` that maps to a `pa.StructType` (at least one attribute, `#text`, or child element field). The outer `xml_tag` key is unwrapped so that the top-level columns of `pa.Schema` are the struct fields of `xml_tag`.
  - If `xml_tag` maps to a scalar, list, or empty marker type with no attributes and no struct fields, schema generation raises a `ValueError` (see Error E3 in Section 3).

### 1.2 Naming & Field Ordering
- **Attributes:** Prefixed with `_` (e.g., `id` $\rightarrow$ `_id`). Namespace declarations (`xmlns`, `xmlns:*`) and `xsi:*` attributes (`xsi:nil`, `xsi:schemaLocation`, `xsi:noNamespaceSchemaLocation`) are stripped.
- **Text content (`#text`):** When a complex type has both attributes and text (`simpleContent` with attributes, or `mixed="true"` with attributes), text is stored under `#text`.
- **Struct Field Order:**
  1. Declared attributes (`_attr`), in XSD declaration order (base type attributes before extension attributes).
  2. `#text` (if applicable).
  3. Child elements, in XSD content-model traversal order (deduplicated on first appearance if the same element name appears in multiple branches of a `<xs:choice>` or `<xs:sequence>`).

### 1.3 Element Classification (Shape)
| XSD Element Type | Cardinality (`overall_max_occurs`) | Python / JSON Value | PyArrow DataType |
| :--- | :--- | :--- | :--- |
| **`xs:simpleType`** (atomic) | `1` | Scalar (`str`, `int`, `float`, `bool`, `date`, `bytes`) or `None` | Scalar `pa.DataType` |
| **`xs:simpleType`** (atomic) | `> 1` or `unbounded` | `list[Scalar \| None]` or `None` (when absent) | `pa.list_(Scalar)` |
| **`xs:simpleType`** (`xs:list`) | `1` | `list[ItemScalar]` or `None` | `pa.list_(ItemScalar)` |
| **`xs:simpleType`** (`xs:list`) | `> 1` or `unbounded` | `list[list[ItemScalar] \| None]` or `None` | `pa.list_(pa.list_(ItemScalar))` |
| **`complexType` (`simpleContent`)** without attributes | `1` / `> 1` | Same as underlying `simpleType` | Same as underlying `simpleType` |
| **`complexType` (`simpleContent`)** with attributes | `1` | `{"_attr": ..., "#text": Scalar}` or `None` | `pa.struct([("_attr", ...), ("#text", Scalar)])` |
| **`complexType` (`mixed="true"`)** without attributes | `1` | Raw inner XML string `str` or `None` | `pa.string()` |
| **`complexType` (`mixed="true"`)** with attributes | `1` | `{"_attr": ..., "#text": str}` or `None` | `pa.struct([("_attr", ...), ("#text", pa.string())])` |
| **`complexType` (attribute-only)** | `1` | `{"_attr": ...}` or `None` | `pa.struct([("_attr", ...)])` |
| **`complexType` (element-only)** | `1` | `{"_attr": ..., "child": ...}` or `None` | `pa.struct([...])` |
| **`complexType` (empty: no attrs, no children, not mixed)** | `1` / `> 1` | `""` (when present) or `None` (when absent) / `["", ...]` | `pa.string()` / `pa.list_(pa.string())` |

*(For any `complexType` row above with `overall_max_occurs > 1` or `unbounded`, the Python value is `list[...]` and the PyArrow type is `pa.list_(...)`. When structural recursion exceeds `max_recursion_depth`, a single recursive element becomes `pa.string()` and a repeated recursive element (`overall_max_occurs > 1`) becomes `pa.list_(pa.string())`, where each truncated item is serialized as a JSON string.)*

### 1.4 Field Nullability (`nullable=False` vs `nullable=True`) & Empty Tag Coercion
- **Attributes (`_attr`):**
  - `nullable=False` when `use="required"` (and the enclosing element is not a `nillable="true"` element that can be nilled without attributes).
  - `nullable=True` when `use="optional"` (default).
- **Child Elements:**
  - `nullable=False` when **all** of the following hold:
    1. The element particle has `minOccurs >= 1`.
    2. The element declaration has `nillable="false"` (default).
    3. Every enclosing compositor/group (`<xs:sequence>`, `<xs:all>`, `<xs:group>`) up to the parent `complexType` has `minOccurs >= 1`, and the element is not conditional inside an `<xs:choice>` branch (unless required in every branch of the choice).
    4. The parent element is not `nillable="true"` (since a parent with `xsi:nil="true"` and attributes, e.g. `<comment type="function" xsi:nil="true"/>`, produces a struct `{"_type": "function"}` whose child elements are `None`).
  - `nullable=True` in all other cases.
- **`#text` Field Inside a Struct with Attributes:**
  - `nullable=True` for `mixed="true"` types, `nillable="true"` elements, and string-based `simpleContent` types (where `<tag attr="x"/>` has no `#text` key).
  - `nullable=False` for non-nillable `simpleContent` types with non-string atomic base types (e.g., `xs:int`, `xs:decimal`), where valid XML always provides lexical text.
- **Empty Standalone Tag Coercion (`<tag/>` $\rightarrow$ `""`):**
  - When a standalone element mapped to `pa.string()` (atomic string `simpleType`, string `simpleContent` without attributes, `mixed="true"` without attributes, or empty `complexType` marker) is **present** in the XML as `<tag/>` or `<tag></tag>` (without `xsi:nil="true"`), its value is coerced to `""` so that present empty elements satisfy `nullable=False` and are distinguishable from omitted optional elements (`None`).

### 1.5 Built-In XSD Scalar Type Mapping & Field Metadata
| XSD Built-In Type(s) | Python Type (Coerced) | PyArrow `DataType` |
| :--- | :--- | :--- |
| `xs:string`, `xs:normalizedString`, `xs:token`, `xs:language`, `xs:Name`, `xs:NCName`, `xs:NMTOKEN`, `xs:ID`, `xs:IDREF`, `xs:ENTITY`, `xs:anyURI`, `xs:QName`, `xs:NOTATION` | `str` | `pa.string()` |
| `xs:boolean` (`"true"`/`"1"` $\rightarrow$ `True`, `"false"`/`"0"` $\rightarrow$ `False`) | `bool` | `pa.bool_()` |
| `xs:byte` | `int` | `pa.int8()` |
| `xs:short` | `int` | `pa.int16()` |
| `xs:int` | `int` | `pa.int32()` |
| `xs:long`, `xs:integer`, `xs:positiveInteger`, `xs:nonNegativeInteger`, `xs:negativeInteger`, `xs:nonPositiveInteger` | `int` | `pa.int64()` |
| `xs:unsignedByte` | `int` | `pa.uint8()` |
| `xs:unsignedShort` | `int` | `pa.uint16()` |
| `xs:unsignedInt` | `int` | `pa.uint32()` |
| `xs:unsignedLong` | `int` | `pa.uint64()` |
| `xs:float`, `xs:double`, `xs:decimal` | `float` | `pa.float64()` |
| `xs:date` | `datetime.date` | `pa.date32()` |
| `xs:dateTime`, `xs:time`, `xs:duration`, `xs:gYear`, `xs:gYearMonth`, `xs:gMonthDay`, `xs:gDay`, `xs:gMonth` | `str` | `pa.string()` |
| `xs:base64Binary` | `bytes` (`base64.b64decode`) | `pa.binary()` |
| `xs:hexBinary` | `bytes` (`bytes.fromhex`) | `pa.binary()` |
| `xs:NMTOKENS`, `xs:IDREFS`, `xs:ENTITIES` | `list[str]` (`.split()`) | `pa.list_(pa.string())` |

- **`xs:restriction`:** Inherits the primitive/built-in base type mapping above.
- **`xs:list`:** Maps to `pa.list_(item_type)` and coerces by splitting whitespace (`value.split()`) and coercing each item.
- **`xs:union`:** If all `memberTypes` map to the same `pa.DataType`, uses that type; otherwise falls back to `pa.string()` (`str`).
- **`pa.Field.metadata`:**
  - **`"xsd_type"`:** Recorded on every leaf scalar `pa.Field` (including `#text` and list `"item"` fields when list elements are scalars) with the prefixed built-in XSD type name (e.g., `"xs:string"`, `"xs:int"`, `"xs:float"`, `"xs:double"`, `"xs:decimal"`, or `"xs:date | xs:gYearMonth | xs:gYear"` for unions). Omitted on `pa.struct(...)`, outer `pa.list_(...)`, empty-marker `pa.string()`, and recursion-truncated JSON `pa.string()` fields.
  - **`"comment"`:** Extracted from `<xs:annotation><xs:documentation>` on the element or attribute declaration (falling back to its referenced `type`'s `<xs:annotation><xs:documentation>` if the declaration itself has none), whitespace-normalized with `" ".join(doc.split())`.

---

## 2. Specification Examples

### Example 1: Scalar Types, Restrictions, Lists (`xs:list`), and Unions (`xs:union`)
Covers built-in numeric/boolean/date/binary scalars, `xs:restriction` enumerations, `xs:list` (`intListType` from [`uniprot.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/uniprot/uniprot_kb/uniprot.xsd#L1122-L1124)), and `xs:union` (`date | gYearMonth | gYear` from [`uniprot.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/uniprot/uniprot_kb/uniprot.xsd#L341-L345)).

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:simpleType name="intListType">
        <xs:list itemType="xs:int"/>
    </xs:simpleType>

    <xs:simpleType name="flexibleDateType">
        <xs:union memberTypes="xs:date xs:gYearMonth xs:gYear"/>
    </xs:simpleType>

    <xs:simpleType name="intUnionType">
        <xs:union memberTypes="xs:int xs:short"/>
    </xs:simpleType>

    <xs:simpleType name="statusEnum">
        <xs:restriction base="xs:string">
            <xs:enumeration value="active"/>
            <xs:enumeration value="archived"/>
        </xs:restriction>
    </xs:simpleType>

    <xs:element name="dataset">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="record" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="flag" type="xs:boolean"/>
                            <xs:element name="smallInt" type="xs:byte" minOccurs="0"/>
                            <xs:element name="medInt" type="xs:short" minOccurs="0"/>
                            <xs:element name="stdInt" type="xs:int" minOccurs="0"/>
                            <xs:element name="bigInt" type="xs:integer" minOccurs="0"/>
                            <xs:element name="ulong" type="xs:unsignedLong" minOccurs="0"/>
                            <xs:element name="f32" type="xs:float" minOccurs="0"/>
                            <xs:element name="f64" type="xs:double" minOccurs="0"/>
                            <xs:element name="dec" type="xs:decimal" minOccurs="0"/>
                            <xs:element name="created" type="xs:date" minOccurs="0"/>
                            <xs:element name="published" type="flexibleDateType" minOccurs="0"/>
                            <xs:element name="code" type="intUnionType" minOccurs="0"/>
                            <xs:element name="scores" type="intListType" minOccurs="0"/>
                            <xs:element name="rawHex" type="xs:hexBinary" minOccurs="0"/>
                            <xs:element name="rawB64" type="xs:base64Binary" minOccurs="0"/>
                        </xs:sequence>
                        <xs:attribute name="id" type="xs:ID" use="required"/>
                        <xs:attribute name="status" type="statusEnum" use="required"/>
                        <xs:attribute name="evidence" type="intListType" use="optional"/>
                        <xs:attribute name="version" type="xs:int" use="required"/>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="record"`)
```xml
<dataset>
    <record id="rec1" status="active" evidence="1 4 9" version="2">
        <flag>true</flag>
        <smallInt>-5</smallInt>
        <medInt>1024</medInt>
        <stdInt>100000</stdInt>
        <bigInt>999999999999</bigInt>
        <ulong>18446744073709551615</ulong>
        <f32>1.5</f32>
        <f64>2.718281828</f64>
        <dec>3.14</dec>
        <created>2025-06-15</created>
        <published>2025-06</published>
        <code >42</code>
        <scores>10 20 30</scores>
        <rawHex>48656c6c6f</rawHex>
        <rawB64>SGVsbG8=</rawB64>
    </record>
    <record id="rec2" status="archived" version="1">
        <flag>0</flag>
    </record>
</dataset>
```

#### Python / JSON Translation
```python
[
    {
        "_id": "rec1",
        "_status": "active",
        "_evidence": [1, 4, 9],
        "_version": 2,
        "flag": True,
        "smallInt": -5,
        "medInt": 1024,
        "stdInt": 100000,
        "bigInt": 999999999999,
        "ulong": 18446744073709551615,
        "f32": 1.5,
        "f64": 2.718281828,
        "dec": 3.14,
        "created": datetime.date(2025, 6, 15),
        "published": "2025-06",
        "code": "42",  # int + short -> different Arrow types (int32 vs int16) -> falls back to str
        "scores": [10, 20, 30],
        "rawHex": b"Hello",
        "rawB64": b"Hello",
    },
    {
        "_id": "rec2",
        "_status": "archived",
        "_version": 1,
        "flag": False,
    },
]
```

#### PyArrow Schema
```python
pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False, metadata={"xsd_type": "xs:ID"}),
        pa.field("_status", pa.string(), nullable=False, metadata={"xsd_type": "xs:string"}),
        pa.field(
            "_evidence",
            pa.list_(pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"})),
            nullable=True,
        ),
        pa.field("_version", pa.int32(), nullable=False, metadata={"xsd_type": "xs:int"}),
        pa.field("flag", pa.bool_(), nullable=False, metadata={"xsd_type": "xs:boolean"}),
        pa.field("smallInt", pa.int8(), nullable=True, metadata={"xsd_type": "xs:byte"}),
        pa.field("medInt", pa.int16(), nullable=True, metadata={"xsd_type": "xs:short"}),
        pa.field("stdInt", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"}),
        pa.field("bigInt", pa.int64(), nullable=True, metadata={"xsd_type": "xs:integer"}),
        pa.field("ulong", pa.uint64(), nullable=True, metadata={"xsd_type": "xs:unsignedLong"}),
        pa.field("f32", pa.float64(), nullable=True, metadata={"xsd_type": "xs:float"}),
        pa.field("f64", pa.float64(), nullable=True, metadata={"xsd_type": "xs:double"}),
        pa.field("dec", pa.float64(), nullable=True, metadata={"xsd_type": "xs:decimal"}),
        pa.field("created", pa.date32(), nullable=True, metadata={"xsd_type": "xs:date"}),
        pa.field("published", pa.string(), nullable=True, metadata={"xsd_type": "xs:date | xs:gYearMonth | xs:gYear"}),
        pa.field("code", pa.string(), nullable=True, metadata={"xsd_type": "xs:int | xs:short"}),
        pa.field(
            "scores",
            pa.list_(pa.field("item", pa.int32(), nullable=True, metadata={"xsd_type": "xs:int"})),
            nullable=True,
        ),
        pa.field("rawHex", pa.binary(), nullable=True, metadata={"xsd_type": "xs:hexBinary"}),
        pa.field("rawB64", pa.binary(), nullable=True, metadata={"xsd_type": "xs:base64Binary"}),
    ]
)
```
*(In Examples 2–9 below, `metadata={"xsd_type": ...}` on scalar fields is omitted for brevity.)*

---

### Example 2: Namespaces & `nillable="true"` (`xsi:nil="true"`)
Covers `targetNamespace`, `elementFormDefault="qualified"`, stripping `xmlns` / `xmlns:xsi`, and the 3 cases of `xsi:nil="true"`:
1. Nillable element with **no other attributes** $\rightarrow$ `None`.
2. Nillable complex element **with attributes** (like `<comment type="function" xsi:nil="true"/>` in [`uniprot.xsd` L47](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/uniprot/uniprot_kb/uniprot.xsd#L47)) $\rightarrow$ `{"_type": "function"}` (with `_xsi:nil` stripped and child fields `null`).
3. Nillable `simpleContent` element **with attributes** $\rightarrow$ `{"_lang": "en"}` (`#text` is `null`).

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema"
           targetNamespace="http://example.org/ns"
           xmlns="http://example.org/ns"
           elementFormDefault="qualified">
    <xs:complexType name="commentType">
        <xs:sequence>
            <xs:element name="text" type="xs:string"/>
        </xs:sequence>
        <xs:attribute name="type" type="xs:string" use="required"/>
    </xs:complexType>

    <xs:complexType name="labelType">
        <xs:simpleContent>
            <xs:extension base="xs:string">
                <xs:attribute name="lang" type="xs:string"/>
            </xs:extension>
        </xs:simpleContent>
    </xs:complexType>

    <xs:element name="root">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="entry" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="title" type="xs:string" nillable="true"/>
                            <xs:element name="label" type="labelType" nillable="true" minOccurs="0"/>
                            <xs:element name="comment" type="commentType" nillable="true" minOccurs="0" maxOccurs="unbounded"/>
                        </xs:sequence>
                        <xs:attribute name="id" type="xs:string" use="required"/>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="entry"`)
```xml
<root xmlns="http://example.org/ns" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <entry id="e1">
        <title xsi:nil="true"/>
        <label lang="en" xsi:nil="true"/>
        <comment type="function" xsi:nil="true"/>
        <comment type="subunit">
            <text>Homodimer</text>
        </comment>
    </entry>
    <entry id="e2">
        <title>Normal Title</title>
        <label xsi:nil="true"/>
    </entry>
</root>
```

#### Python / JSON Translation
```python
[
    {
        "_id": "e1",
        "title": None,
        "label": {"_lang": "en"},
        "comment": [
            {"_type": "function"},
            {"_type": "subunit", "text": "Homodimer"},
        ],
    },
    {
        "_id": "e2",
        "title": "Normal Title",
        "label": None,
    },
]
```

#### PyArrow Schema
```python
pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False),
        pa.field("title", pa.string(), nullable=True),
        pa.field(
            "label",
            pa.struct(
                [
                    pa.field("_lang", pa.string(), nullable=True),
                    pa.field("#text", pa.string(), nullable=True),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "comment",
            pa.list_(
                pa.struct(
                    [
                        pa.field("_type", pa.string(), nullable=False),
                        pa.field("text", pa.string(), nullable=True),
                    ]
                )
            ),
            nullable=True,
        ),
    ]
)
```

---

### Example 3: Mixed Content Kept as Raw Text
Covers `mixed="true"` elements with inline HTML/formatting tags (`<b>`, `<i>`, `<sup>`, `<sub>`, `<u>`), both with and without attributes on the mixed element itself (as in [`pubmed.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/pubmed/pubmed.xsd) `ArticleTitle`, `AbstractText`, and [`10_combinations/schema.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/10_combinations/schema.xsd) `ParaType`, `NoteType`, `RemarkType`).

- Child elements declared inside a `mixed="true"` complexType are **not** emitted as struct fields; instead, the mixed element's entire inner XML is captured as a raw text string.
- Character entities in text nodes (`&amp;`, `&lt;`, `&gt;`) remain XML-escaped in the captured raw text string so that literal `<` or `>` is distinguishable from inline HTML tags like `<b>`.
- If the `mixed="true"` type has **no attributes**, its PyArrow type is `pa.string()` (and an empty tag `<title/>` without `xsi:nil="true"` coerces to `""`).
- If the `mixed="true"` type **has attributes**, its PyArrow type is `pa.struct([("_attr", ...), ("#text", pa.string())])`.

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:group name="inlineFormatting">
        <xs:choice>
            <xs:element name="b" type="formattedText"/>
            <xs:element name="i" type="formattedText"/>
            <xs:element name="sup" type="formattedText"/>
        </xs:choice>
    </xs:group>

    <xs:complexType name="formattedText" mixed="true">
        <xs:group ref="inlineFormatting" minOccurs="0" maxOccurs="unbounded"/>
    </xs:complexType>

    <xs:complexType name="abstractTextType" mixed="true">
        <xs:group ref="inlineFormatting" minOccurs="0" maxOccurs="unbounded"/>
        <xs:attribute name="Label" type="xs:string"/>
    </xs:complexType>

    <xs:element name="articles">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="article" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="title" type="formattedText"/>
                            <xs:element name="abstractText" type="abstractTextType" minOccurs="0" maxOccurs="unbounded"/>
                        </xs:sequence>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="article"`)
```xml
<articles>
    <article>
        <title> Ca<sup>2+</sup> signaling in <i>E. <b>coli</b></i> cells </title>
        <abstractText Label="BACKGROUND">Studied <b>bold</b> &amp; <i>italic</i> text.</abstractText>
        <abstractText>Unlabelled paragraph with <b>tags</b>.</abstractText>
        <abstractText Label="EMPTY"/>
        <abstractText/>
    </article>
    <article>
        <title/>
    </article>
</articles>
```

#### Python / JSON Translation
```python
[
    {
        "title": "Ca<sup>2+</sup> signaling in <i>E. <b>coli</b></i> cells",
        "abstractText": [
            {"_Label": "BACKGROUND", "#text": "Studied <b>bold</b> &amp; <i>italic</i> text."},
            {"#text": "Unlabelled paragraph with <b>tags</b>."},
            {"_Label": "EMPTY"},
            None,
        ],
    },
    {
        "title": "",
    },
]
```

#### PyArrow Schema
```python
pa.schema(
    [
        pa.field("title", pa.string(), nullable=False),
        pa.field(
            "abstractText",
            pa.list_(
                pa.struct(
                    [
                        pa.field("_Label", pa.string(), nullable=True),
                        pa.field("#text", pa.string(), nullable=True),
                    ]
                )
            ),
            nullable=True,
        ),
    ]
)
```

---

### Example 4: Structural Recursion with Configurable `max_recursion_depth`
Covers both:
1. **Type self-recursion** (`TreeType` $\rightarrow$ `branch: TreeType`)
2. **Element-ref self-recursion** (`Section` $\rightarrow$ `<xs:element ref="Section"/>` as in [`pubmed.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/pubmed/pubmed.xsd))
3. **Mutual recursion** (`NodeA` $\rightarrow$ `NodeB` $\rightarrow$ `NodeA`)

With `max_recursion_depth = 2`:
- Occurrence 1 (top-level `Section` or `tree`) = depth 1 (`pa.struct`)
- Occurrence 2 (child `Section` or `branch`) = depth 2 (`pa.struct`)
- Occurrence 3 (grandchild `Section` or `branch` beyond `max_recursion_depth=2`) = flattened to JSON string (`pa.string()` when `overall_max_occurs == 1`, or `pa.list_(pa.string())` when `overall_max_occurs > 1`).

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:complexType name="TreeType">
        <xs:sequence>
            <xs:element name="branch" type="TreeType" minOccurs="0" maxOccurs="unbounded"/>
        </xs:sequence>
        <xs:attribute name="label" type="xs:string"/>
    </xs:complexType>

    <xs:element name="doc">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="record" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="tree" type="TreeType" minOccurs="0"/>
                        </xs:sequence>
                        <xs:attribute name="id" type="xs:string" use="required"/>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="record"`, `max_recursion_depth=2`)
```xml
<doc>
    <record id="r1">
        <tree label="L1">
            <branch label="L2">
                <branch label="L3">
                    <branch label="L4"/>
                </branch>
                <branch label="L3_sibling"/>
            </branch>
        </tree>
    </record>
</doc>
```

#### Python / JSON Translation (`max_recursion_depth=2`)
```python
[
    {
        "_id": "r1",
        "tree": {
            "_label": "L1",
            "branch": [
                {
                    "_label": "L2",
                    "branch": [
                        '{"_label": "L3", "branch": [{"_label": "L4"}]}',
                        '{"_label": "L3_sibling"}',
                    ],
                }
            ],
        },
    }
]
```

#### PyArrow Schema (`max_recursion_depth=2`)
```python
pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False),
        pa.field(
            "tree",
            pa.struct(
                [  # depth 1 of TreeType
                    pa.field("_label", pa.string(), nullable=True),
                    pa.field(
                        "branch",
                        pa.list_(
                            pa.struct(
                                [  # depth 2 of TreeType
                                    pa.field("_label", pa.string(), nullable=True),
                                    pa.field(
                                        "branch", pa.list_(pa.string()), nullable=True
                                    ),  # depth > 2 -> JSON string per item
                                ]
                            )
                        ),
                        nullable=True,
                    ),
                ]
            ),
            nullable=True,
        ),
    ]
)
```

---

### Example 5: Model Groups (`xs:group`), Attribute Groups (`xs:attributeGroup`), `<xs:all>`, and Compositors with Duplicate Element Names
Covers:
- `<xs:group ref="..."/>` (including with `maxOccurs="unbounded"`)
- `<xs:attributeGroup ref="..."/>`
- Unrepeated `<xs:choice>` where the **same element name** appears in multiple branches (e.g., `(StartPage, EndPage?, MedlinePgn?) | MedlinePgn` from [`pubmed.xsd` L672–679](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/pubmed/pubmed.xsd#L672-L679)) $\rightarrow$ `MedlinePgn` has `overall_max_occurs == 1`, so it is a **single** `pa.string()` field emitted once.
- `<xs:sequence>` where the **same element name** appears twice (`<xs:element name="interactant"/><xs:element name="interactant"/>`) $\rightarrow$ `overall_max_occurs == 2`, so it is a **list** field emitted once.
- `<xs:all>` compositor.

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:attributeGroup name="auditAttrs">
        <xs:attribute name="createdBy" type="xs:string"/>
        <xs:attribute name="version" type="xs:int"/>
    </xs:attributeGroup>

    <xs:group name="kvGroup">
        <xs:sequence>
            <xs:element name="k" type="xs:string"/>
            <xs:element name="v" type="xs:int"/>
        </xs:sequence>
    </xs:group>

    <xs:complexType name="paginationType">
        <xs:choice>
            <xs:sequence>
                <xs:element name="StartPage" type="xs:string"/>
                <xs:element name="EndPage" type="xs:string" minOccurs="0"/>
                <xs:element name="MedlinePgn" type="xs:string" minOccurs="0"/>
            </xs:sequence>
            <xs:element name="MedlinePgn" type="xs:string"/>
        </xs:choice>
    </xs:complexType>

    <xs:complexType name="metadataAllType">
        <xs:all>
            <xs:element name="author" type="xs:string"/>
            <xs:element name="year" type="xs:int" minOccurs="0"/>
        </xs:all>
    </xs:complexType>

    <xs:element name="dataset">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="item" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="pagination" type="paginationType" minOccurs="0"/>
                            <xs:element name="meta" type="metadataAllType" minOccurs="0"/>
                            <xs:element name="pairs" minOccurs="0">
                                <xs:complexType>
                                    <xs:group ref="kvGroup" maxOccurs="unbounded"/>
                                </xs:complexType>
                            </xs:element>
                        </xs:sequence>
                        <xs:attributeGroup ref="auditAttrs"/>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="item"`)
```xml
<dataset>
    <item createdBy="alice" version="1">
        <pagination>
            <StartPage>10</StartPage>
            <EndPage>20</EndPage>
            <MedlinePgn>10-20</MedlinePgn>
        </pagination>
        <meta>
            <year>2025</year>
            <author>Bob</author>
        </meta>
        <pairs>
            <k>a</k>
            <v>1</v>
            <k>b</k>
            <v>2</v>
        </pairs>
    </item>
    <item>
        <pagination>
            <MedlinePgn>e12345</MedlinePgn>
        </pagination>
    </item>
</dataset>
```

#### Python / JSON Translation
```python
[
    {
        "_createdBy": "alice",
        "_version": 1,
        "pagination": {
            "StartPage": "10",
            "EndPage": "20",
            "MedlinePgn": "10-20",
        },
        "meta": {
            "author": "Bob",
            "year": 2025,
        },
        "pairs": {
            "k": ["a", "b"],
            "v": [1, 2],
        },
    },
    {
        "pagination": {
            "MedlinePgn": "e12345",
        },
    },
]
```

#### PyArrow Schema
```python
pa.schema(
    [
        pa.field("_createdBy", pa.string(), nullable=True),
        pa.field("_version", pa.int32(), nullable=True),
        pa.field(
            "pagination",
            pa.struct(
                [
                    pa.field("StartPage", pa.string(), nullable=True),
                    pa.field("EndPage", pa.string(), nullable=True),
                    pa.field("MedlinePgn", pa.string(), nullable=True),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "meta",
            pa.struct(
                [
                    pa.field("author", pa.string(), nullable=False),
                    pa.field("year", pa.int32(), nullable=True),
                ]
            ),
            nullable=True,
        ),
        pa.field(
            "pairs",
            pa.struct(
                [
                    pa.field("k", pa.list_(pa.string()), nullable=False),
                    pa.field("v", pa.list_(pa.int32()), nullable=False),
                ]
            ),
            nullable=True,
        ),
    ]
)
```

---

### Example 6: Context-Dependent Same `(parent, child)` Names with Different Types
In [`06_context_dependent_paths/schema.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/06_context_dependent_paths/schema.xsd) and [`10_combinations/schema.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/10_combinations/schema.xsd#L195-L226), `left/node/leaf` is a single `xs:string` while `right/node/leaf` is an unbounded `LabelType` (complex `simpleContent` with attribute `lang`).
Because the full ancestor path from `xml_tag` (`("root", "left", "node")` vs `("root", "right", "node")`) distinguishes them, each gets its exact PyArrow type instead of merging `("node", "leaf")` globally.

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:element name="root">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="left">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="node">
                                <xs:complexType>
                                    <xs:sequence>
                                        <xs:element name="leaf" type="xs:string"/>
                                    </xs:sequence>
                                </xs:complexType>
                            </xs:element>
                        </xs:sequence>
                    </xs:complexType>
                </xs:element>
                <xs:element name="right">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="node">
                                <xs:complexType>
                                    <xs:sequence>
                                        <xs:element name="leaf" maxOccurs="unbounded">
                                            <xs:complexType>
                                                <xs:simpleContent>
                                                    <xs:extension base="xs:string">
                                                        <xs:attribute name="lang" type="xs:string"/>
                                                    </xs:extension>
                                                </xs:simpleContent>
                                            </xs:complexType>
                                        </xs:element>
                                    </xs:sequence>
                                </xs:complexType>
                            </xs:element>
                        </xs:sequence>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="root"`)
```xml
<root>
    <left>
        <node>
            <leaf>l1</leaf>
        </node>
    </left>
    <right>
        <node>
            <leaf>r1</leaf>
            <leaf lang="en">r2</leaf>
        </node>
    </right>
</root>
```

#### Python / JSON Translation
```python
[
    {
        "left": {
            "node": {
                "leaf": "l1",
            }
        },
        "right": {
            "node": {
                "leaf": [
                    {"#text": "r1"},
                    {"_lang": "en", "#text": "r2"},
                ]
            }
        },
    }
]
```

#### PyArrow Schema
```python
pa.schema(
    [
        pa.field(
            "left",
            pa.struct(
                [
                    pa.field(
                        "node",
                        pa.struct(
                            [
                                pa.field("leaf", pa.string(), nullable=False),
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
                                        pa.struct(
                                            [
                                                pa.field("_lang", pa.string(), nullable=True),
                                                pa.field("#text", pa.string(), nullable=True),
                                            ]
                                        )
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
```

---

### Example 7: Multiple Root Elements Sharing the Same Child `xml_tag` (`uniref.xsd` Pattern)
In [`uniref.xsd`](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/uniprot/uniref/uniref.xsd), four root elements (`UniRef100`, `UniRef90`, `UniRef50`, `UniRef`) each declare a child `<xs:element name="entry" type="entryType"/>`. Because all four `entry` declarations share the same `entryType`, resolving `xml_tag="entry"` succeeds unambiguously.

---

### Example 8: Field Nullability & Empty String Tag Coercion
Covers:
1. Required attributes (`use="required"` $\rightarrow$ `nullable=False`) vs optional attributes (`use="optional"` $\rightarrow$ `nullable=True`).
2. Required non-nillable elements (`minOccurs="1"`, `nillable="false"` $\rightarrow$ `nullable=False`), with present empty string tags `<summary/>` coerced to `""`.
3. Elements inside `<xs:choice>` branches (`nullable=True`).
4. Children of `nillable="true"` complex elements (`nullable=True` so `<comment type="function" xsi:nil="true"/>` can omit child `text`).

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:complexType name="commentType">
        <xs:sequence>
            <xs:element name="text" type="xs:string"/>
        </xs:sequence>
        <xs:attribute name="type" type="xs:string" use="required"/>
    </xs:complexType>

    <xs:element name="root">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="entry" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="name" type="xs:string"/>
                            <xs:element name="summary" type="xs:string"/>
                            <xs:element name="optionalNote" type="xs:string" minOccurs="0"/>
                            <xs:choice>
                                <xs:element name="geneId" type="xs:int"/>
                                <xs:element name="proteinId" type="xs:string"/>
                            </xs:choice>
                            <xs:element name="comment" type="commentType" nillable="true"/>
                        </xs:sequence>
                        <xs:attribute name="id" type="xs:string" use="required"/>
                        <xs:attribute name="source" type="xs:string" use="optional"/>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="entry"`)
```xml
<root xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
    <entry id="e1">
        <name>BRCA1</name>
        <summary/>
        <geneId>672</geneId>
        <comment type="function" xsi:nil="true"/>
    </entry>
</root>
```

#### Python / JSON Translation
```python
[
    {
        "_id": "e1",
        "name": "BRCA1",
        "summary": "",  # coerced from None -> "" to satisfy nullable=False
        "geneId": 672,
        "comment": {"_type": "function"},
    }
]
```

#### PyArrow Schema
```python
pa.schema(
    [
        pa.field("_id", pa.string(), nullable=False),  # use="required"
        pa.field("_source", pa.string(), nullable=True),  # use="optional"
        pa.field("name", pa.string(), nullable=False),  # minOccurs=1, non-nillable
        pa.field("summary", pa.string(), nullable=False),  # minOccurs=1, non-nillable
        pa.field("optionalNote", pa.string(), nullable=True),  # minOccurs=0
        pa.field("geneId", pa.int32(), nullable=True),  # inside xs:choice -> nullable=True
        pa.field("proteinId", pa.string(), nullable=True),  # inside xs:choice -> nullable=True
        pa.field(
            "comment",
            pa.struct(
                [
                    pa.field("_type", pa.string(), nullable=False),  # use="required"
                    pa.field("text", pa.string(), nullable=True),  # parent is nillable="true" -> nullable=True
                ]
            ),
            nullable=True,
        ),  # nillable="true" -> nullable=True
    ]
)
```

---

### Example 9: Purely Empty `complexType` Marker Elements
Covers DTD-derived empty marker elements like [`pubmed.xsd` L1312](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/pubmed/pubmed.xsd#L1312) (`<BookInProcess/>`) and [`10_combinations/schema.xsd` L11–12](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/10_combinations/schema.xsd#L11-L12) (`EmptyType`):
- Mapped to `pa.string()` (or `pa.list_(pa.string())` when `overall_max_occurs > 1`).
- Coerced to `""` when present in the XML, and `None` when omitted.

#### XSD
```xml
<?xml version="1.0" encoding="UTF-8"?>
<xs:schema xmlns:xs="http://www.w3.org/2001/XMLSchema">
    <xs:complexType name="EmptyType"/>

    <xs:element name="books">
        <xs:complexType>
            <xs:sequence>
                <xs:element name="book" maxOccurs="unbounded">
                    <xs:complexType>
                        <xs:sequence>
                            <xs:element name="title" type="xs:string"/>
                            <xs:element name="BookInProcess" type="EmptyType" minOccurs="0"/>
                            <xs:element name="flag" type="EmptyType" minOccurs="0" maxOccurs="unbounded"/>
                        </xs:sequence>
                    </xs:complexType>
                </xs:element>
            </xs:sequence>
        </xs:complexType>
    </xs:element>
</xs:schema>
```

#### Sample XML (`xml_tag="book"`)
```xml
<books>
    <book>
        <title>In-Process Book</title>
        <BookInProcess/>
        <flag/>
        <flag/>
    </book>
    <book>
        <title>Completed Book</title>
    </book>
</books>
```

#### Python / JSON Translation
```python
[
    {
        "title": "In-Process Book",
        "BookInProcess": "",
        "flag": ["", ""],
    },
    {
        "title": "Completed Book",
    },
]
```

#### PyArrow Schema
```python
pa.schema(
    [
        pa.field("title", pa.string(), nullable=False),
        pa.field("BookInProcess", pa.string(), nullable=True),
        pa.field("flag", pa.list_(pa.string()), nullable=True),
    ]
)
```

---

## 3. Error Specification (Unsupported / Invalid XSD Features)

The converter raises a descriptive `ValueError` during schema generation when encountering any of the following:

| # | Error Case | Example XSD Pattern | Reason |
| :--- | :--- | :--- | :--- |
| **E1** | **`xml_tag` Not Found** | `<xs:schema/>` or `xml_tag="nonexistent"` | Target element does not exist in the XSD. |
| **E2** | **Ambiguous `xml_tag` Type** | Root `A` has `<entry type="Type1"/>` and root `B` has `<entry type="Type2"/>` | Cannot determine a single `pa.Schema` for `xml_tag="entry"` when declarations at the same resolution depth have different types. |
| **E3** | **Non-Struct Top-Level `xml_tag`** | `<xs:element name="leaf" type="xs:string"/>` or empty `<xs:complexType/>` with no attributes/children | A `pa.Table` schema requires top-level columns (at least one attribute, child element, or `#text` field). |
| **E4** | **Element Wildcard (`<xs:any>`)** | `<xs:sequence><xs:any processContents="lax"/></xs:sequence>` | Arbitrary undeclared child elements cannot be represented as fixed fields in a `pa.struct`. |
| **E5** | **Attribute Wildcard (`<xs:anyAttribute>`)** | `<xs:anyAttribute processContents="skip"/>` | Arbitrary undeclared attributes cannot be represented as fixed fields in a `pa.struct`. |
| **E6** | **Abstract Element or Abstract Type (`abstract="true"`)** | `<xs:complexType name="Base" abstract="true"/>` used as an element's type, or `<xs:element name="item" abstract="true"/>` | Concrete runtime instances require `xsi:type` or `substitutionGroup` polymorphism with varying schemas. |
| **E7** | **Substitution Group with Heterogeneous Types** | `<xs:element name="sub" substitutionGroup="head" type="OtherType"/>` | Substituted elements introduce different element names/types not declared in the parent's content model. |
| **E8** | **Prefixed Attribute / Child Name Collision** | Element has both `<xs:attribute name="id"/>` (mapped to `_id`) and `<xs:element name="_id"/>` | Both map to field name `"_id"` in the same `pa.struct`. |
| **E9** | **Non-JSON-Serializable Scalar Type Under Recursion** | Recursive element subtree contains `xs:date`, `xs:dateTime`, `xs:hexBinary`, or `xs:base64Binary` | Elements truncated beyond `max_recursion_depth` are serialized to JSON (`json.dumps`), which does not natively serialize `datetime.date` or `bytes` without lossy/ambiguous string encoding. |

---

## 4. Generated Specification Fixtures, Real-World Examples & Tests

### 4.1 Specification Fixtures (`tests/data/xsd/spec/`)

| Suite | `xml_tag` | XSD | Expected PyArrow Schema | Sample XML | Expected JSON |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`01_scalars_and_simple_types`** | `record` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/01_scalars_and_simple_types/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/01_scalars_and_simple_types/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/01_scalars_and_simple_types/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/01_scalars_and_simple_types/expected.json) |
| **`02_namespaces_and_nillable`** | `entry` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/02_namespaces_and_nillable/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/02_namespaces_and_nillable/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/02_namespaces_and_nillable/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/02_namespaces_and_nillable/expected.json) |
| **`03_mixed_content_raw_text`** | `article` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/03_mixed_content_raw_text/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/03_mixed_content_raw_text/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/03_mixed_content_raw_text/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/03_mixed_content_raw_text/expected.json) |
| **`04_recursion_depth`** | `record` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/04_recursion_depth/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/04_recursion_depth/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/04_recursion_depth/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/04_recursion_depth/expected.json), [expected_depth_1.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/04_recursion_depth/expected_depth_1.json) |
| **`05_groups_and_compositors`** | `item` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/05_groups_and_compositors/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/05_groups_and_compositors/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/05_groups_and_compositors/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/05_groups_and_compositors/expected.json) |
| **`06_context_dependent_paths`** | `root` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/06_context_dependent_paths/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/06_context_dependent_paths/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/06_context_dependent_paths/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/06_context_dependent_paths/expected.json) |
| **`07_multi_root_shared_type`** | `entry` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/07_multi_root_shared_type/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/07_multi_root_shared_type/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/07_multi_root_shared_type/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/07_multi_root_shared_type/expected.json) |
| **`08_nullability_and_empty_tags`** | `entry` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/08_nullability_and_empty_tags/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/08_nullability_and_empty_tags/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/08_nullability_and_empty_tags/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/08_nullability_and_empty_tags/expected.json) |
| **`09_empty_marker_elements`** | `book` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/09_empty_marker_elements/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/09_empty_marker_elements/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/09_empty_marker_elements/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/09_empty_marker_elements/expected.json) |
| **`10_combinations`** | `record` | [schema.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/10_combinations/schema.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/10_combinations/expected_schema.py) | [data.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/10_combinations/data.xml) | [expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/10_combinations/expected.json) |

### 4.2 Error Specification Fixtures (`tests/data/xsd/spec/errors/`)

- [ambiguous_xml_tag.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/ambiguous_xml_tag.xsd) (E2)
- [scalar_top_level.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/scalar_top_level.xsd) (E3)
- [empty_marker_top_level.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/empty_marker_top_level.xsd) (E3)
- [wildcard_any_element.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/wildcard_any_element.xsd) (E4)
- [wildcard_any_attribute.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/wildcard_any_attribute.xsd) (E5)
- [abstract_element.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/abstract_element.xsd) (E6)
- [abstract_type.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/abstract_type.xsd) (E6)
- [substitution_group.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/substitution_group.xsd) (E7)
- [attr_child_name_collision.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/attr_child_name_collision.xsd) (E8)
- [nested_date_under_recursion.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/nested_date_under_recursion.xsd) (E9)
- [nested_binary_under_recursion.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/spec/errors/nested_binary_under_recursion.xsd) (E9)

### 4.3 Real-World Minimal & Maximal Examples (`tests/data/xsd/examples/`)

| Dataset | `xml_tag` | XSD | Expected PyArrow Schema | Minimal XML & Expected JSON | Maximal XML & Expected JSON |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **UniRef** | `entry` | [uniref.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/uniprot/uniref/uniref.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniref/expected_schema.py) | [minimal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniref/minimal.xml), [minimal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniref/minimal_expected.json) | [maximal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniref/maximal.xml), [maximal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniref/maximal_expected.json) |
| **UniProtKB** | `entry` | [uniprot.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/uniprot/uniprot_kb/uniprot.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniprot/expected_schema.py) | [minimal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniprot/minimal.xml), [minimal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniprot/minimal_expected.json) | [maximal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniprot/maximal.xml), [maximal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/uniprot/maximal_expected.json) |
| **PubMed (`PubmedArticle`)** | `PubmedArticle` | [pubmed.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/pubmed/pubmed.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_article/expected_schema.py) | [minimal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_article/minimal.xml), [minimal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_article/minimal_expected.json) | [maximal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_article/maximal.xml), [maximal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_article/maximal_expected.json) |
| **PubMed (`PubmedBookArticle`)** | `PubmedBookArticle` | [pubmed.xsd](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/pubmed/pubmed.xsd) | [expected_schema.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_book_article/expected_schema.py) | [minimal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_book_article/minimal.xml), [minimal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_book_article/minimal_expected.json) | [maximal.xml](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_book_article/maximal.xml), [maximal_expected.json](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/data/xsd/examples/pubmed_book_article/maximal_expected.json) |

### 4.4 Specification Test Suite

- [test_xsd_spec.py](file:///Users/gwg/code/kbase/cdm-data-loaders-no-batcher/tests/cdm_data_loaders/readers/test_xsd_spec.py): Validates all fixtures against their XSDs, tests `generate_pyarrow_schema`, tests error cases E1–E9, tests `convert_xml_element` + PyArrow Table / Parquet round-trips, and tests `process_xml_file_to_dict` yielding `pa.Table` batches.

