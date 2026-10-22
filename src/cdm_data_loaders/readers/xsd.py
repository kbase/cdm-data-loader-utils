"""Functions for reading XSD files, generating PyArrow schemas, and converting XML elements."""

import base64
import datetime as dt
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final, Literal

import pyarrow as pa
import xmlschema
from lxml.etree import Element
from xmlschema.validators import (
    XsdAnyAttribute,
    XsdAnyElement,
    XsdAttribute,
    XsdComplexType,
    XsdElement,
    XsdGroup,
    XsdSimpleType,
    XsdType,
)

from cdm_data_loaders.utils.xml_utils import local_name as _local_name

DEFAULT_MAX_RECURSION_DEPTH: Final[int] = 2

XSD_NAMESPACE: Final[str] = "{http://www.w3.org/2001/XMLSchema}"
XSI_NAMESPACE: Final[str] = "{http://www.w3.org/2001/XMLSchema-instance}"
XSI_NIL: Final[str] = f"{XSI_NAMESPACE}nil"

FORBIDDEN_RECURSIVE_XSD_TYPES: Final[frozenset[str]] = frozenset(
    {
        "xs:date",
        "xs:dateTime",
        "xs:hexBinary",
        "xs:base64Binary",
    }
)

BUILTIN_LIST_ITEM_TYPES: Final[Mapping[str, str]] = {
    "xs:NMTOKENS": "xs:NMTOKEN",
    "xs:IDREFS": "xs:IDREF",
    "xs:ENTITIES": "xs:ENTITY",
}

SCALAR_ARROW_TYPES: Final[Mapping[str, pa.DataType]] = {
    "xs:string": pa.string(),
    "xs:normalizedString": pa.string(),
    "xs:token": pa.string(),
    "xs:language": pa.string(),
    "xs:Name": pa.string(),
    "xs:NCName": pa.string(),
    "xs:NMTOKEN": pa.string(),
    "xs:ID": pa.string(),
    "xs:IDREF": pa.string(),
    "xs:ENTITY": pa.string(),
    "xs:anyURI": pa.string(),
    "xs:QName": pa.string(),
    "xs:NOTATION": pa.string(),
    "xs:dateTime": pa.string(),
    "xs:time": pa.string(),
    "xs:duration": pa.string(),
    "xs:gYear": pa.string(),
    "xs:gYearMonth": pa.string(),
    "xs:gMonthDay": pa.string(),
    "xs:gDay": pa.string(),
    "xs:gMonth": pa.string(),
    "xs:boolean": pa.bool_(),
    "xs:byte": pa.int8(),
    "xs:short": pa.int16(),
    "xs:int": pa.int32(),
    "xs:long": pa.int64(),
    "xs:integer": pa.int64(),
    "xs:positiveInteger": pa.int64(),
    "xs:nonNegativeInteger": pa.int64(),
    "xs:negativeInteger": pa.int64(),
    "xs:nonPositiveInteger": pa.int64(),
    "xs:unsignedByte": pa.uint8(),
    "xs:unsignedShort": pa.uint16(),
    "xs:unsignedInt": pa.uint32(),
    "xs:unsignedLong": pa.uint64(),
    "xs:float": pa.float64(),
    "xs:double": pa.float64(),
    "xs:decimal": pa.float64(),
    "xs:date": pa.date32(),
    "xs:base64Binary": pa.binary(),
    "xs:hexBinary": pa.binary(),
}


def load_schema(schema_path: str | Path) -> xmlschema.XMLSchema:
    """Load an XML schema document into an XMLSchema object.

    :param schema_path: path to the schema file
    :type  schema_path: str | Path
    :return: loaded schema
    :rtype : xmlschema.XMLSchema
    """
    return xmlschema.XMLSchema(str(schema_path))


def _extract_doc_comment(component: object, fallback_type: object = None) -> str | None:
    """Extract whitespace-normalized ``<xs:documentation>`` text from an XSD component or its type."""
    candidates = [component]
    ref = getattr(component, "ref", None)
    if ref is not None and ref is not component:
        candidates.append(ref)
    if fallback_type is not None:
        candidates.append(fallback_type)

    for candidate in candidates:
        annotation = getattr(candidate, "annotation", None)
        docs = getattr(annotation, "documentation", None)
        if not docs:
            continue
        parts = [" ".join(doc.text.split()) for doc in docs if getattr(doc, "text", None) and doc.text.strip()]
        if parts:
            return " ".join(parts)
    return None


def _build_field_metadata(
    *,
    xsd_type: str | None = None,
    comment: str | None = None,
) -> dict[str, str] | None:
    """Construct the ``pa.Field`` metadata dictionary from optional comment and XSD type strings."""
    metadata: dict[str, str] = {}
    if comment:
        metadata["comment"] = comment
    if xsd_type:
        metadata["xsd_type"] = xsd_type
    return metadata or None


def _builtin_xsd_name(simple_type: XsdType) -> str:
    """Resolve an atomic simpleType or restriction to its built-in ``xs:<name>`` identifier."""
    cur: Any = simple_type
    while cur is not None:
        name = getattr(cur, "name", None)
        if isinstance(name, str) and name.startswith(XSD_NAMESPACE):
            return f"xs:{cur.local_name}"
        cur = getattr(cur, "base_type", None)
    return "xs:string"


def _unwrap_member_types(simple_type: XsdType) -> list[XsdSimpleType]:
    """Return the member types of an ``xs:union`` simpleType (unwrapping restrictions)."""
    cur: Any = simple_type
    while cur is not None:
        members = getattr(cur, "member_types", None)
        if members:
            return list(members)
        cur = getattr(cur, "base_type", None)
    return []


def _unwrap_item_type(simple_type: XsdType) -> XsdSimpleType | None:
    """Return the item type of an ``xs:list`` simpleType (unwrapping restrictions)."""
    cur: Any = simple_type
    while cur is not None:
        item_type = getattr(cur, "item_type", None)
        if item_type is not None:
            return item_type
        cur = getattr(cur, "base_type", None)
    return None


def _resolve_scalar_type(simple_type: XsdType) -> tuple[pa.DataType, str]:
    """Map a non-list XSD simpleType (atomic or union) to ``(pa.DataType, xsd_type_str)``."""
    if simple_type.is_union():
        members = _unwrap_member_types(simple_type)
        member_pairs = [_resolve_scalar_type(m) for m in members]
        xsd_type_str = " | ".join(xsd_name for _, xsd_name in member_pairs)
        unique_arrow_types = {arrow_t for arrow_t, _ in member_pairs}
        arrow_type = next(iter(unique_arrow_types)) if len(unique_arrow_types) == 1 else pa.string()
        return arrow_type, xsd_type_str

    xsd_name = _builtin_xsd_name(simple_type)
    arrow_type = SCALAR_ARROW_TYPES.get(xsd_name, pa.string())
    return arrow_type, xsd_name


def _simple_type_to_arrow(simple_type: XsdType) -> tuple[pa.DataType, str | None]:
    """Map any XSD simpleType (atomic, union, or list) to ``(pa.DataType, leaf_xsd_type_or_None)``.

    When ``simple_type`` is a list type (``xs:list``, ``xs:NMTOKENS``, ``xs:IDREFS``, ``xs:ENTITIES``),
    the returned ``pa.DataType`` is ``pa.list_(pa.field("item", item_arrow_type, nullable=True, metadata={"xsd_type": ...}))``
    and the second element is ``None`` because ``"xsd_type"`` is attached to the inner ``"item"`` field.
    """
    builtin_name = _builtin_xsd_name(simple_type)
    if builtin_name in BUILTIN_LIST_ITEM_TYPES:
        item_xsd = BUILTIN_LIST_ITEM_TYPES[builtin_name]
        item_field = pa.field("item", pa.string(), nullable=True, metadata={"xsd_type": item_xsd})
        return pa.list_(item_field), None

    if simple_type.is_list():
        item_type = _unwrap_item_type(simple_type)
        if item_type is not None:
            item_arrow, item_xsd = _resolve_scalar_type(item_type)
        else:
            item_arrow, item_xsd = pa.string(), "xs:string"
        item_field = pa.field("item", item_arrow, nullable=True, metadata={"xsd_type": item_xsd})
        return pa.list_(item_field), None

    return _resolve_scalar_type(simple_type)


def _iter_simple_builtin_names(simple_type: XsdType) -> list[str]:
    """Return all constituent built-in ``xs:*`` type names for an atomic, list, or union simpleType."""
    builtin_name = _builtin_xsd_name(simple_type)
    if builtin_name in BUILTIN_LIST_ITEM_TYPES:
        return [BUILTIN_LIST_ITEM_TYPES[builtin_name]]
    if simple_type.is_list():
        item_type = _unwrap_item_type(simple_type)
        return _iter_simple_builtin_names(item_type) if item_type is not None else ["xs:string"]
    if simple_type.is_union():
        names: list[str] = []
        for member in _unwrap_member_types(simple_type):
            names.extend(_iter_simple_builtin_names(member))
        return names
    return [builtin_name]


@dataclass(frozen=True, slots=True)
class _TraversalContext:
    """Immutable traversal context for XSD schema generation and XML-to-dict conversion."""

    schema: xmlschema.XMLSchema
    type_counts: Mapping[int, int]
    max_recursion_depth: int | None


_LeafKind = Literal["simple", "mixed", "empty_marker"]


def _check_complex_type_forbidden_recursion(root_type: XsdType, current: XsdComplexType) -> list[XsdType]:
    """Validate attributes on ``current`` and return child ``XsdType``s to inspect for forbidden recursive scalars."""
    for attr in current.attributes.values():
        if isinstance(attr, XsdAttribute) and attr.type is not None:
            for xsd_name in _iter_simple_builtin_names(attr.type):
                if xsd_name in FORBIDDEN_RECURSIVE_XSD_TYPES:
                    msg = (
                        f"Unsupported XSD feature: recursive type {root_type.local_name or root_type.name!r} "
                        f"contains attribute @{attr.local_name} of type {xsd_name!r}"
                    )
                    raise ValueError(msg)

    if current.has_simple_content() and isinstance(current.content, XsdType):
        return [current.content]
    if not current.mixed and current.model_group is not None:
        return [
            child.type
            for child in current.model_group.iter_elements()
            if isinstance(child, XsdElement) and child.type is not None and child.max_occurs != 0
        ]
    return []


def _check_no_forbidden_types_under_recursion(root_type: XsdType) -> None:
    """Raise ValueError if a recursive type or any descendant type uses date, dateTime, or binary scalars."""
    visited: set[int] = set()
    stack: list[XsdType] = [root_type]

    while stack:
        current = stack.pop()
        if id(current) in visited:
            continue
        visited.add(id(current))

        if not isinstance(current, XsdComplexType):
            for xsd_name in _iter_simple_builtin_names(current):
                if xsd_name in FORBIDDEN_RECURSIVE_XSD_TYPES:
                    msg = (
                        f"Unsupported XSD feature: recursive type {root_type.local_name or root_type.name!r} "
                        f"contains {xsd_name!r}, which cannot be serialized to JSON when recursion is truncated"
                    )
                    raise ValueError(msg)
            continue

        stack.extend(_check_complex_type_forbidden_recursion(root_type, current))


def _group_min_occurs_by_name(group: XsdGroup, local_name: str) -> int:
    """Compute minimum occurrences of elements with ``local_name`` across a model group."""
    if group.min_occurs == 0 or group.max_occurs == 0:
        return 0
    target_group = group.ref if group.ref is not None else group
    counts: list[int] = []
    for item in target_group.content:
        if isinstance(item, XsdGroup):
            counts.append(_group_min_occurs_by_name(item, local_name))
        elif isinstance(item, XsdElement) and item.max_occurs != 0 and item.local_name == local_name:
            counts.append(item.min_occurs)
        else:
            counts.append(0)
    if not counts:
        return 0
    branch_min = min(counts) if target_group.model == "choice" else sum(counts)
    return group.min_occurs * branch_min


def _particle_max_occurs_by_name(item: object, local_name: str) -> int | None:
    """Compute maximum occurrences (None for unbounded) of ``local_name`` for a single group particle."""
    if isinstance(item, XsdGroup):
        return _group_max_occurs_by_name(item, local_name)
    if isinstance(item, XsdElement) and item.max_occurs != 0 and item.local_name == local_name:
        return item.max_occurs
    return 0


def _group_max_occurs_by_name(group: XsdGroup, local_name: str) -> int | None:
    """Compute maximum occurrences (None for unbounded) of elements with ``local_name`` across a model group."""
    if group.max_occurs == 0:
        return 0
    target_group = group.ref if group.ref is not None else group
    branch_maxes = [_particle_max_occurs_by_name(item, local_name) for item in target_group.content]
    if not branch_maxes:
        return 0
    if any(m is None for m in branch_maxes):
        inner_max: int | None = None
    elif target_group.model == "choice":
        inner_max = max(m for m in branch_maxes if m is not None)
    else:
        inner_max = sum(m for m in branch_maxes if m is not None)
    if inner_max == 0:
        return 0
    if inner_max is None or group.max_occurs is None:
        return None
    return group.max_occurs * inner_max


def _is_any_element(item: object) -> bool:
    """Return True if ``item`` is an ``xs:any`` wildcard particle."""
    return isinstance(item, XsdAnyElement)


def _is_any_attribute(attr: object) -> bool:
    """Return True if ``attr`` is an ``xs:anyAttribute`` wildcard declaration."""
    return isinstance(attr, XsdAnyAttribute)


def _validate_group_particles(schema: xmlschema.XMLSchema, group: XsdGroup) -> None:
    """Validate that a model group contains no wildcards, abstract elements/types, or substitution groups."""
    target_group = group.ref if group.ref is not None else group
    for item in target_group.content:
        if _is_any_element(item):
            msg = "Unsupported XSD feature: xs:any wildcard element is not supported"
            raise ValueError(msg)
        if isinstance(item, XsdGroup):
            _validate_group_particles(schema, item)
        elif isinstance(item, XsdElement):
            target_el = item.ref if item.ref is not None else item
            if getattr(item, "abstract", False) or getattr(target_el, "abstract", False):
                msg = f"Unsupported XSD feature: abstract element {item.local_name!r} is not supported"
                raise ValueError(msg)
            if (
                item.name in schema.substitution_groups
                or target_el.name in schema.substitution_groups
                or getattr(item, "substitution_group", None) is not None
                or getattr(target_el, "substitution_group", None) is not None
            ):
                msg = f"Unsupported XSD feature: substitutionGroup on element {item.local_name!r} is not supported"
                raise ValueError(msg)
            if item.type is not None and getattr(item.type, "abstract", False):
                msg = f"Unsupported XSD feature: abstract type {item.type.local_name!r} is not supported"
                raise ValueError(msg)


def _collect_unique_children(schema: xmlschema.XMLSchema, complex_type: XsdComplexType) -> list[XsdElement]:
    """Validate and return unique child ``XsdElement`` declarations in content-model order."""
    model_group = complex_type.model_group
    if model_group is None:
        return []
    _validate_group_particles(schema, model_group)
    seen_names: set[str] = set()
    children: list[XsdElement] = []
    for child in model_group.iter_elements():
        if not isinstance(child, XsdElement):
            continue
        max_occurs = _group_max_occurs_by_name(model_group, child.local_name)
        if max_occurs == 0 or child.local_name in seen_names:
            continue
        seen_names.add(child.local_name)
        children.append(child)
    return children


def _collect_attributes(complex_type: XsdComplexType) -> list[XsdAttribute]:
    """Validate and return declared attributes on ``complex_type`` in declaration order."""
    attrs: list[XsdAttribute] = []
    for attr in complex_type.attributes.values():
        if _is_any_attribute(attr):
            msg = "Unsupported XSD feature: xs:anyAttribute wildcard is not supported"
            raise ValueError(msg)
        if isinstance(attr, XsdAttribute) and getattr(attr, "use", None) != "prohibited":
            if attr.local_name is None:
                continue
            attrs.append(attr)
    return attrs


def _is_empty_marker_complex_type(complex_type: XsdComplexType) -> bool:
    """Return True if ``complex_type`` has no attributes, no child elements, not simpleContent, and not mixed."""
    if complex_type.has_simple_content() or complex_type.mixed:
        return False
    if any(
        isinstance(a, XsdAttribute) and getattr(a, "use", None) != "prohibited"
        for a in complex_type.attributes.values()
    ):
        return False
    if complex_type.model_group is None:
        return True
    return not any(
        isinstance(el, XsdElement) and _group_max_occurs_by_name(complex_type.model_group, el.local_name) != 0
        for el in complex_type.model_group.iter_elements()
    )


def _resolve_child_cardinality(
    model_group: XsdGroup,
    child_decl: XsdElement,
    *,
    parent_nillable: bool = False,
) -> tuple[bool, bool]:
    """Return ``(is_list, nullable)`` for ``child_decl`` within ``model_group``."""
    child_name = child_decl.local_name or ""
    max_occurs = _group_max_occurs_by_name(model_group, child_name)
    min_occurs = _group_min_occurs_by_name(model_group, child_name)
    is_list = max_occurs is None or max_occurs > 1
    nullable = parent_nillable or child_decl.nillable or min_occurs < 1
    return is_list, nullable


def _build_complex_struct_fields(
    ctx: _TraversalContext,
    complex_type: XsdComplexType,
    *,
    parent_nillable: bool,
) -> list[pa.Field]:
    """Build the ordered list of ``pa.Field``s for a struct-mapped ``XsdComplexType``."""
    fields: list[pa.Field] = []
    attr_field_names: set[str] = set()

    for attr in _collect_attributes(complex_type):
        field_name = f"_{attr.local_name}"
        attr_field_names.add(field_name)
        attr_type = attr.type
        if attr_type is not None:
            pa_type, xsd_type_str = _simple_type_to_arrow(attr_type)
        else:
            pa_type, xsd_type_str = pa.string(), "xs:string"
        comment = _extract_doc_comment(attr, attr_type)
        nullable = parent_nillable or getattr(attr, "use", "optional") != "required"
        fields.append(
            pa.field(
                field_name,
                pa_type,
                nullable=nullable,
                metadata=_build_field_metadata(xsd_type=xsd_type_str, comment=comment),
            )
        )

    if complex_type.has_simple_content():
        content_type = complex_type.content
        if isinstance(content_type, XsdType):
            pa_type, xsd_type_str = _simple_type_to_arrow(content_type)
        else:
            pa_type, xsd_type_str = pa.string(), "xs:string"
        fields.append(
            pa.field(
                "#text",
                pa_type,
                nullable=True,
                metadata=_build_field_metadata(xsd_type=xsd_type_str),
            )
        )
        return fields

    if complex_type.mixed:
        fields.append(
            pa.field(
                "#text",
                pa.string(),
                nullable=True,
                metadata=_build_field_metadata(xsd_type="xs:string"),
            )
        )
        return fields

    model_group = complex_type.model_group
    if model_group is None:
        return fields

    for child in _collect_unique_children(ctx.schema, complex_type):
        child_name = child.local_name
        if child_name in attr_field_names:
            msg = (
                f"Field name collision in {complex_type.local_name or 'complexType'!r}: "
                f"attribute and child element both map to {child_name!r}"
            )
            raise ValueError(msg)

        is_list, nullable = _resolve_child_cardinality(model_group, child, parent_nillable=parent_nillable)

        fields.append(
            _build_element_field(
                ctx,
                child,
                is_list=is_list,
                nullable=nullable,
            )
        )

    return fields


def _classify_leaf_element_type(elem_type: XsdType | None) -> tuple[_LeafKind, XsdType | None] | None:
    """Classify ``elem_type`` as a non-struct leaf ``(kind, simple_type)`` or ``None`` if it maps to a struct."""
    if not isinstance(elem_type, XsdComplexType):
        return "simple", elem_type

    attrs = _collect_attributes(elem_type)
    if elem_type.has_simple_content() and not attrs:
        content_type = elem_type.content if isinstance(elem_type.content, XsdType) else None
        return "simple", content_type

    if elem_type.mixed and not attrs:
        return "mixed", None

    if _is_empty_marker_complex_type(elem_type):
        return "empty_marker", None

    return None


def _resolve_non_struct_element_type(elem_type: XsdType | None) -> tuple[pa.DataType, str | None] | None:
    """Return ``(pa_type, xsd_type_or_None)`` if ``elem_type`` maps to a leaf scalar rather than a struct."""
    leaf_info = _classify_leaf_element_type(elem_type)
    if leaf_info is None:
        return None
    kind, simple_type = leaf_info
    if kind == "simple":
        return _simple_type_to_arrow(simple_type) if simple_type is not None else (pa.string(), "xs:string")
    if kind == "mixed":
        return pa.string(), "xs:string"
    return pa.string(), None


def _build_element_field(
    ctx: _TraversalContext,
    element: XsdElement,
    *,
    is_list: bool,
    nullable: bool,
) -> pa.Field:
    """Build a ``pa.Field`` for a child ``XsdElement`` given its cardinality and nullability."""
    elem_type = element.type
    comment = _extract_doc_comment(element, elem_type)

    leaf_type_info = _resolve_non_struct_element_type(elem_type)
    if leaf_type_info is not None:
        return _wrap_scalar_or_list_field(
            element.local_name,
            leaf_type_info,
            comment=comment,
            is_list=is_list,
            nullable=nullable,
        )

    assert isinstance(elem_type, XsdComplexType)  # noqa: S101
    type_id = id(elem_type)
    if type_id in ctx.type_counts:
        _check_no_forbidden_types_under_recursion(elem_type)

    next_count = ctx.type_counts.get(type_id, 0) + 1
    if ctx.max_recursion_depth is not None and next_count > ctx.max_recursion_depth:
        return _wrap_scalar_or_list_field(
            element.local_name,
            (pa.string(), None),
            comment=comment,
            is_list=is_list,
            nullable=nullable,
        )

    next_ctx = _TraversalContext(
        schema=ctx.schema,
        type_counts={**ctx.type_counts, type_id: next_count},
        max_recursion_depth=ctx.max_recursion_depth,
    )
    struct_fields = _build_complex_struct_fields(
        next_ctx,
        elem_type,
        parent_nillable=element.nillable,
    )
    struct_type = pa.struct(struct_fields)
    field_type = pa.list_(pa.field("item", struct_type, nullable=True)) if is_list else struct_type
    return pa.field(
        element.local_name,
        field_type,
        nullable=nullable,
        metadata=_build_field_metadata(comment=comment),
    )


def _wrap_scalar_or_list_field(
    name: str,
    type_info: tuple[pa.DataType, str | None],
    *,
    comment: str | None,
    is_list: bool,
    nullable: bool,
) -> pa.Field:
    """Create a scalar or list ``pa.Field`` with ``xsd_type`` on leaf scalar fields and ``comment`` on the outer field."""
    pa_type, xsd_type = type_info
    if is_list:
        item_field = pa.field("item", pa_type, nullable=True, metadata=_build_field_metadata(xsd_type=xsd_type))
        return pa.field(
            name,
            pa.list_(item_field),
            nullable=nullable,
            metadata=_build_field_metadata(comment=comment),
        )
    return pa.field(
        name,
        pa_type,
        nullable=nullable,
        metadata=_build_field_metadata(xsd_type=xsd_type, comment=comment),
    )


def _resolve_target_element(schema: xmlschema.XMLSchema, xml_tag: str) -> XsdElement:
    """Resolve ``xml_tag`` (local name or Clark notation) to a unique ``XsdElement`` declaration."""
    if xml_tag.startswith("{") and "}" in xml_tag:
        ns_filter: str | None = xml_tag
        local_target = xml_tag.rsplit("}", maxsplit=1)[-1]
    else:
        ns_filter = None
        local_target = xml_tag

    def _matches(el: XsdElement) -> bool:
        if el.local_name != local_target:
            return False
        return ns_filter is None or el.name == ns_filter

    # Depth 0: global elements
    depth0 = [el for el in schema.elements.values() if isinstance(el, XsdElement) and _matches(el)]
    if depth0:
        return _deduplicate_target_elements(depth0, xml_tag)

    # BFS from depth 1 onwards through complexType model groups
    current_level: list[XsdElement] = [el for el in schema.elements.values() if isinstance(el, XsdElement)]
    visited_types: set[int] = set()
    while current_level:
        matches: list[XsdElement] = []
        next_level: list[XsdElement] = []
        for parent in current_level:
            p_type = parent.type
            if not isinstance(p_type, XsdComplexType) or id(p_type) in visited_types:
                continue
            visited_types.add(id(p_type))
            if p_type.model_group is None:
                continue
            for child in p_type.model_group.iter_elements():
                if not isinstance(child, XsdElement):
                    continue
                if _matches(child):
                    matches.append(child)
                next_level.append(child)
        if matches:
            return _deduplicate_target_elements(matches, xml_tag)
        current_level = next_level

    msg = f"Element {xml_tag!r} not found in XSD schema"
    raise ValueError(msg)


def _deduplicate_target_elements(matches: list[XsdElement], xml_tag: str) -> XsdElement:
    """Verify all matching ``XsdElement`` declarations at the resolution depth share a single ``XsdType``."""
    by_type: dict[int, XsdElement] = {id(el.type): el for el in matches}
    if len(by_type) > 1:
        msg = f"Ambiguous xml_tag {xml_tag!r}: found {len(by_type)} declarations with different types"
        raise ValueError(msg)
    target = next(iter(by_type.values()))
    target_ref = target.ref if target.ref is not None else target
    if getattr(target, "abstract", False) or getattr(target_ref, "abstract", False):
        msg = f"Unsupported XSD feature: abstract element {xml_tag!r} is not supported"
        raise ValueError(msg)
    if target.type is not None and getattr(target.type, "abstract", False):
        msg = f"Unsupported XSD feature: abstract type on {xml_tag!r} is not supported"
        raise ValueError(msg)
    return target


def _require_complex_type(elem_type: XsdType | None, xml_tag: str) -> XsdComplexType:
    """Ensure that the resolved target element has an ``XsdComplexType``."""
    if elem_type is None or not elem_type.is_complex():
        msg = f"Top-level xml_tag {xml_tag!r} must have a complex type with attributes or child elements"
        raise ValueError(msg)
    assert isinstance(elem_type, XsdComplexType)  # noqa: S101
    return elem_type


def _prepare_target_traversal(
    schema: xmlschema.XMLSchema,
    xml_tag: str,
    *,
    max_recursion_depth: int,
) -> tuple[XsdElement, XsdComplexType, _TraversalContext]:
    """Resolve ``xml_tag`` to its ``XsdElement``, validated ``XsdComplexType``, and initial ``_TraversalContext``."""
    target_elem = _resolve_target_element(schema, xml_tag)
    elem_type = _require_complex_type(target_elem.type, xml_tag)
    ctx = _TraversalContext(
        schema=schema,
        type_counts={id(elem_type): 1},
        max_recursion_depth=max_recursion_depth,
    )
    return target_elem, elem_type, ctx


def generate_pyarrow_schema(
    schema: xmlschema.XMLSchema,
    xml_tag: str,
    *,
    max_recursion_depth: int = DEFAULT_MAX_RECURSION_DEPTH,
) -> pa.Schema:
    """Generate a PyArrow schema for the target XML element from an XSD schema.

    :param schema: loaded XMLSchema instance
    :type  schema: xmlschema.XMLSchema
    :param xml_tag: target XML element name (local name or Clark-notation ``{ns}local``)
    :type  xml_tag: str
    :param max_recursion_depth: maximum structural recursion depth before JSON serialization
    :type  max_recursion_depth: int
    :return: PyArrow schema whose top-level fields correspond to the target element's fields
    :rtype : pa.Schema
    """
    target_elem, elem_type, ctx = _prepare_target_traversal(
        schema,
        xml_tag,
        max_recursion_depth=max_recursion_depth,
    )

    if _is_empty_marker_complex_type(elem_type):
        msg = f"Top-level xml_tag {xml_tag!r} has an empty complex type with no attributes or child elements"
        raise ValueError(msg)

    attrs = _collect_attributes(elem_type)
    if (elem_type.has_simple_content() or elem_type.mixed) and not attrs:
        msg = f"Top-level xml_tag {xml_tag!r} must have a complex type that resolves to a struct"
        raise ValueError(msg)

    fields = _build_complex_struct_fields(
        ctx,
        elem_type,
        parent_nillable=target_elem.nillable,
    )
    if not fields:
        msg = f"Top-level xml_tag {xml_tag!r} must have a complex type with at least one field"
        raise ValueError(msg)
    return pa.schema(fields)


def _is_ignored_attr(attr_name: str) -> bool:
    """Return True if an XML attribute is an ``xmlns`` or ``xsi:*`` attribute that should be stripped."""
    return attr_name.startswith((XSI_NAMESPACE, "xsi:", "xmlns", "{http://www.w3.org/2000/xmlns/}"))


def _is_xsi_nil(element: Element) -> bool:
    """Return True if ``element`` carries ``xsi:nil='true'`` or ``xsi:nil='1'``."""
    val = element.attrib.get(XSI_NIL) or element.attrib.get("xsi:nil")
    return val is not None and val.strip() in {"true", "1"}


def _escape_xml_text(text: str) -> str:
    """Escape XML character entities ``&``, ``<``, and ``>`` in mixed-content text."""
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _escape_xml_attr(value: str) -> str:
    """Escape XML character entities and double quotes in mixed-content attribute values."""
    return _escape_xml_text(value).replace('"', "&quot;")


def _serialize_mixed_node(node: Element) -> str:
    """Serialize an inline child element inside a ``mixed='true'`` element without namespace declarations."""
    local = _local_name(node.tag)
    attrs = "".join(
        f' {_local_name(k)}="{_escape_xml_attr(v)}"' for k, v in node.attrib.items() if not _is_ignored_attr(k)
    )
    inner = _serialize_inner_xml(node)
    if not inner and len(node) == 0 and node.text is None:
        return f"<{local}{attrs}/>"
    return f"<{local}{attrs}>{inner}</{local}>"


def _serialize_inner_xml(element: Element) -> str:
    """Serialize the raw inner XML (text and inline child elements) of a ``mixed='true'`` element."""
    parts: list[str] = []
    if element.text:
        parts.append(_escape_xml_text(element.text))
    for child in element:
        if isinstance(child.tag, str):
            parts.append(_serialize_mixed_node(child))
        if child.tail:
            parts.append(_escape_xml_text(child.tail))
    return "".join(parts)


def _coerce_resolved_scalar(stripped: str, arrow_type: pa.DataType, xsd_name: str) -> Any:  # noqa: ANN401
    """Coerce a stripped lexical string using its resolved ``arrow_type`` and ``xsd_name``."""
    if pa.types.is_boolean(arrow_type):
        return stripped in {"true", "1"}
    if pa.types.is_integer(arrow_type):
        return int(stripped)
    if pa.types.is_floating(arrow_type):
        return float(stripped)
    if pa.types.is_date32(arrow_type):
        return dt.date.fromisoformat(stripped[:10])
    if pa.types.is_binary(arrow_type):
        return bytes.fromhex(stripped) if xsd_name == "xs:hexBinary" else base64.b64decode(stripped)
    return stripped


def _coerce_scalar_text(raw_text: str, simple_type: XsdType | None) -> Any:  # noqa: ANN401
    """Coerce a non-empty lexical XML string according to its XSD ``simple_type``."""
    if simple_type is None:
        return raw_text

    builtin_name = _builtin_xsd_name(simple_type)
    if builtin_name in BUILTIN_LIST_ITEM_TYPES:
        return raw_text.split()

    if simple_type.is_list():
        item_type = _unwrap_item_type(simple_type)
        return [_coerce_scalar_text(token, item_type) for token in raw_text.split()]

    arrow_type, xsd_name = _resolve_scalar_type(simple_type)
    return _coerce_resolved_scalar(raw_text.strip(), arrow_type, xsd_name)


def _extract_struct_attributes(element: Element, complex_type: XsdComplexType) -> dict[str, Any]:
    """Extract and coerce declared XML attributes from ``element`` into ``_<name>`` dictionary entries."""
    declared_attrs = {attr.local_name: attr for attr in _collect_attributes(complex_type) if attr.local_name}
    raw_attrs = {_local_name(k): v for k, v in element.attrib.items() if not _is_ignored_attr(k)}
    return {
        f"_{attr_name}": _coerce_scalar_text(raw_attrs[attr_name], attr_decl.type)
        for attr_name, attr_decl in declared_attrs.items()
        if attr_name in raw_attrs
    }


def _group_child_elements_by_name(element: Element) -> dict[str, list[Element]]:
    """Group direct XML child elements by their local tag name in document order."""
    children_by_name: dict[str, list[Element]] = {}
    for child_el in element:
        if isinstance(child_el.tag, str):
            children_by_name.setdefault(_local_name(child_el.tag), []).append(child_el)
    return children_by_name


def _convert_complex_struct(
    ctx: _TraversalContext,
    element: Element,
    complex_type: XsdComplexType,
) -> dict[str, Any] | None:
    """Convert an XML ``element`` with a struct-mapped ``complex_type`` into a sparse dictionary."""
    result = _extract_struct_attributes(element, complex_type)

    if _is_xsi_nil(element):
        return result or None

    if complex_type.has_simple_content():
        text_val = (element.text or "").strip()
        if text_val:
            content_type = complex_type.content if isinstance(complex_type.content, XsdType) else None
            result["#text"] = _coerce_scalar_text(text_val, content_type)
        return result or None

    if complex_type.mixed:
        inner_xml = _serialize_inner_xml(element).strip()
        if inner_xml:
            result["#text"] = inner_xml
        return result or None

    model_group = complex_type.model_group
    if model_group is None:
        return result or None

    children_by_name = _group_child_elements_by_name(element)

    for child_decl in _collect_unique_children(ctx.schema, complex_type):
        child_name = child_decl.local_name
        matching_elements = children_by_name.get(child_name)
        if not matching_elements:
            continue

        is_list, nullable = _resolve_child_cardinality(model_group, child_decl)

        if is_list:
            result[child_name] = [
                _convert_child_occurrence(ctx, el, child_decl, nullable=True) for el in matching_elements
            ]
        else:
            val = _convert_child_occurrence(ctx, matching_elements[0], child_decl, nullable=nullable)
            if val is not None:
                result[child_name] = val

    return result or None


def _convert_text_or_empty_leaf(
    element: Element,
    *,
    is_mixed: bool,
    simple_type: XsdType | None,
    nullable: bool,
) -> Any:  # noqa: ANN401
    """Convert an XML element that maps to a scalar or mixed-text leaf field."""
    if _is_xsi_nil(element):
        return None
    raw_val = _serialize_inner_xml(element).strip() if is_mixed else (element.text or "").strip()
    if not raw_val:
        return "" if not nullable else None
    return raw_val if is_mixed else _coerce_scalar_text(raw_val, simple_type)


def _convert_child_occurrence(
    ctx: _TraversalContext,
    element: Element,
    xsd_element: XsdElement,
    *,
    nullable: bool,
) -> Any:  # noqa: ANN401
    """Convert a single XML element occurrence according to ``xsd_element``."""
    elem_type = xsd_element.type
    leaf_info = _classify_leaf_element_type(elem_type)
    if leaf_info is not None:
        kind, simple_type = leaf_info
        if kind == "empty_marker":
            return None if _is_xsi_nil(element) else ""
        return _convert_text_or_empty_leaf(
            element,
            is_mixed=(kind == "mixed"),
            simple_type=simple_type,
            nullable=nullable,
        )

    assert isinstance(elem_type, XsdComplexType)  # noqa: S101
    type_id = id(elem_type)
    next_count = ctx.type_counts.get(type_id, 0) + 1
    if ctx.max_recursion_depth is not None and next_count > ctx.max_recursion_depth:
        unbounded_ctx = _TraversalContext(schema=ctx.schema, type_counts={}, max_recursion_depth=None)
        nested = _convert_complex_struct(unbounded_ctx, element, elem_type)
        return json.dumps(nested) if nested is not None else None

    next_ctx = _TraversalContext(
        schema=ctx.schema,
        type_counts={**ctx.type_counts, type_id: next_count},
        max_recursion_depth=ctx.max_recursion_depth,
    )
    return _convert_complex_struct(next_ctx, element, elem_type)


def convert_xml_element(
    element: Element,
    schema: xmlschema.XMLSchema,
    xml_tag: str,
    *,
    max_recursion_depth: int = DEFAULT_MAX_RECURSION_DEPTH,
) -> dict[str, Any]:
    """Convert an lxml Element matching ``xml_tag`` into a typed dictionary according to the XSD.

    :param element: XML element to convert
    :type  element: Element
    :param schema: loaded XMLSchema instance
    :type  schema: xmlschema.XMLSchema
    :param xml_tag: target XML element name (local name or Clark-notation ``{ns}local``)
    :type  xml_tag: str
    :param max_recursion_depth: maximum structural recursion depth before JSON serialization
    :type  max_recursion_depth: int
    :return: sparse dictionary of coerced field values matching the generated PyArrow schema
    :rtype : dict[str, Any]
    """
    _, elem_type, ctx = _prepare_target_traversal(
        schema,
        xml_tag,
        max_recursion_depth=max_recursion_depth,
    )
    converted = _convert_complex_struct(ctx, element, elem_type)
    return converted or {}
