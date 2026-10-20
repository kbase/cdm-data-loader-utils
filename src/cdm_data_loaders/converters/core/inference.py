"""Target-independent JSON Schema type inference and decimal scale helpers."""

from decimal import Decimal
from typing import Any, Final, Literal

from cdm_data_loaders.converters.core.ir import TypedNode

IMPLICIT_OBJECT_KEYWORDS: Final[frozenset[str]] = frozenset(
    {
        "properties",
        "patternProperties",
        "additionalProperties",
        "unevaluatedProperties",
        "required",
        "propertyNames",
        "minProperties",
        "maxProperties",
        "dependentSchemas",
        "dependentRequired",
        "dependencies",
    }
)
IMPLICIT_ARRAY_KEYWORDS: Final[frozenset[str]] = frozenset(
    {
        "items",
        "prefixItems",
        "additionalItems",
        "unevaluatedItems",
        "contains",
        "minItems",
        "maxItems",
        "uniqueItems",
        "minContains",
        "maxContains",
    }
)
IMPLICIT_STRING_KEYWORDS: Final[frozenset[str]] = frozenset(
    {
        "pattern",
        "minLength",
        "maxLength",
        "format",
        "contentEncoding",
        "contentMediaType",
        "contentSchema",
    }
)
IMPLICIT_NUMBER_KEYWORDS: Final[frozenset[str]] = frozenset(
    {"minimum", "maximum", "exclusiveMinimum", "exclusiveMaximum", "multipleOf"}
)


def infer_implicit_type(schema: dict[str, Any]) -> str | None:
    """Approximate a type-less schema from its type-specific keywords.

    Precedence is object, array, string, number. These keywords constrain
    values of their own type but do not exclude other JSON types.
    """
    if schema.keys() & IMPLICIT_OBJECT_KEYWORDS:
        return "object"
    if schema.keys() & IMPLICIT_ARRAY_KEYWORDS:
        return "array"
    if schema.keys() & IMPLICIT_STRING_KEYWORDS:
        return "string"
    if schema.keys() & IMPLICIT_NUMBER_KEYWORDS:
        return "number"
    return None


def json_type_from_enum(values: list[Any]) -> Literal["boolean", "integer", "number", "string"]:
    """Infer a scalar JSON type, using string for empty or heterogeneous enums."""
    if not values:
        return "string"
    if all(isinstance(value, bool) for value in values):
        return "boolean"
    if all(isinstance(value, int) and not isinstance(value, bool) for value in values):
        return "integer"
    if all(isinstance(value, (int, float)) and not isinstance(value, bool) for value in values):
        return "number"
    return "string"


def is_unconstrained_node(node: TypedNode) -> bool:
    """Identify false or empty JSON schemas that carry no type information.

    :param node: the schema node to inspect
    :type node: TypedNode
    :return: True for a ``never`` node or an unknown node with no keywords, metadata or combiners
    :rtype: bool
    """
    if node.type == "never":
        return True
    return (
        node.type == "unknown"
        and node.declared_type is None
        and node.inferred_type is None
        and not node.source_keywords
        and not node.constraints
        and not node.annotations
        and not node.extensions
        and node.one_of is None
        and node.any_of is None
        and not node.schema_keywords
        and not node.schema_maps
    )


def decimal_places(value: float | Decimal) -> int:
    """Count fractional digits in the decimal representation; non-finite values return zero."""
    exponent = Decimal(str(value)).as_tuple().exponent
    if isinstance(exponent, int):
        return max(-exponent, 0)
    return 0


def decimal_pattern(precision: int, scale: int) -> str:
    """Build the decimal string pattern from integer precision and scale."""
    digits_before = precision - scale
    if scale > 0:
        if digits_before == 0:
            return rf"^-?0(\.\d{{1,{scale}}})?$"
        return rf"^-?\d{{1,{digits_before}}}(\.\d{{1,{scale}}})?$"
    return rf"^-?\d{{1,{digits_before}}}$"


def resolve_union_type(node: TypedNode, fallback: str = "text") -> str:
    """Collapse a union or combiner into a single best-fit scalar type.

    Precedence:
    1. Explicit node type
    2. explicit declared_type (if single non-null)
    3. enum-based inference
    4. inferred_type
    5. first branch of oneOf
    6. first branch of anyOf
    7. fallback
    """
    if node.type not in {"unknown", "any", "never"}:
        return node.type
    declaration = node.declared_type
    if isinstance(declaration, tuple):
        non_null = tuple(kind for kind in declaration if kind != "null")
        if len(non_null) == 1:
            return non_null[0]

    if "enum" in node.constraints:
        enum_values = node.constraints["enum"]
        if isinstance(enum_values, tuple):
            return json_type_from_enum(list(enum_values))

    if node.inferred_type is not None:
        return node.inferred_type

    for branches in (node.one_of, node.any_of):
        if branches:
            # Recurse into the first branch to see if it's a simple type
            # This is a simplification; real logic might need to check if the branch is a scalar
            first = branches[0]
            if first.type not in {"object", "array", "map"}:
                return first.type

    return declaration if isinstance(declaration, str) else fallback
