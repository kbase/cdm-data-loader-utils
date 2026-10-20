"""Emit PyArrow schema objects from immutable typed schema nodes."""

import json
import logging
from collections.abc import Mapping
from decimal import Decimal
from typing import Final

import pyarrow as pa
from frozendict import frozendict
from jsonschema import Draft202012Validator
from pydantic import BaseModel, ConfigDict, Field, field_validator

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.inference import decimal_places, is_unconstrained_node, json_type_from_enum
from cdm_data_loaders.converters.core.ir import Field as SchemaField
from cdm_data_loaders.converters.core.ir import NodeHints, SchemaDocument, TypedNode
from cdm_data_loaders.converters.core.ir_values import mutable_value
from cdm_data_loaders.converters.core.metadata import DEFAULT_METADATA_KEYWORDS, ConversionContext

logger = logging.getLogger(__name__)

INT8_WIDTH: Final = 8
INT16_WIDTH: Final = 16
INT32_WIDTH: Final = 32
INT32_MIN: Final = -2_147_483_648
INT32_MAX: Final = 2_147_483_647
DECIMAL128_MAX_PRECISION: Final = 38
DECIMAL256_MAX_PRECISION: Final = 76
WEI_PRECISION: Final = DECIMAL256_MAX_PRECISION
MAX_MILLIS_DIGITS: Final = 3
MAX_MICROS_DIGITS: Final = 6
DEFAULT_FORMAT_MAP: Final[frozendict[str, pa.DataType]] = frozendict(
    {
        "date": pa.date32(),
        "date-time": pa.timestamp("us", tz="UTC"),
        **dict.fromkeys(
            (
                "time",
                "duration",
                "email",
                "hostname",
                "idn-email",
                "idn-hostname",
                "ipv4",
                "ipv6",
                "iri-reference",
                "iri",
                "uri-reference",
                "uri",
                "uuid",
            ),
            pa.string(),
        ),
    }
)


def merge_format_map(value: object) -> object:
    """Validate format overrides and merge them into the built-in mapping.

    :param value: a mapping from JSON Schema format names to PyArrow types
    :type value: object
    :return: the merged mapping, or ``value`` unchanged when it is not a mapping
    :rtype: object
    :raises ValueError: if a key is not a string or a value is not a PyArrow ``DataType``
    """
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) or not isinstance(item, pa.DataType) for key, item in value.items()):
            msg = "format_map must map strings to PyArrow DataType instances"
            raise ValueError(msg)
        return frozendict({**DEFAULT_FORMAT_MAP, **value})
    return value


class PyArrowEmitterError(ConversionError):
    """A typed schema cannot be represented by the configured Arrow policy."""


class PyArrowEmitter(BaseModel):
    """Convert typed schema nodes into PyArrow schemas and empty tables."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    format_map: frozendict[str, pa.DataType] = Field(default_factory=lambda: DEFAULT_FORMAT_MAP)
    treat_unknown_as_string: bool = True
    emit_unions_as_json: bool = False
    extra_metadata_keywords: frozenset[str] = Field(default_factory=frozenset)
    max_nesting: int = Field(default=10, ge=1)

    @field_validator("format_map", mode="before")
    @classmethod
    def _merge_with_builtin_format_map(cls, value: object) -> object:
        """Merge nonempty format overrides into the built-in mapping."""
        return merge_format_map(value) if value else value

    def build_context(self, validator_cls: type) -> ConversionContext:
        """Build draft-specific metadata policy and report invalid requested keywords.

        :param validator_cls: the ``jsonschema`` validator class for the source draft
        :type validator_cls: type
        :return: the metadata selection context
        :rtype: ConversionContext
        """
        context = ConversionContext(validator_cls, self.extra_metadata_keywords)
        if context.invalid_extra_metadata_keywords:
            logger.warning(
                "Ignoring invalid extra_metadata_keywords %s: these are neither recognised JSON Schema "
                "keywords for the detected draft (%s) nor vendor extensions prefixed with 'x-'.",
                sorted(context.invalid_extra_metadata_keywords),
                validator_cls.__name__,
            )
        return context

    def emit(self, document: SchemaDocument) -> pa.Schema:
        """Emit a table schema, requiring the root to be a struct.

        :param document: the schema document to convert
        :type document: SchemaDocument
        :return: the Arrow schema
        :rtype: pa.Schema
        :raises PyArrowEmitterError: if the root does not resolve to a struct
        """
        validator_cls = document.jsonschema_validator_cls or Draft202012Validator
        result = self.emit_node(document.root, self.build_context(validator_cls))
        if pa.types.is_struct(result):
            return pa.schema([result.field(index) for index in range(result.num_fields)])
        msg = (
            "Root JSON Schema did not resolve to a struct.\nThis can happen if it has no "
            "'properties' but does have 'patternProperties' or 'additionalProperties', which map to "
            "a map type instead. An Arrow table requires a struct at the top level of its schema."
        )
        raise PyArrowEmitterError(msg)

    def emit_table(self, document: SchemaDocument) -> pa.Table:
        """Emit an empty table with the document's schema.

        :param document: the schema document to convert
        :type document: SchemaDocument
        :return: a zero-row table
        :rtype: pa.Table
        :raises PyArrowEmitterError: if the root does not resolve to a struct
        """
        return self.emit(document).empty_table()

    def emit_node(self, node: TypedNode, ctx: ConversionContext | None = None, depth: int = 0) -> pa.DataType:
        """Emit a scalar or container node without a table-root restriction.

        :param node: the node to convert
        :type node: TypedNode
        :param ctx: draft-specific metadata context; Draft 2020-12 when omitted
        :type ctx: ConversionContext | None
        :param depth: current nesting depth
        :type depth: int
        :return: the Arrow type
        :rtype: pa.DataType
        :raises PyArrowEmitterError: if the node cannot be represented under the configured policy
        """
        if ctx is None:
            ctx = self.build_context(Draft202012Validator)
        logical = self._emit_logical(node)
        if logical is not None:
            return logical
        kind = self._type_name(node)
        if kind in {"object", "array", "map"}:
            return {"object": self._emit_object, "array": self._emit_array, "map": self._emit_map}[kind](
                node, ctx, depth
            )
        scalars = {"boolean": pa.bool_(), "null": pa.null(), "collapsed-union": pa.string()}
        if kind in scalars:
            return scalars[kind]
        if kind == "enum":
            values = node.constraints["enum"]
            if isinstance(values, tuple):
                return infer_type_from_enum(list(values))
        converters = {"string": self._emit_string, "integer": self._emit_integer, "number": self._emit_number}
        return converters[kind](node) if kind in converters else self._emit_unknown(node, ctx)

    def build_metadata(
        self, node: TypedNode, ctx: ConversionContext, field: SchemaField | None = None
    ) -> dict[str, str]:
        """Select field metadata and encode it as Arrow's string key/value pairs.

        :param node: the node whose constraints, annotations and extensions are searched
        :type node: TypedNode
        :param ctx: draft-specific metadata context
        :type ctx: ConversionContext
        :param field: the owning field, whose annotations and extensions take precedence
        :type field: SchemaField | None
        :return: ``jsonschema`` (JSON object), ``comment`` and ``original_union`` (JSON array) entries
        :rtype: dict[str, str]
        """
        values = node.aggregate_metadata(field)
        selected = DEFAULT_METADATA_KEYWORDS | ctx.allowed_extra_metadata_keywords
        jsonschema = {key: mutable_value(values[key]) for key in selected if key in values}
        metadata: dict[str, str] = {}
        if jsonschema:
            metadata["jsonschema"] = _dump_json(jsonschema)
        description = values.get("description")
        if isinstance(description, str):
            metadata["comment"] = description
        if self.emit_unions_as_json:
            declaration = node.declared_type
            if isinstance(declaration, tuple):
                non_null = [kind for kind in declaration if kind != "null"]
                if len(non_null) > 1:
                    metadata["original_union"] = _dump_json(non_null)
        return metadata

    def _type_name(self, node: TypedNode) -> str:
        """Select declared, enum or inferred types without flattening the IR."""
        declaration = node.declared_type
        if isinstance(declaration, tuple):
            non_null = [kind for kind in declaration if kind != "null"]
            if not non_null:
                return "null"
            if len(non_null) == 1:
                return non_null[0]
            if not self.treat_unknown_as_string:
                msg = f"Unsupported multi-type union: {list(declaration)!r}"
                raise PyArrowEmitterError(msg)
            if self.emit_unions_as_json:
                logger.info("Emitting multi-type union %r as JSON string.", list(declaration))
            else:
                logger.warning("Collapsing multi-type union %r to string.", list(declaration))
            return "collapsed-union"
        if declaration is not None:
            return declaration
        if "enum" in node.constraints:
            return "enum"
        return node.inferred_type or node.type

    def _emit_unknown(self, node: TypedNode, ctx: ConversionContext) -> pa.DataType:
        """Approximate combiners or apply strict unsupported-type handling."""
        if node.type in {"any", "never"} and (node.provenance is None or node.provenance.source == "json-schema"):
            message = (
                f"Boolean JSON Schema `{node.type == 'any'}` has no Arrow equivalent "
                "(unconstrained 'any' type or unsatisfiable 'never' type)."
            )
        else:
            for keyword, branches in (("oneOf", node.one_of), ("anyOf", node.any_of)):
                if branches:
                    logger.warning(
                        "Approximating '%s' by using only its first branch for type inference; Arrow has no union "
                        "type compatible with table storage.",
                        keyword,
                    )
                    return self.emit_node(branches[0], ctx)
            for keyword in ("not", "if", "then", "else"):
                if keyword in node.schema_keywords:
                    logger.warning(
                        "Ignoring unsupported conditional keyword '%s'; it has no effect on the resulting Arrow type.",
                        keyword,
                    )
            message = f"Unsupported/unknown schema type: {node.declared_type or node.type!r}"
        if not self.treat_unknown_as_string:
            raise PyArrowEmitterError(message)
        logger.warning("%s Falling back to string.", message)
        return pa.string()

    def _emit_object(self, node: TypedNode, ctx: ConversionContext, depth: int = 0) -> pa.DataType:
        """Emit fixed fields, a dynamic map, or an empty struct."""
        if depth >= self.max_nesting:
            logger.warning("Max nesting limit (%d) reached; collapsing object to string.", self.max_nesting)
            return pa.string()

        additional = node.additional_properties
        schema_additional = additional is not None and additional.type not in {"any", "never"}
        if node.properties:
            if node.pattern_properties or schema_additional:
                logger.warning(
                    "Object schema declares both fixed 'properties' and dynamic 'patternProperties' / "
                    "'additionalProperties'; an Arrow struct requires fixed field names, so the "
                    "dynamic keys are ignored."
                )
            return pa.struct(
                [
                    pa.field(
                        prop.name,
                        self.emit_node(prop.node, ctx, depth + 1),
                        nullable=not prop.required,
                        metadata=self.build_metadata(prop.node, ctx, prop) or None,
                    )
                    for prop in node.properties
                ]
            )
        value_node = additional if schema_additional else None
        if node.pattern_properties:
            if len(node.pattern_properties) > 1:
                logger.warning(
                    "Multiple 'patternProperties' patterns found; only the first pattern's schema is "
                    "used as the map value type, since an Arrow map has a single value type."
                )
            value_node = next(iter(node.pattern_properties.values()))
        if value_node is not None:
            return pa.map_(pa.string(), pa.field("value", self.emit_node(value_node, ctx, depth + 1), nullable=True))
        logger.warning(
            "Object schema has no 'properties', 'patternProperties', or schema-valued "
            "'additionalProperties'; converting to an empty struct."
        )
        return pa.struct([])

    def _emit_map(self, node: TypedNode, ctx: ConversionContext, depth: int = 0) -> pa.MapType:
        """Emit a map with independently typed keys and values."""
        if node.key_type is None or node.value_type is None:
            msg = "Map nodes require both key_type and value_type"
            raise PyArrowEmitterError(msg)
        value = pa.field(
            "value", self.emit_node(node.value_type, ctx, depth + 1), nullable=_element_nullable(node.value_type)
        )
        return pa.map_(self.emit_node(node.key_type, ctx, depth + 1), value)

    def _emit_array(self, node: TypedNode, ctx: ConversionContext, depth: int = 0) -> pa.ListType:
        """Emit homogeneous arrays, approximating tuple schemas by their first slot."""
        element = node.items
        if node.prefix_items is not None:
            logger.warning(
                "Array schema uses 'prefixItems' (tuple validation); Arrow lists are homogeneous, so "
                "only the first tuple slot's type is used as the list's element type."
            )
            element = node.prefix_items[0] if node.prefix_items else None
        elif isinstance(element, tuple):
            logger.warning(
                "Array schema uses tuple-style 'items' (list of schemas); Arrow lists are "
                "homogeneous, so only the first item's type is used as the list's element type."
            )
            element = element[0] if element else None
        if element is None or is_unconstrained_node(element):
            return pa.list_(pa.field("item", pa.string(), nullable=True))
        return pa.list_(pa.field("item", self.emit_node(element, ctx, depth + 1), nullable=_element_nullable(element)))

    def _emit_string(self, node: TypedNode) -> pa.DataType:
        """Map recognized formats to Arrow types."""
        value = node.constraints.get("format")
        return self.format_map.get(value, pa.string()) if isinstance(value, str) else pa.string()

    @staticmethod
    def _emit_integer(node: TypedNode) -> pa.DataType:
        """Use int32 only when both bounds fit."""
        minimum = node.constraints.get("minimum", node.constraints.get("exclusiveMinimum"))
        maximum = node.constraints.get("maximum", node.constraints.get("exclusiveMaximum"))
        if (
            isinstance(minimum, (int, float, Decimal))
            and isinstance(maximum, (int, float, Decimal))
            and minimum >= INT32_MIN
            and maximum <= INT32_MAX
        ):
            return pa.int32()
        return pa.int64()

    @staticmethod
    def _emit_number(node: TypedNode) -> pa.DataType:
        """Derive decimal scale from numeric JSON multipleOf values."""
        multiple = node.constraints.get("multipleOf")
        if isinstance(multiple, (int, float, Decimal)):
            scale = decimal_places(multiple)
            if scale > DECIMAL128_MAX_PRECISION:
                msg = f"Unsupported Arrow decimal precision/scale: ({DECIMAL128_MAX_PRECISION}, {scale})"
                raise PyArrowEmitterError(msg)
            return pa.decimal128(DECIMAL128_MAX_PRECISION, scale)
        return pa.float64()

    @staticmethod
    def _emit_logical(node: TypedNode) -> pa.DataType | None:
        """Preserve physical widths and logical source types without a JSON encoding."""
        if node.provenance is not None and node.provenance.source == "json-schema":
            return None
        hints = node.hints
        logical = hints.logical_type
        if logical in {"decimal", "wei"}:
            return _decimal_type(hints, wei=logical == "wei")
        if logical in {"timestamp", "timestamp-with-tz", "timestamp-without-tz"}:
            unit = _timestamp_unit(hints.precision)
            naive = hints.timezone is False or logical == "timestamp-without-tz"
            return pa.timestamp(unit) if naive else pa.timestamp(unit, tz="UTC")
        if node.type == "integer" and hints.bit_width is not None:
            return _integer_type(hints.bit_width)
        if node.type == "number" and (logical == "float" or hints.bit_width == INT32_WIDTH):
            return pa.float32()
        fixed = pa.binary(hints.length) if hints.length is not None else pa.binary()
        return {"date": pa.date32(), "time": pa.time64("us"), "binary": pa.binary(), "fixed": fixed}.get(logical)


def infer_type_from_enum(values: list[object]) -> pa.DataType:
    """Map homogeneous JSON enum values to Arrow scalar types.

    :param values: the enum members
    :type values: list[object]
    :return: ``bool_``, ``int64``, ``float64`` or ``string``
    :rtype: pa.DataType
    """
    return {
        "boolean": pa.bool_(),
        "integer": pa.int64(),
        "number": pa.float64(),
        "string": pa.string(),
    }[json_type_from_enum(values)]


def _element_nullable(node: TypedNode) -> bool:
    """Honor non-null element facts from typed sources; JSON item schemas stay nullable."""
    from_json = node.provenance is None or node.provenance.source == "json-schema"
    return from_json or node.nullable is not False


def _decimal_type(hints: NodeHints, *, wei: bool) -> pa.DataType:
    """Select decimal128 or decimal256 from precision and scale hints."""
    default_precision = WEI_PRECISION if wei else DECIMAL128_MAX_PRECISION
    precision = hints.precision if hints.precision is not None else default_precision
    scale = hints.scale if hints.scale is not None else 0
    if not 1 <= precision <= DECIMAL256_MAX_PRECISION or scale > precision:
        msg = f"Unsupported Arrow decimal precision/scale: ({precision}, {scale})"
        raise PyArrowEmitterError(msg)
    return pa.decimal128(precision, scale) if precision <= DECIMAL128_MAX_PRECISION else pa.decimal256(precision, scale)


def _timestamp_unit(precision: int | None) -> str:
    """Convert a count of fractional second digits to an Arrow time unit."""
    if precision is None:
        return "us"
    if precision == 0:
        return "s"
    if precision <= MAX_MILLIS_DIGITS:
        return "ms"
    return "us" if precision <= MAX_MICROS_DIGITS else "ns"


def _integer_type(bit_width: int) -> pa.DataType:
    """Select the narrowest signed Arrow integer holding the given bit width."""
    if bit_width <= INT8_WIDTH:
        return pa.int8()
    if bit_width <= INT16_WIDTH:
        return pa.int16()
    return pa.int32() if bit_width <= INT32_WIDTH else pa.int64()


def _json_default(value: object) -> str:
    """Encode exact decimals as strings; reject all other unsupported objects."""
    if isinstance(value, Decimal):
        return str(value)
    msg = f"Object of type {type(value).__name__} is not JSON serializable"
    raise TypeError(msg)


def _dump_json(value: object) -> str:
    """Serialize metadata without converting decimals to binary floats."""
    return json.dumps(value, default=_json_default)
