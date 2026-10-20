"""Emit PySpark schema objects from immutable typed schema nodes."""

import logging
from collections.abc import Mapping
from decimal import Decimal
from typing import Final

from frozendict import frozendict
from jsonschema import Draft202012Validator
from pydantic import BaseModel, ConfigDict, Field, field_validator
from pyspark.sql.types import (
    ArrayType,
    BinaryType,
    BooleanType,
    DataType,
    DateType,
    DecimalType,
    DoubleType,
    FloatType,
    IntegerType,
    LongType,
    MapType,
    NullType,
    StringType,
    StructField,
    StructType,
    TimestampNTZType,
    TimestampType,
)

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.inference import decimal_places, is_unconstrained_node, json_type_from_enum
from cdm_data_loaders.converters.core.ir import Field as SchemaField
from cdm_data_loaders.converters.core.ir import SchemaDocument, TypedNode
from cdm_data_loaders.converters.core.ir_values import mutable_value
from cdm_data_loaders.converters.core.metadata import DEFAULT_METADATA_KEYWORDS, ConversionContext

logger = logging.getLogger(__name__)

INT32_MIN: Final = -2_147_483_648
INT32_MAX: Final = 2_147_483_647
INT32_WIDTH: Final = 32
MAX_DECIMAL_PRECISION: Final = 38
DEFAULT_FORMAT_MAP: Final[frozendict[str, DataType]] = frozendict(
    {
        "date": DateType(),
        "date-time": TimestampType(),
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
            StringType(),
        ),
    }
)


def merge_format_map(value: object) -> object:
    """Validate format overrides and merge them into the built-in mapping."""
    if isinstance(value, Mapping):
        if any(not isinstance(key, str) or not isinstance(item, DataType) for key, item in value.items()):
            msg = "format_map must map strings to PySpark DataType instances"
            raise ValueError(msg)
        return frozendict({**DEFAULT_FORMAT_MAP, **value})
    return value


class PySparkEmitterError(ConversionError):
    """A typed schema cannot be represented by the configured Spark policy."""


class PySparkEmitter(BaseModel):
    """Convert typed schema nodes into pure Python Spark schema objects."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    format_map: frozendict[str, DataType] = Field(default_factory=lambda: DEFAULT_FORMAT_MAP)
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
        """Build draft-specific metadata policy and report invalid requested keywords."""
        context = ConversionContext(validator_cls, self.extra_metadata_keywords)
        if context.invalid_extra_metadata_keywords:
            logger.warning(
                "Ignoring invalid extra_metadata_keywords %s: these are neither recognised JSON Schema "
                "keywords for the detected draft (%s) nor vendor extensions prefixed with 'x-'.",
                sorted(context.invalid_extra_metadata_keywords),
                validator_cls.__name__,
            )
        return context

    def emit(self, document: SchemaDocument) -> StructType:
        """Emit a table schema, requiring the result to be a struct."""
        validator_cls = document.jsonschema_validator_cls or Draft202012Validator
        result = self.emit_node(document.root, self.build_context(validator_cls))
        if isinstance(result, StructType):
            return result
        msg = (
            "Root JSON Schema did not resolve to a StructType.\nThis can happen if it has no "
            "'properties' but does have 'patternProperties' or 'additionalProperties', which map to "
            "MapType instead. PySpark requires a StructType at the top level of the table schema."
        )
        raise PySparkEmitterError(msg)

    def emit_node(self, node: TypedNode, ctx: ConversionContext | None = None, depth: int = 0) -> DataType:
        """Emit a scalar or container node without a table-root restriction."""
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
        scalars = {"boolean": BooleanType(), "null": NullType(), "collapsed-union": StringType()}
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
    ) -> dict[str, object]:
        """Copy selected field metadata into independent mutable containers."""
        values = node.aggregate_metadata(field)
        selected = DEFAULT_METADATA_KEYWORDS | ctx.allowed_extra_metadata_keywords
        jsonschema = {key: mutable_value(values[key]) for key in selected if key in values}
        metadata: dict[str, object] = {}
        if jsonschema:
            metadata["jsonschema"] = jsonschema
        if "description" in values:
            metadata["comment"] = mutable_value(values["description"])
        if self.emit_unions_as_json:
            declaration = node.declared_type
            if isinstance(declaration, tuple):
                non_null = [kind for kind in declaration if kind != "null"]
                if len(non_null) > 1:
                    metadata["original_union"] = list(non_null)
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
                raise PySparkEmitterError(msg)
            if self.emit_unions_as_json:
                logger.info("Emitting multi-type union %r as JSON string.", list(declaration))
            else:
                logger.warning("Collapsing multi-type union %r to StringType.", list(declaration))
            return "collapsed-union"
        if declaration is not None:
            return declaration
        if "enum" in node.constraints:
            return "enum"
        return node.inferred_type or node.type

    def _emit_unknown(self, node: TypedNode, ctx: ConversionContext) -> DataType:
        """Approximate combiners or apply strict unsupported-type handling."""
        if node.type in {"any", "never"} and (node.provenance is None or node.provenance.source == "json-schema"):
            message = (
                f"Boolean JSON Schema `{node.type == 'any'}` has no PySpark equivalent "
                "(unconstrained 'any' type or unsatisfiable 'never' type)."
            )
        else:
            for keyword, branches in (("oneOf", node.one_of), ("anyOf", node.any_of)):
                if branches:
                    logger.warning(
                        "Approximating '%s' by using only its first branch for type inference; PySpark has no union type.",
                        keyword,
                    )
                    return self.emit_node(branches[0], ctx)
            for keyword in ("not", "if", "then", "else"):
                if keyword in node.schema_keywords:
                    logger.warning(
                        "Ignoring unsupported conditional keyword '%s'; it has no effect on the resulting PySpark type.",
                        keyword,
                    )
            message = f"Unsupported/unknown schema type: {node.declared_type or node.type!r}"
        if not self.treat_unknown_as_string:
            raise PySparkEmitterError(message)
        logger.warning("%s Falling back to StringType.", message)
        return StringType()

    def _emit_object(self, node: TypedNode, ctx: ConversionContext, depth: int = 0) -> DataType:
        """Emit fixed fields, a dynamic map, or an empty struct."""
        if depth >= self.max_nesting:
            logger.warning("Max nesting limit (%d) reached; collapsing object to StringType.", self.max_nesting)
            return StringType()

        additional = node.additional_properties
        schema_additional = additional is not None and additional.type not in {"any", "never"}
        if node.properties:
            if node.pattern_properties or schema_additional:
                logger.warning(
                    "Object schema declares both fixed 'properties' and dynamic 'patternProperties' / "
                    "'additionalProperties'; PySpark's StructType requires fixed field names, so the "
                    "dynamic keys are ignored."
                )
            return StructType(
                [
                    StructField(
                        prop.name,
                        self.emit_node(prop.node, ctx, depth + 1),
                        nullable=not prop.required,
                        metadata=self.build_metadata(prop.node, ctx, prop),
                    )
                    for prop in node.properties
                ]
            )
        value_node = additional if schema_additional else None
        if node.pattern_properties:
            if len(node.pattern_properties) > 1:
                logger.warning(
                    "Multiple 'patternProperties' patterns found; only the first pattern's schema is "
                    "used as the MapType value type, since PySpark's MapType has a single value type."
                )
            value_node = next(iter(node.pattern_properties.values()))
        if value_node is not None:
            return MapType(StringType(), self.emit_node(value_node, ctx, depth + 1), valueContainsNull=True)
        logger.warning(
            "Object schema has no 'properties', 'patternProperties', or schema-valued "
            "'additionalProperties'; converting to an empty StructType()."
        )
        return StructType([])

    def _emit_map(self, node: TypedNode, ctx: ConversionContext, depth: int = 0) -> MapType:
        """Emit a map with independently typed keys and values."""
        if node.key_type is None or node.value_type is None:
            msg = "Map nodes require both key_type and value_type"
            raise PySparkEmitterError(msg)
        return MapType(
            self.emit_node(node.key_type, ctx, depth + 1),
            self.emit_node(node.value_type, ctx, depth + 1),
            valueContainsNull=True,
        )

    def _emit_array(self, node: TypedNode, ctx: ConversionContext, depth: int = 0) -> ArrayType:
        """Emit homogeneous arrays, approximating tuple schemas by their first slot."""
        element = node.items
        if node.prefix_items is not None:
            logger.warning(
                "Array schema uses 'prefixItems' (tuple validation); PySpark arrays are homogeneous, so "
                "only the first tuple slot's type is used as the array's element type."
            )
            element = node.prefix_items[0] if node.prefix_items else None
        elif isinstance(element, tuple):
            logger.warning(
                "Array schema uses tuple-style 'items' (list of schemas); PySpark arrays are "
                "homogeneous, so only the first item's type is used as the array's element type."
            )
            element = element[0] if element else None
        if element is None or is_unconstrained_node(element):
            element_type = StringType()
        else:
            element_type = self.emit_node(element, ctx, depth + 1)
        return ArrayType(element_type, containsNull=True)

    def _emit_string(self, node: TypedNode) -> DataType:
        """Map recognized formats to Spark types."""
        value = node.constraints.get("format")
        return self.format_map.get(value, StringType()) if isinstance(value, str) else StringType()

    @staticmethod
    def _emit_integer(node: TypedNode) -> DataType:
        """Use int32 only when both bounds fit."""
        minimum = node.constraints.get("minimum", node.constraints.get("exclusiveMinimum"))
        maximum = node.constraints.get("maximum", node.constraints.get("exclusiveMaximum"))
        if (
            isinstance(minimum, (int, float, Decimal))
            and isinstance(maximum, (int, float, Decimal))
            and minimum >= INT32_MIN
            and maximum <= INT32_MAX
        ):
            return IntegerType()
        return LongType()

    @staticmethod
    def _emit_number(node: TypedNode) -> DataType:
        """Derive decimal scale from numeric JSON multipleOf values."""
        multiple = node.constraints.get("multipleOf")
        if isinstance(multiple, (int, float, Decimal)):
            scale = decimal_places(multiple)
            if scale > MAX_DECIMAL_PRECISION:
                msg = f"Unsupported Spark decimal precision/scale: ({MAX_DECIMAL_PRECISION}, {scale})"
                raise PySparkEmitterError(msg)
            return DecimalType(MAX_DECIMAL_PRECISION, scale)
        return DoubleType()

    @staticmethod
    def _emit_logical(node: TypedNode) -> DataType | None:
        """Preserve physical widths and logical source types without a JSON encoding."""
        if node.provenance is not None and node.provenance.source == "json-schema":
            return None
        hints = node.hints
        logical = hints.logical_type
        if logical in {"decimal", "wei"}:
            precision = hints.precision if hints.precision is not None else MAX_DECIMAL_PRECISION
            scale = hints.scale if hints.scale is not None else 0
            if not 1 <= precision <= MAX_DECIMAL_PRECISION or scale > precision:
                msg = f"Unsupported Spark decimal precision/scale: ({precision}, {scale})"
                raise PySparkEmitterError(msg)
            return DecimalType(precision, scale)
        if logical in {"timestamp", "timestamp-with-tz", "timestamp-without-tz"}:
            return (
                TimestampNTZType() if hints.timezone is False or logical == "timestamp-without-tz" else TimestampType()
            )
        if node.type == "integer" and hints.bit_width is not None:
            return IntegerType() if hints.bit_width <= INT32_WIDTH else LongType()
        if node.type == "number" and (logical == "float" or hints.bit_width == INT32_WIDTH):
            return FloatType()
        return {"date": DateType(), "binary": BinaryType(), "fixed": BinaryType()}.get(logical)


def infer_type_from_enum(values: list[object]) -> DataType:
    """Map homogeneous JSON enum values to Spark scalar types."""
    return {
        "boolean": BooleanType(),
        "integer": LongType(),
        "number": DoubleType(),
        "string": StringType(),
    }[json_type_from_enum(values)]
