"""Emit JSON Schema from typed schema facts and scoped metadata."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Final

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.inference import decimal_pattern
from cdm_data_loaders.converters.core.ir import Field, SchemaDocument, TypedNode
from cdm_data_loaders.converters.core.ir_values import Value, mutable_value

JSON_SCHEMA_DIALECT: Final = "https://json-schema.org/draft/2020-12/schema"
DLT_TYPES: Final[dict[str, dict[str, Any]]] = {
    "text": {"type": "string"},
    "bigint": {"type": "integer"},
    "double": {"type": "number"},
    "decimal": {"type": "string", "pattern": r"^-?\d+(\.\d+)?$", "x-dlt": {"data_type": "decimal"}},
    "bool": {"type": "boolean"},
    "timestamp": {"type": "string", "format": "date-time", "x-dlt": {"data_type": "timestamp"}},
    "date": {"type": "string", "format": "date"},
    "time": {"type": "string", "format": "time"},
    "binary": {"type": "string", "contentEncoding": "base64"},
    "json": {},
    "wei": {"type": "integer", "x-dlt": {"data_type": "wei"}},
}
DLT_KNOWN: Final = frozenset({"name", "data_type", "nullable", "description", "precision", "scale"})
ICEBERG_LOGICAL_TYPES: Final = frozenset(
    {
        "boolean",
        "int",
        "long",
        "float",
        "double",
        "string",
        "date",
        "time",
        "uuid",
        "timestamp-without-tz",
        "timestamp-with-tz",
        "decimal",
        "fixed",
    }
)
ICEBERG_ELEMENT_KEYS: Final = frozenset({"element_id", "element_required", "key_id", "value_id", "value_required"})


def _values(values: Mapping[str, Value]) -> dict[str, Any]:
    """Copy literal containers without changing numeric values."""
    return {key: mutable_value(value) for key, value in values.items()}


def _nullable(schema: dict[str, Any], *, nullable: bool) -> dict[str, Any]:
    """Wrap a value schema without moving its metadata."""
    return {"anyOf": [schema, {"type": "null"}]} if nullable else schema


def _schema_object(*, schema: dict[str, Any] | bool) -> dict[str, Any]:
    """Represent a boolean schema as an object when metadata must be attached."""
    if isinstance(schema, bool):
        return {} if schema else {"not": {}}
    return schema


@dataclass(frozen=True, slots=True)
class JsonSchemaEmitter:
    """Render typed trees using the source's JSON Schema compatibility policy."""

    preserve_unknown_hints: bool = True

    def __post_init__(self) -> None:
        """Reject nonboolean hint policies."""
        if type(self.preserve_unknown_hints) is not bool:
            msg = "preserve_unknown_hints must be a boolean"
            raise TypeError(msg)

    def emit(self, document: SchemaDocument) -> dict[str, Any]:
        """Return an independently owned JSON Schema document."""
        provenance = document.provenance or document.root.provenance
        source = provenance.source if provenance is not None else None
        if source == "dlt":
            result = self._dlt_node(document.root)
        elif source == "iceberg":
            result = self._iceberg_node(document.root)
        else:
            return self._json_document(document)
        annotations = _values(document.annotations)
        if source == "iceberg" and not annotations.get("description"):
            annotations.pop("description", None)
        result.update(annotations)
        result.update(_values(document.extensions))
        result["$schema"] = document.dialect or JSON_SCHEMA_DIALECT
        if document.identifier is not None:
            result["$id"] = document.identifier
        elif source == "dlt":
            result["$id"] = f"urn:dlt:{document.name}"
        if document.name is not None:
            result["title"] = document.name
        return result

    def emit_into(self, document: SchemaDocument, schema: dict[str, Any]) -> dict[str, Any]:
        """Add a document's schema to a namespace document's ``$defs`` without overwriting existing entries.

        :param document: Named, table-rooted document to convert.
        :param schema: Destination namespace document, updated only after successful conversion.
        :returns: The destination namespace document.
        """
        if document.name is None:
            msg = "A namespaced JSON Schema entry requires a named document"
            raise ConversionError(msg)
        table_schema = self.emit(document)
        table_schema.pop("$schema", None)
        schema.setdefault("$defs", {})[document.name] = table_schema
        return schema

    def emit_node(self, node: TypedNode) -> dict[str, Any] | bool:
        """Emit a standalone node using its source's rendering policy."""
        source = node.provenance.source if node.provenance is not None else None
        if source == "iceberg":
            return self._iceberg_node(node)
        if source == "dlt":
            return self._dlt_node(node)
        return self._json_node(node)

    def emit_field(self, field: Field) -> dict[str, Any] | bool:
        """Emit a field, retaining source-specific metadata placement."""
        source = field.node.provenance.source if field.node.provenance is not None else None
        if source == "iceberg":
            return self._iceberg_field(field)
        if source == "dlt":
            return self._dlt_node(field.node)
        return self._json_field(field)

    def _json_document(self, document: SchemaDocument) -> dict[str, Any]:
        """Emit structural JSON facts without changing the source dialect."""
        result = _schema_object(schema=self._json_node(document.root))
        result.update(_values(document.annotations))
        result.update(_values(document.extensions))
        if document.dialect is not None:
            result["$schema"] = document.dialect
        if document.identifier is not None:
            result["$id"] = document.identifier
        if document.name is not None:
            result.setdefault("title", document.name)
        return result

    def _json_node(self, node: TypedNode) -> dict[str, Any] | bool:
        """Preserve declarations separately from inferred types and literal values."""
        from_json = node.provenance is not None and node.provenance.source == "json-schema"
        result = _values(node.constraints)
        result.update(_values(node.annotations))
        result.update(_values(node.extensions))
        if node.declared_type is not None:
            result["type"] = list(node.declared_type) if isinstance(node.declared_type, tuple) else node.declared_type
        elif not from_json and node.type not in {"any", "never", "unknown", "map"}:
            result["type"] = node.type
        elif not from_json and node.type == "unknown":
            msg = "Cannot emit an unknown type without a JSON Schema declaration"
            raise ConversionError(msg)
        self._json_children(node, result)
        if node.type == "never":
            result["not"] = {}
        if from_json:
            if not node.source_keywords and (
                (node.type == "any" and not result) or (node.type == "never" and result == {"not": {}})
            ):
                return node.type == "any"
            result = {**{key: result[key] for key in node.source_keywords if key in result}, **result}
        else:
            result = _nullable(result, nullable=node.nullable is True and node.type not in {"any", "never", "null"})
        return result

    def _json_children(self, node: TypedNode, result: dict[str, Any]) -> None:
        """Render schema positions without interpreting literal keyword payloads."""
        if node.properties or "properties" in node.source_keywords:
            result["properties"] = {field.name: self._json_field(field) for field in node.properties}
        if node.required_names is not None:
            result["required"] = list(node.required_names)
        elif any(field.required for field in node.properties):
            result["required"] = [field.name for field in node.properties if field.required]
        if node.pattern_properties or "patternProperties" in node.source_keywords:
            result["patternProperties"] = {
                key: self._json_node(child) for key, child in node.pattern_properties.items()
            }
        if node.items is not None:
            result["items"] = (
                [self._json_node(child) for child in node.items]
                if isinstance(node.items, tuple)
                else self._json_node(node.items)
            )
        if node.additional_properties is not None:
            result["additionalProperties"] = self._json_node(node.additional_properties)
        for keyword, children in (("prefixItems", node.prefix_items), ("anyOf", node.any_of), ("oneOf", node.one_of)):
            if children is not None:
                result[keyword] = [self._json_node(child) for child in children]
        result.update({key: self._json_node(child) for key, child in node.schema_keywords.items()})
        for keyword, children in node.schema_maps.items():
            result.setdefault(keyword, {}).update({key: self._json_node(child) for key, child in children.items()})
        if node.type == "map":
            result.update(self._json_map(node))

    def _json_field(self, field: Field) -> dict[str, Any] | bool:
        """Attach field metadata without changing an unannotated boolean schema."""
        schema = self._json_node(field.node)
        if not field.annotations and not field.extensions:
            return schema
        result = _schema_object(schema=schema)
        result.update(_values(field.annotations))
        result.update(_values(field.extensions))
        return result

    def _json_map(self, node: TypedNode) -> dict[str, Any]:
        """Represent arbitrary typed keys as key/value pairs."""
        if node.key_type is None or node.value_type is None:
            msg = "Map nodes require key and value types"
            raise ConversionError(msg)
        return {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {"key": self._json_node(node.key_type), "value": self._json_node(node.value_type)},
                "required": ["key", "value"],
                "additionalProperties": False,
            },
        }

    def _dlt_node(self, node: TypedNode) -> dict[str, Any]:
        """Render table structure separately from nullable columns."""
        if node.type == "object":
            result: dict[str, Any] = {
                "type": "object",
                "properties": {field.name: self._dlt_node(field.node) for field in node.properties},
                "additionalProperties": False,
            }
            required = node.required_names
            if required is None:
                required = tuple(field.name for field in node.properties if field.required)
            if required:
                result["required"] = list(required)
        elif node.type == "array" and isinstance(node.items, TypedNode):
            result = {"type": "array", "items": self._dlt_node(node.items)}
        else:
            return self._dlt_scalar(node)
        result.update(_values(node.annotations))
        result.update(_values(node.extensions))
        return result

    def _dlt_scalar(self, node: TypedNode) -> dict[str, Any]:
        """Select logical templates and keep column metadata inside null unions."""
        logical = node.hints.logical_type
        if node.type == "unknown" or (logical is not None and logical not in DLT_TYPES):
            path = node.provenance.path if node.provenance is not None else ()
            name = path[-1] if path else "<unknown>"
            msg = f"Column {name!r} has unknown dlt data_type {logical!r}."
            raise ConversionError(msg)
        result = deepcopy(DLT_TYPES[logical]) if logical is not None else {}
        if node.hints.precision is not None or node.hints.scale is not None:
            hints = result.setdefault("x-dlt", {})
            hints["precision"] = node.hints.precision
            if node.hints.scale is not None:
                hints["scale"] = node.hints.scale
        result.update(_values(node.constraints))
        result.update(_values(node.annotations))
        metadata = _values(node.extensions)
        source_hints = metadata.pop("x-dlt", {})
        if self.preserve_unknown_hints:
            hints = {key: value for key, value in source_hints.items() if key not in DLT_KNOWN and value is not None}
            if hints:
                result.setdefault("x-dlt", {}).update(hints)
        result.update(metadata)
        return _nullable(result, nullable=bool(node.nullable and result))

    def _iceberg_node(self, node: TypedNode) -> dict[str, Any]:
        """Wrap nullable types and elements before applying field metadata."""
        result = self._iceberg_type(node)
        result.update(_values(node.constraints))
        result.update(_values(node.annotations))
        extensions = _values(node.extensions)
        metadata = extensions.pop("x-iceberg", {})
        metadata = {key: value for key, value in metadata.items() if key not in ICEBERG_ELEMENT_KEYS}
        if metadata:
            result.setdefault("x-iceberg", {}).update(metadata)
        result.update(extensions)
        return _nullable(result, nullable=node.nullable is True)

    def _iceberg_type(self, node: TypedNode) -> dict[str, Any]:
        """Render structural Iceberg nodes directly from their typed children."""
        if node.type == "object":
            result: dict[str, Any] = {
                "type": "object",
                "properties": {field.name: self._iceberg_field(field) for field in node.properties},
                "additionalProperties": False,
            }
            required = node.required_names
            if required is None:
                required = tuple(field.name for field in node.properties if field.required)
            if required:
                result["required"] = list(required)
        elif node.type == "array" and isinstance(node.items, TypedNode):
            result = {"type": "array", "items": self._iceberg_node(node.items)}
        elif node.type == "map" and node.key_type is not None and node.value_type is not None:
            result = {
                "type": "array",
                "description": "Iceberg map rendered as key/value pairs (JSON object keys must be strings).",
                "items": {
                    "type": "object",
                    "properties": {
                        "key": self._iceberg_node(node.key_type),
                        "value": self._iceberg_node(node.value_type),
                    },
                    "required": ["key", "value"],
                    "additionalProperties": False,
                },
            }
        else:
            result = self._iceberg_scalar(node)
        return result

    def _iceberg_scalar(self, node: TypedNode) -> dict[str, Any]:
        """Apply physical bounds and source logical encodings to scalar nodes."""
        hints = node.hints
        logical = hints.logical_type
        if logical not in ICEBERG_LOGICAL_TYPES:
            msg = f"No JSON Schema mapping for Iceberg type: {logical or node.type}"
            raise ConversionError(msg)
        result: dict[str, Any] = {"type": node.type}
        if node.type == "integer" and hints.bit_width:
            bound = 2 ** (hints.bit_width - 1)
            result.update(minimum=-bound, maximum=bound - 1)
        if logical in {"date", "time", "uuid", "timestamp-with-tz"}:
            result["format"] = "date-time" if logical == "timestamp-with-tz" else logical
        elif logical == "timestamp-without-tz":
            result.update(
                description="ISO 8601 timestamp with no timezone attached.",
                **{"x-iceberg": {"logical_type": logical}},
            )
        elif logical == "decimal":
            if hints.precision is None or hints.scale is None:
                msg = "Iceberg decimals require precision and scale"
                raise ConversionError(msg)
            result = {
                "type": "string",
                "pattern": decimal_pattern(hints.precision, hints.scale),
                "x-iceberg": {"logical_type": logical, "precision": hints.precision, "scale": hints.scale},
            }
        elif logical == "fixed":
            if hints.length is None:
                msg = "Iceberg fixed types require length"
                raise ConversionError(msg)
            result = {
                "type": "string",
                "contentEncoding": "base64",
                "x-iceberg": {"logical_type": logical, "length": hints.length},
            }
        return result

    def _iceberg_field(self, field: Field) -> dict[str, Any]:
        """Place descriptions, IDs and non-null defaults outside optional types."""
        result = self._iceberg_node(field.node)
        annotations = _values(field.annotations)
        if not annotations.get("description"):
            annotations.pop("description", None)
        result.update(annotations)
        extensions = _values(field.extensions)
        metadata = extensions.pop("x-iceberg", {})
        metadata = {
            key: value
            for key, value in metadata.items()
            if key not in {"initial_default", "write_default"} or value is not None
        }
        if metadata:
            result.setdefault("x-iceberg", {}).update(metadata)
        result.update(extensions)
        return result
