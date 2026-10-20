"""Read dereferenced JSON Schema without applying target conversion policy."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any, Final, cast

from jsonschema.validators import validator_for

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry, Extensions
from cdm_data_loaders.converters.core.guards import reject_unresolved_references, require_schema_keyword
from cdm_data_loaders.converters.core.inference import infer_implicit_type
from cdm_data_loaders.converters.core.ir import Field, NodeType, Provenance, SchemaDocument, TypedNode

TYPES: Final = frozenset({"object", "array", "string", "integer", "number", "boolean", "null"})
ANNOTATIONS: Final = frozenset(
    {
        "$schema",
        "$id",
        "$anchor",
        "$dynamicAnchor",
        "$comment",
        "$vocabulary",
        "title",
        "description",
        "default",
        "examples",
        "readOnly",
        "writeOnly",
        "deprecated",
    }
)
SCHEMA_KEYS: Final = frozenset(
    {
        "additionalItems",
        "unevaluatedItems",
        "contains",
        "unevaluatedProperties",
        "propertyNames",
        "not",
        "if",
        "then",
        "else",
        "contentSchema",
    }
)
SCHEMA_MAPS: Final = frozenset({"$defs", "definitions", "dependentSchemas"})
STRUCTURE: Final = (
    SCHEMA_KEYS
    | SCHEMA_MAPS
    | {
        "type",
        "properties",
        "required",
        "items",
        "prefixItems",
        "patternProperties",
        "additionalProperties",
        "anyOf",
        "oneOf",
        "dependencies",
    }
)


@dataclass(frozen=True, slots=True)
class JsonSchemaReader:
    """Read validated object-compatible documents using an explicit extension registry."""

    extension_registry: ExtensionRegistry = DEFAULT_EXTENSIONS
    error_type: type[ConversionError] = ConversionError
    converter_name: str = "JsonSchemaReader"
    dereference_function: str = "dereference_schema"

    def read(self, source: Mapping[str, Any]) -> SchemaDocument:
        """Read a document; validation and reference resolution are upstream steps."""
        if not isinstance(source, Mapping):
            msg = "A JSON Schema document must be an object"
            raise self.error_type(msg)
        schema = dict(source)
        require_schema_keyword(schema, self.error_type, converter_name=self.converter_name)
        declaration = schema.get("type")
        if (
            declaration is not None
            and declaration != "object"
            and not (isinstance(declaration, list) and "object" in declaration)
        ):
            msg = "Root schema must be object-compatible"
            raise self.error_type(msg)
        node = self._node(schema, ())
        return SchemaDocument(
            root=node,
            dialect=schema["$schema"],
            identifier=schema.get("$id"),
            annotations=node.annotations,
            extensions=node.extensions,
            provenance=Provenance("json-schema"),
            jsonschema_validator_cls=validator_for(schema),
        )

    def read_node(self, source: Mapping[str, Any] | bool) -> TypedNode:  # noqa: FBT001
        """Read a schema fragment without document-level dialect or root guards."""
        return self._node(source, ())

    def _node(self, schema: object, path: tuple[str | int, ...]) -> TypedNode:
        """Traverse only schema-valued keywords, leaving literal data untouched."""
        provenance = Provenance("json-schema", path)
        if isinstance(schema, bool):
            return TypedNode(type="any" if schema else "never", nullable=schema, provenance=provenance)
        if not isinstance(schema, Mapping):
            msg = f"Expected a schema object or boolean at {path}"
            raise self.error_type(msg)
        reject_unresolved_references(
            dict(schema),
            self.error_type,
            converter_name=self.converter_name,
            dereference_function=self.dereference_function,
        )
        if "$dynamicRef" in schema or "$recursiveRef" in schema:
            msg = f"Unresolved dynamic reference at {path}"
            raise self.error_type(msg)
        declaration = schema.get("type")
        inferred = infer_implicit_type(dict(schema))
        kind = declaration if isinstance(declaration, str) and declaration in TYPES else "unknown"
        nullable = None
        if declaration is not None:
            nullable = "null" in declaration if isinstance(declaration, list) else declaration == "null"
        constraints = {
            key: value
            for key, value in schema.items()
            if key not in STRUCTURE and key not in ANNOTATIONS and not key.startswith("x-")
        }
        dependencies = schema.get("dependencies", {})
        literal_dependencies = {key: value for key, value in dependencies.items() if isinstance(value, list)}
        if literal_dependencies:
            constraints["dependencies"] = literal_dependencies
        schema_maps = {
            key: {name: self._node(child, (*path, key, name)) for name, child in schema[key].items()}
            for key in SCHEMA_MAPS
            if key in schema
        }
        if "dependencies" in schema:
            schema_maps["dependencies"] = {
                name: self._node(child, (*path, "dependencies", name))
                for name, child in dependencies.items()
                if not isinstance(child, list)
            }
        required = schema.get("required", [])
        items = schema.get("items")
        return TypedNode(
            type=cast("NodeType", kind),
            nullable=nullable,
            declared_type=tuple(declaration) if isinstance(declaration, list) else declaration,
            inferred_type=inferred,
            properties=tuple(
                Field(name, self._node(child, (*path, "properties", name)), required=name in required)
                for name, child in schema.get("properties", {}).items()
            ),
            required_names=tuple(required) if "required" in schema else None,
            items=(
                tuple(self._node(child, (*path, "items", index)) for index, child in enumerate(items))
                if isinstance(items, list)
                else self._node(items, (*path, "items"))
            )
            if "items" in schema
            else None,
            prefix_items=self._branches(schema, "prefixItems", path),
            pattern_properties={
                name: self._node(child, (*path, "patternProperties", name))
                for name, child in schema.get("patternProperties", {}).items()
            },
            additional_properties=self._node(schema["additionalProperties"], (*path, "additionalProperties"))
            if "additionalProperties" in schema
            else None,
            any_of=self._branches(schema, "anyOf", path),
            one_of=self._branches(schema, "oneOf", path),
            schema_keywords={key: self._node(schema[key], (*path, key)) for key in SCHEMA_KEYS if key in schema},
            schema_maps=schema_maps,
            constraints=constraints,
            annotations={key: value for key, value in schema.items() if key in ANNOTATIONS},
            extensions=Extensions(
                {key: value for key, value in schema.items() if key.startswith("x-")},
                self.extension_registry,
            ),
            provenance=provenance,
            source_keywords=tuple(schema),
        )

    def _branches(
        self,
        schema: Mapping[str, Any],
        keyword: str,
        path: tuple[str | int, ...],
    ) -> tuple[TypedNode, ...] | None:
        """Retain branch order and distinguish an absent keyword from an empty list."""
        if keyword not in schema:
            return None
        return tuple(self._node(child, (*path, keyword, index)) for index, child in enumerate(schema[keyword]))
