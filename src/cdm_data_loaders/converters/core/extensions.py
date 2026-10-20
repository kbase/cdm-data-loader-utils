"""Explicit, immutable extension namespaces validated using JSON Schema."""

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from decimal import Decimal
from fractions import Fraction
from typing import Any, Final

from frozendict import frozendict
from jsonschema import Draft202012Validator
from jsonschema.exceptions import ValidationError

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.ir_values import Value, freeze_mapping, freeze_value, mutable_value
from cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.custom_metaschema import X_XSV_CONFIG_SCHEMA


class ExtensionError(ConversionError):
    """An extension namespace or payload violates its registered contract."""


def _validation_value(value: Value) -> object:
    """Copy validation inputs with exact arithmetic and unchanged numeric types."""
    if isinstance(value, Mapping):
        return {key: _validation_value(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_validation_value(item) for item in value]
    if isinstance(value, float):
        number = Fraction(str(value))
        return number.numerator if number.denominator == 1 else number
    if isinstance(value, Decimal):
        return Fraction(value)
    return value


@dataclass(frozen=True, slots=True)
class ExtensionSpec:
    """An exact namespace spelling and its complete validation schema."""

    namespace: str
    schema: Mapping[str, Value] | bool

    def __post_init__(self) -> None:
        """Validate the namespace and own a valid draft 2020-12 schema."""
        if not isinstance(self.namespace, str) or not self.namespace.startswith("x-") or self.namespace == "x-":
            msg = "Extension namespaces must start with 'x-' and have a name"
            raise ExtensionError(msg)
        if not isinstance(self.schema, (Mapping, bool)):
            msg = "Extension validation schemas must be mappings or booleans"
            raise ExtensionError(msg)
        owned = freeze_value(self.schema)
        Draft202012Validator.check_schema(mutable_value(owned))
        object.__setattr__(self, "schema", owned)


@dataclass(frozen=True, slots=True)
class ExtensionRegistry:
    """A closed set of namespace contracts; additions produce a new registry."""

    specs: tuple[ExtensionSpec, ...] = ()

    def __post_init__(self) -> None:
        """Reject duplicate contracts and isolate the declaration sequence."""
        object.__setattr__(self, "specs", tuple(self.specs))
        if any(not isinstance(spec, ExtensionSpec) for spec in self.specs):
            msg = "Registry entries must be ExtensionSpec instances"
            raise TypeError(msg)
        if len({spec.namespace for spec in self.specs}) != len(self.specs):
            msg = "Duplicate extension namespace"
            raise ExtensionError(msg)

    def register(self, spec: ExtensionSpec, *, replace: bool = False) -> "ExtensionRegistry":
        """Add a contract, requiring an explicit flag to replace an existing one."""
        retained = tuple(existing for existing in self.specs if existing.namespace != spec.namespace)
        if len(retained) != len(self.specs) and not replace:
            msg = f"Already registered: {spec.namespace}"
            raise ExtensionError(msg)
        return ExtensionRegistry((*retained, spec))

    def validate(self, values: Mapping[str, object]) -> Mapping[str, Value]:
        """Validate every namespace and return an isolated immutable payload."""
        owned = freeze_mapping(values)
        specs = {spec.namespace: spec for spec in self.specs}
        for namespace, value in owned.items():
            if namespace not in specs:
                msg = f"Unregistered extension namespace: {namespace}"
                raise ExtensionError(msg)
            try:
                Draft202012Validator(_validation_value(specs[namespace].schema)).validate(_validation_value(value))
            except ValidationError as error:
                msg = f"Invalid {namespace} at {list(error.path)}: {error.message}"
                raise ExtensionError(msg) from error
        return owned

    def extend_payload(self, namespace: str, properties: Mapping[str, object]) -> "ExtensionRegistry":
        """Register additional exact payload keys on a closed object namespace."""
        for spec in self.specs:
            if spec.namespace == namespace:
                schema = mutable_value(spec.schema)
                if not isinstance(schema, dict) or schema.get("type") != "object":
                    msg = f"Not an object namespace: {namespace}"
                    raise ExtensionError(msg)
                existing = schema.setdefault("properties", {})
                if existing.keys() & properties.keys():
                    msg = "Payload keys are already registered"
                    raise ExtensionError(msg)
                existing.update(properties)
                return self.register(ExtensionSpec(namespace, schema), replace=True)
        msg = f"Unregistered extension namespace: {namespace}"
        raise ExtensionError(msg)


def _object_schema(properties: Mapping[str, object]) -> dict[str, Any]:
    """Build a closed namespace payload schema."""
    return {"type": "object", "properties": dict(properties), "additionalProperties": False}


_DLT_PROPERTIES: Final = {
    **{
        name: {"type": ["boolean", "null"]}
        for name in (
            "nullable",
            "primary_key",
            "unique",
            "foreign_key",
            "sort",
            "cluster",
            "partition",
            "merge_key",
            "row_key",
            "root_key",
            "variant",
            "timezone",
        )
    },
    **{name: {"type": ["integer", "null"], "minimum": 0} for name in ("precision", "scale")},
    "name": {"type": "string"},
    "description": {"type": ["string", "null"]},
    "data_type": {"type": ["string", "null"]},
}
_ICEBERG_PROPERTIES: Final = {
    **{
        name: {"type": "integer"}
        for name in (
            "field_id",
            "schema_id",
            "format_version",
            "element_id",
            "key_id",
            "value_id",
        )
    },
    **{name: {"type": "integer", "minimum": 0} for name in ("precision", "scale", "length")},
    **{name: {"type": "boolean"} for name in ("required", "element_required", "value_required")},
    **{name: {"type": "string"} for name in ("logical_type", "location", "generated_at")},
    "initial_default": True,
    "write_default": True,
    "identifier": {"type": "array", "items": {"type": "string"}},
    "identifier_field_ids": {"type": "array", "items": {"type": "integer"}},
    "current_snapshot_id": {"type": ["integer", "null"]},
    "properties": {"type": "object", "additionalProperties": {"type": "string"}},
    "partition_spec": {
        "type": "array",
        "items": {
            **_object_schema(
                {"source_id": {"type": "integer"}, "name": {"type": "string"}, "transform": {"type": "string"}}
            ),
            "required": ["source_id", "name", "transform"],
        },
    },
}
DEFAULT_EXTENSIONS: Final = ExtensionRegistry(
    (
        ExtensionSpec("x-dlt", _object_schema(_DLT_PROPERTIES)),
        ExtensionSpec("x-iceberg", _object_schema(_ICEBERG_PROPERTIES)),
        ExtensionSpec("x-xsv-config", X_XSV_CONFIG_SCHEMA),
        ExtensionSpec("x-file-glob", {"type": "string"}),
        ExtensionSpec("x-delimiter", {"type": "string", "minLength": 1, "maxLength": 1}),
        ExtensionSpec("x-dlt-prefix", {"type": "string"}),
        ExtensionSpec("x-dlt-split", {"type": "string"}),
        ExtensionSpec("x-pii", {"type": "boolean"}),
    )
)


@dataclass(frozen=True, slots=True)
class Extensions(Mapping[str, Value]):
    """A validated immutable mapping that retains its validation contract."""

    payload: Mapping[str, Value] = field(default_factory=frozendict)
    registry: ExtensionRegistry = field(default=DEFAULT_EXTENSIONS, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Validate normal construction and own all nested payloads."""
        object.__setattr__(self, "payload", self.registry.validate(self.payload))

    def __getitem__(self, key: str) -> Value:
        """Return a namespace payload."""
        return self.payload[key]

    def __iter__(self) -> Iterator[str]:
        """Iterate namespace names in declaration order."""
        return iter(self.payload)

    def __len__(self) -> int:
        """Return the number of namespaces."""
        return len(self.payload)


def validated_extensions(values: Mapping[str, Value]) -> Extensions:
    """Retain an already validated custom contract or apply the default registry."""
    return values if isinstance(values, Extensions) else Extensions(values)
