"""Convert dereferenced JSON Schema through the shared reader and dlt emitter."""

from typing import TYPE_CHECKING, Any

from dlt.common.data_types.typing import TDataType
from dlt.common.schema.typing import TStoredSchema
from pydantic import BaseModel, ConfigDict, Field, InstanceOf, field_validator

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry
from cdm_data_loaders.converters.core.guards import (
    reject_unresolved_references,
    require_object_root,
    require_schema_keyword,
)
from cdm_data_loaders.converters.core.io import load_schema_file, load_schema_text
from cdm_data_loaders.converters.core.ir import SchemaDocument
from cdm_data_loaders.converters.emitters.dlt import (
    DltEmitter,
    _scalar_type,
)
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader

if TYPE_CHECKING:
    from dlt.common.schema import Schema


# FIXME: dead code
def _data_type_from_enum(values: list[Any]) -> TDataType:
    """Infer an enum's dlt type using the shared node conversion."""
    return _scalar_type(JsonSchemaReader().read_node({"enum": values}))


class JSONSchemaToDltError(ConversionError):
    """Raised when a JSON Schema construct cannot be converted."""


class InvalidJSONSchemaError(JSONSchemaToDltError):
    """Raised when the input document is not a valid JSON Schema."""


class JSONSchemaToDlt(BaseModel):
    """Read a dereferenced JSON Schema and emit dlt tables."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    schema_name: str
    skip_nested_types: bool = False
    max_nesting: int = Field(default=10, ge=1)
    write_disposition: str | None = None
    flatten_scalars: bool = True
    extension_registry: InstanceOf[ExtensionRegistry] = DEFAULT_EXTENSIONS

    @field_validator("schema_name")
    @classmethod
    def _validate_schema_name(cls, value: str) -> str:
        """Reject empty schema names."""
        if not value.strip():
            msg = "schema_name must be a non-empty string"
            raise ValueError(msg)
        return value

    def _emitter(self) -> DltEmitter:
        """Pass target policy to the emitter."""
        return DltEmitter(**self.model_dump(exclude={"extension_registry"}))

    def _read(self, schema: dict[str, Any]) -> SchemaDocument:
        """Preserve facade guards and error identity around the reader."""
        require_schema_keyword(schema, InvalidJSONSchemaError, converter_name="JSONSchemaToDlt")
        require_object_root(schema, JSONSchemaToDltError, target="a dlt table")
        try:
            return JsonSchemaReader(
                self.extension_registry,
                error_type=JSONSchemaToDltError,
                converter_name="JSONSchemaToDlt",
                dereference_function="jsonschema_to_dlt.dereferencing.dereference_schema",
            ).read(schema)
        except ConversionError as error:
            if isinstance(error, JSONSchemaToDltError):
                raise
            raise JSONSchemaToDltError(str(error)) from error

    def convert(self, schema: dict[str, Any]) -> TStoredSchema:
        """Convert an object document into a dlt stored schema."""
        return self._emitter().emit(self._read(schema))

    def convert_from_string(self, schema_str: str) -> TStoredSchema:
        """Convert a JSON string, not YAML."""
        return self.convert(load_schema_text(schema_str))

    def convert_from_file(self, path: str) -> TStoredSchema:
        """Convert a JSON or YAML file selected by its extension."""
        return self.convert(load_schema_file(path))

    def to_yaml(self, schema: dict[str, Any]) -> str:
        """Emit a dlt stored schema as YAML."""
        return self._emitter().to_yaml(self._read(schema))

    def to_schema(self, schema: dict[str, Any]) -> "Schema":
        """Emit a live dlt Schema with parent tables installed first."""
        return self._emitter().to_schema(self._read(schema))

    @staticmethod
    def _reject_unresolved_references(schema: dict[str, Any]) -> None:
        """Reject unresolved references with the direction-specific diagnostic."""
        reject_unresolved_references(
            schema,
            JSONSchemaToDltError,
            converter_name="JSONSchemaToDlt",
            dereference_function="jsonschema_to_dlt.dereferencing.dereference_schema",
        )
