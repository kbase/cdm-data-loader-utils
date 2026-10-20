"""Convert dlt schemas through the typed reader and JSON Schema emitter."""

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, InstanceOf

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry
from cdm_data_loaders.converters.core.io import load_schema_file, load_schema_text
from cdm_data_loaders.converters.emitters.json_schema import JsonSchemaEmitter
from cdm_data_loaders.converters.readers.dlt import DltReader


class DltToJSONSchemaError(ConversionError):
    """Raised when a dlt schema construct cannot be converted."""


class DltToJSONSchema(BaseModel):
    """Read root table documents and emit their JSON Schema equivalents."""

    model_config = ConfigDict(frozen=True, arbitrary_types_allowed=True)

    child_table_mode: Literal["object", "array"] = "object"
    include_dlt_columns: bool = False
    include_variant_columns: bool = False
    preserve_unknown_hints: bool = True
    extension_registry: InstanceOf[ExtensionRegistry] = DEFAULT_EXTENSIONS

    def convert(self, schema: dict[str, Any]) -> dict[str, dict[str, Any]]:
        """Convert a stored schema or table mapping into root JSON documents."""
        reader = DltReader(
            include_dlt_columns=self.include_dlt_columns,
            include_variant_columns=self.include_variant_columns,
            child_table_mode=self.child_table_mode,
            extension_registry=self.extension_registry,
        )
        emitter = JsonSchemaEmitter(preserve_unknown_hints=self.preserve_unknown_hints)
        try:
            return {name: emitter.emit(document) for name, document in reader.read(schema).items()}
        except ConversionError as error:
            raise DltToJSONSchemaError(str(error)) from error

    def convert_from_string(self, schema_str: str) -> dict[str, dict[str, Any]]:
        """Convert a JSON or YAML string."""
        return self.convert(load_schema_text(schema_str, allow_yaml=True))

    def convert_from_file(self, path: str) -> dict[str, dict[str, Any]]:
        """Convert a JSON or YAML file."""
        return self.convert(load_schema_file(path))
