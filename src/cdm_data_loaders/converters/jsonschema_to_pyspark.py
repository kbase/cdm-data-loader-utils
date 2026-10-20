"""Convert JSON Schema through the shared reader and PySpark emitter."""

from typing import Any, cast

from frozendict import frozendict
from pydantic import BaseModel, ConfigDict, Field, InstanceOf, field_validator
from pyspark.sql.types import ArrayType, DataType, StructType

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry
from cdm_data_loaders.converters.core.guards import (
    reject_unresolved_references,
    require_object_root,
    require_schema_keyword,
)
from cdm_data_loaders.converters.core.io import load_schema_file, load_schema_text
from cdm_data_loaders.converters.emitters.pyspark import (
    DEFAULT_FORMAT_MAP,
    ConversionContext,
    PySparkEmitter,
    merge_format_map,
)
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader


class JSONSchemaToPySparkError(ConversionError):
    """Raised when a JSON Schema construct cannot be converted."""


class InvalidJSONSchemaError(JSONSchemaToPySparkError):
    """Raised when the input document is not a valid JSON Schema."""


class JSONSchemaToPySpark(BaseModel):
    """Read dereferenced JSON Schema and emit pure Python Spark schema objects."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    format_map: frozendict[str, DataType] = Field(default_factory=lambda: DEFAULT_FORMAT_MAP)
    treat_unknown_as_string: bool = True
    extra_metadata_keywords: frozenset[str] = Field(default_factory=frozenset)
    extension_registry: InstanceOf[ExtensionRegistry] = DEFAULT_EXTENSIONS

    @field_validator("format_map", mode="before")
    @classmethod
    def _merge_with_builtin_format_map(cls, value: object) -> object:
        """Delegate nonempty format overrides to the emitter."""
        return merge_format_map(value) if value else value

    @field_validator("extra_metadata_keywords", mode="before")
    @classmethod
    def _coerce_extra_metadata_keywords(cls, value: object) -> object:
        """Normalize supported metadata keyword collections."""
        return frozenset(value) if isinstance(value, (set, list, tuple)) else value

    def _reader(self) -> JsonSchemaReader:
        """Configure reader diagnostics for this public facade."""
        return JsonSchemaReader(
            self.extension_registry,
            error_type=JSONSchemaToPySparkError,
            converter_name="JSONSchemaToPySpark",
            dereference_function="jsonschema_to_pyspark.dereferencing.dereference_schema",
        )

    def _emitter(self) -> PySparkEmitter:
        """Pass target policy to the emitter without serializing arbitrary types."""
        return PySparkEmitter(
            format_map=self.format_map,
            treat_unknown_as_string=self.treat_unknown_as_string,
            extra_metadata_keywords=self.extra_metadata_keywords,
        )

    def _build_context(self, validator_cls: type) -> ConversionContext:
        """Build the emitter's draft-specific metadata context."""
        return self._emitter().build_context(validator_cls)

    def convert(self, schema: dict[str, Any]) -> StructType:
        """Convert an object document, preserving direction-specific guards."""
        require_schema_keyword(schema, InvalidJSONSchemaError, converter_name="JSONSchemaToPySpark")
        require_object_root(schema, JSONSchemaToPySparkError, target="a StructType")
        try:
            return self._emitter().emit(self._reader().read(schema))
        except ConversionError as error:
            raise JSONSchemaToPySparkError(str(error)) from error

    def convert_from_string(self, schema_str: str) -> StructType:
        """Convert a JSON string."""
        return self.convert(load_schema_text(schema_str))

    def convert_from_file(self, path: str) -> StructType:
        """Convert a JSON or YAML file."""
        return self.convert(load_schema_file(path))

    @staticmethod
    def _reject_unresolved_references(schema: dict[str, Any]) -> None:
        """Reject unresolved references with the public converter's diagnostic."""
        reject_unresolved_references(
            schema,
            JSONSchemaToPySparkError,
            converter_name="JSONSchemaToPySpark",
            dereference_function="jsonschema_to_pyspark.dereferencing.dereference_schema",
        )

    def _build_metadata(self, schema: dict[str, Any], ctx: ConversionContext) -> dict[str, Any]:
        """Delegate metadata selection to the emitter."""
        return self._emitter().build_metadata(self._reader().read_node(schema), ctx)

    def _convert_type(self, schema: dict[str, Any] | bool, ctx: ConversionContext) -> DataType:  # noqa: FBT001
        """Read and emit a fragment without document-level root guards."""
        try:
            return self._emitter().emit_node(self._reader().read_node(schema), ctx)
        except ConversionError as error:
            raise JSONSchemaToPySparkError(str(error)) from error

    def _convert_property(
        self,
        prop_schema: dict[str, Any] | bool,  # noqa: FBT001
        ctx: ConversionContext,
    ) -> tuple[DataType, dict[str, Any]]:
        """Delegate property type and metadata emission from one reader node."""
        try:
            node = self._reader().read_node(prop_schema)
            emitter = self._emitter()
            return emitter.emit_node(node, ctx), emitter.build_metadata(node, ctx)
        except ConversionError as error:
            raise JSONSchemaToPySparkError(str(error)) from error

    def _convert_boolean_schema(self, schema: bool) -> DataType:  # noqa: FBT001
        """Emit a boolean schema using the configured unknown-type policy."""
        try:
            return self._emitter().emit_node(self._reader().read_node(schema))
        except ConversionError as error:
            raise JSONSchemaToPySparkError(str(error)) from error

    def _dispatch_type(self, schema: dict[str, Any], ctx: ConversionContext) -> DataType:
        """Delegate fragment dispatch to the emitter."""
        return self._convert_type(schema, ctx)

    def _convert_object(self, schema: dict[str, Any], ctx: ConversionContext) -> DataType:
        """Emit an object fragment."""
        return self._convert_type({**schema, "type": "object"}, ctx)

    def _convert_array(self, schema: dict[str, Any], ctx: ConversionContext) -> ArrayType:
        """Emit an array fragment."""
        return cast("ArrayType", self._convert_type({**schema, "type": "array"}, ctx))

    def _convert_string(self, schema: dict[str, Any]) -> DataType:
        """Emit a string fragment using the configured format map."""
        return self._emitter().emit_node(self._reader().read_node({**schema, "type": "string"}))

    @staticmethod
    def _convert_integer(schema: dict[str, Any]) -> DataType:
        """Emit integer bounds through the shared reader and emitter."""
        return PySparkEmitter().emit_node(JsonSchemaReader().read_node({**schema, "type": "integer"}))

    @staticmethod
    def _convert_number(schema: dict[str, Any]) -> DataType:
        """Emit number precision through the shared reader and emitter."""
        return PySparkEmitter().emit_node(JsonSchemaReader().read_node({**schema, "type": "number"}))
