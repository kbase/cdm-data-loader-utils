"""Convert JSON Schema through the shared reader and PyArrow emitter."""

from typing import Any

import pyarrow as pa
from frozendict import frozendict
from pydantic import BaseModel, ConfigDict, Field, InstanceOf, field_validator

from cdm_data_loaders.converters.core.errors import ConversionError
from cdm_data_loaders.converters.core.extensions import DEFAULT_EXTENSIONS, ExtensionRegistry
from cdm_data_loaders.converters.core.guards import require_object_root, require_schema_keyword
from cdm_data_loaders.converters.core.io import load_schema_file, load_schema_text
from cdm_data_loaders.converters.emitters.pyarrow import DEFAULT_FORMAT_MAP, PyArrowEmitter, merge_format_map
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader


class JSONSchemaToPyArrowError(ConversionError):
    """Raised when a JSON Schema construct cannot be converted."""


class InvalidJSONSchemaError(JSONSchemaToPyArrowError):
    """Raised when the input document is not a valid JSON Schema."""


class JSONSchemaToPyArrow(BaseModel):
    """Read dereferenced JSON Schema and emit PyArrow schemas and empty tables."""

    model_config = ConfigDict(arbitrary_types_allowed=True, frozen=True)

    format_map: frozendict[str, pa.DataType] = Field(default_factory=lambda: DEFAULT_FORMAT_MAP)
    treat_unknown_as_string: bool = True
    emit_unions_as_json: bool = False
    extra_metadata_keywords: frozenset[str] = Field(default_factory=frozenset)
    max_nesting: int = Field(default=10, ge=1)
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
            error_type=JSONSchemaToPyArrowError,
            converter_name="JSONSchemaToPyArrow",
            dereference_function="cdm_data_loaders.utils.jsonschema.dereferencer.dereference_schema",
        )

    def _emitter(self) -> PyArrowEmitter:
        """Pass target policy to the emitter."""
        return PyArrowEmitter(
            format_map=self.format_map,
            treat_unknown_as_string=self.treat_unknown_as_string,
            emit_unions_as_json=self.emit_unions_as_json,
            extra_metadata_keywords=self.extra_metadata_keywords,
            max_nesting=self.max_nesting,
        )

    def convert(self, schema: dict[str, Any]) -> pa.Schema:
        """Convert a dereferenced JSON Schema object document to an Arrow schema.

        :param schema: the JSON Schema, which must declare ``$schema`` and have an object root
        :type schema: dict[str, Any]
        :return: the Arrow schema
        :rtype: pa.Schema
        :raises InvalidJSONSchemaError: if ``$schema`` is missing
        :raises JSONSchemaToPyArrowError: if the root is not an object or a construct cannot be converted
        """
        require_schema_keyword(schema, InvalidJSONSchemaError, converter_name="JSONSchemaToPyArrow")
        require_object_root(schema, JSONSchemaToPyArrowError, target="a struct schema")
        try:
            return self._emitter().emit(self._reader().read(schema))
        except ConversionError as error:
            raise JSONSchemaToPyArrowError(str(error)) from error

    def convert_to_table(self, schema: dict[str, Any]) -> pa.Table:
        """Convert a dereferenced JSON Schema object document to an empty Arrow table.

        :param schema: the JSON Schema, which must declare ``$schema`` and have an object root
        :type schema: dict[str, Any]
        :return: a zero-row table with the converted schema
        :rtype: pa.Table
        :raises InvalidJSONSchemaError: if ``$schema`` is missing
        :raises JSONSchemaToPyArrowError: if the root is not an object or a construct cannot be converted
        """
        return self.convert(schema).empty_table()

    def convert_from_string(self, schema_str: str) -> pa.Schema:
        """Convert a JSON string to an Arrow schema.

        :param schema_str: the JSON Schema text
        :type schema_str: str
        :return: the Arrow schema
        :rtype: pa.Schema
        """
        return self.convert(load_schema_text(schema_str))

    def convert_from_file(self, path: str) -> pa.Schema:
        """Convert a JSON or YAML file to an Arrow schema.

        :param path: path to the JSON Schema file
        :type path: str
        :return: the Arrow schema
        :rtype: pa.Schema
        """
        return self.convert(load_schema_file(path))
