"""JSON Schema input guards that raise the caller's direction-specific error."""

from typing import Any

from cdm_data_loaders.converters.core.errors import ConversionError


def require_schema_keyword(
    schema: dict[str, Any], error: type[ConversionError], *, converter_name: str = "Converter"
) -> None:
    """Reject a schema with no top-level '$schema' keyword.

    :param schema: the input JSON Schema document
    :type schema: dict[str, Any]
    :param error: the direction-specific error class to raise
    :type error: type[ConversionError]
    :param converter_name: converter name used in the diagnostic
    :type converter_name: str
    :raises ConversionError: if the schema has no '$schema' keyword
    """
    if not schema.get("$schema"):
        err_msg = (
            f"Input JSON Schema is missing a '$schema' keyword. {converter_name} requires schemas "
            "to explicitly declare their dialect via '$schema'; it will not assume a default."
        )
        raise error(err_msg)


def require_object_root(schema: dict[str, Any], error: type[ConversionError], target: str) -> None:
    """Reject a schema whose declared root type isn't literally 'object' or absent.

    This is intentionally stricter than `readers.json_schema.JsonSchemaReader`,
    which also accepts a list-form root type (e.g. `["object", "null"]`) when
    used directly. Facade converters call this guard first, so a nullable
    root type is only reachable through direct reader use, not through
    `convert()`.

    :param schema: the input JSON Schema document
    :type schema: dict[str, Any]
    :param error: the direction-specific error class to raise
    :type error: type[ConversionError]
    :param target: what the root maps to in the error message (e.g. 'a dlt table')
    :type target: str
    :raises ConversionError: if the root 'type' is declared and is not 'object'
    """
    if schema.get("type") not in (None, "object"):
        err_msg = f"Root schema must be of type 'object' to map to {target}, got: {schema.get('type')!r}"
        raise error(err_msg)


def reject_unresolved_references(
    schema: dict[str, Any],
    error: type[ConversionError],
    *,
    converter_name: str = "Converter",
    dereference_function: str = "dereference_schema",
) -> None:
    """Guard against un-dereferenced schemas reaching the type-dispatch logic.

    $ref (including references into external JSON Schema documents) and allOf
    resolution are intentionally not performed here; run
    `dereference_schema()` first.

    :param schema: the schema fragment to check for unresolved '$ref'/'allOf'
    :type schema: dict[str, Any]
    :param error: the direction-specific error class to raise
    :type error: type[ConversionError]
    :param converter_name: converter name used in the diagnostic
    :type converter_name: str
    :param dereference_function: qualified dereferencing function name
    :type dereference_function: str
    :raises ConversionError: if the schema still contains '$ref' or 'allOf'
    """
    if "$ref" in schema:
        err_msg = (
            f"Encountered an unresolved $ref {schema['$ref']!r}. {converter_name} requires a fully "
            "dereferenced schema -- this includes references to external JSON Schema documents. Use "
            f"`{dereference_function}()` to resolve all $refs before "
            "calling `convert()`."
        )
        raise error(err_msg)
    if "allOf" in schema:
        err_msg = (
            f"Encountered an unmerged 'allOf'. {converter_name} requires a fully dereferenced schema "
            f"with 'allOf' already merged. Use `{dereference_function}()` "
            "before calling `convert()`."
        )
        raise error(err_msg)
