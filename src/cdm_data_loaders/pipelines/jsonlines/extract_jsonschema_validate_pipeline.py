"""Load JSONL records validated against an external JSON Schema registry.

The module selected by ``schema_files_module`` must define ``ENTITY_SCHEMAS``:
a mapping of table names to JSON Schema dictionaries, each declaring ``$schema``.
Files are read from one input subdirectory per table. Valid records are loaded
unchanged; parsing and validation failures go to ``<table>_rejected`` with
source_file, line_no, raw_record, and error_detail. Format checks are enabled.

Each JSON Schema also sets its table's schema: it is dereferenced, converted to
dlt columns through the JSON Schema to PyArrow converter, and applied as a dlt
table variant, so destination types come from the schema rather than inference.
"""

import json
from importlib import import_module
from logging import Logger, getLogger
from typing import Any

from dlt.common.libs.pyarrow import py_arrow_to_table_schema_columns
from dlt.common.pipeline import LoadInfo
from dlt.common.schema.typing import TTableSchemaColumns
from dlt.extract.resource import DltResource
from jsonschema import FormatChecker
from jsonschema.exceptions import SchemaError
from jsonschema.validators import validator_for

from cdm_data_loaders.converters.jsonschema_to_pyarrow import JSONSchemaToPyArrow
from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.pipelines.jsonlines.base import (
    ParseErrorResolver,
    run_registry_ingest,
)
from cdm_data_loaders.pipelines.jsonlines.base import (
    build_entity_resource as build_shared_entity_resource,
)
from cdm_data_loaders.pipelines.jsonlines.settings import JSONSCHEMA_PIPELINE_NAME, JsonlJsonschemaIngestSettings
from cdm_data_loaders.utils.jsonschema.dereferencer import dereference_schema

logger: Logger = getLogger(__name__)


def load_entity_schemas(module_path: str) -> dict[str, dict[str, Any]]:
    """Import and validate a non-empty ENTITY_SCHEMAS registry.

    Raises ValueError for malformed registries, missing draft declarations, or
    unsupported drafts, TypeError for invalid schema types, and RuntimeError
    listing invalid schemas by table.

    :param module_path: dotted import path of the module defining ENTITY_SCHEMAS
    :type module_path: str
    :return: mapping of table name to JSONSchema model
    :rtype: dict[str, dict[str, Any]]
    :raises ValueError: if ENTITY_SCHEMAS is missing, empty, or not a dict
    """
    module = import_module(module_path)
    registry = getattr(module, "ENTITY_SCHEMAS", None)
    if not isinstance(registry, dict) or not registry:
        err_msg = f"{module_path} does not define a non-empty ENTITY_SCHEMAS"
        raise ValueError(err_msg)
    errors: dict[str, str] = {}
    for table_name, schema in registry.items():
        if not isinstance(table_name, str) or not table_name.strip():
            err_msg = "ENTITY_SCHEMAS keys must be non-empty table names"
            raise ValueError(err_msg)
        if not isinstance(schema, dict):
            err_msg = f"JSON Schema for {table_name} must be a dict"
            raise TypeError(err_msg)
        if "$schema" not in schema:
            err_msg = f"JSON Schema for {table_name} is missing the $schema keyword"
            raise ValueError(err_msg)
        if not isinstance(schema["$schema"], str):
            err_msg = f"JSON Schema for {table_name} must declare a $schema string"
            raise TypeError(err_msg)
        validator = validator_for(schema)
        declared = schema["$schema"]
        metaschema = getattr(validator, "META_SCHEMA", None)
        declared_supported = metaschema is None or "$schema" not in metaschema or metaschema["$schema"] == declared
        if not declared_supported:
            err_msg = f"Unsupported JSON Schema draft for {table_name}: {declared}"
            raise ValueError(err_msg)
        try:
            validator.check_schema(schema)
        except SchemaError as error:
            logger.exception("Error validating JSON Schema")
            errors[table_name] = str(error)

    if errors:
        details = "".join(f"{table_name}: {error}\n" for table_name, error in errors.items())
        err_msg = f"The following errors were found when validating the schemas supplied in {module_path}:\n{details}"
        raise RuntimeError(err_msg)
    return registry


def _make_resolve_error(schema: dict[str, Any]) -> ParseErrorResolver:
    """Build the per-record JSON Schema validation closure.

    Schema failures report every broken rule as a JSON list of messages.
    """
    validator = validator_for(schema)(schema, format_checker=FormatChecker())

    def _validate(record: dict[str, Any]) -> str | None:
        errors = sorted(validator.iter_errors(record), key=str)
        if errors:
            return json.dumps([error.message for error in errors])
        return None

    return _validate


def table_columns_from_schema(schema: dict[str, Any]) -> TTableSchemaColumns:
    """Convert a JSON Schema to dlt column schemas for table hints.

    The schema is dereferenced first, then converted through the JSON Schema to
    PyArrow converter and translated back to dlt columns. Nested-type hints are
    dropped: these pipelines yield Python objects, and dlt's parquet writer
    serializes nested values to strings, so nested types must stay json.
    """
    arrow_schema = JSONSchemaToPyArrow().convert(dereference_schema(schema))
    columns = py_arrow_to_table_schema_columns(arrow_schema)
    return {
        name: {key: value for key, value in column.items() if key != "x-nested-type"}
        for name, column in columns.items()
    }


def build_entity_resource(
    table_name: str, schema: dict[str, Any], settings: JsonlJsonschemaIngestSettings
) -> DltResource:
    """Read plain or gzip JSONL files for one entity and route validated pages.

    A missing entity directory produces no rows. Nested values are preserved
    rather than normalized into child tables. The schema sets the entity
    table's dlt columns after dereferencing.

    :param table_name: table name
    :type  table_name: str
    :param schema: JSONSchema for the table
    :type  schema: dict[str, Any]
    :param settings: pipeline configuration
    :type  settings: JsonlJsonschemaIngestSettings
    :return: resource yielding validated and rejected records for this entity
    :rtype: DltResource
    """
    return build_shared_entity_resource(
        table_name,
        settings,
        resolve_error=_make_resolve_error(schema),
        table_columns=table_columns_from_schema(schema),
    )


def run_jsonlines_pipeline(settings: JsonlJsonschemaIngestSettings) -> LoadInfo | None:
    """Validate the registry and load requested tables using the shared runner.

    :param settings: pipeline configuration
    :type settings:  JsonlJsonschemaIngestSettings
    :raises ValueError: if settings.table_names includes a name not in ENTITY_SCHEMAS
    :return: load information for the pipeline
    :rtype:  LoadInfo | None
    """
    table_to_schema = load_entity_schemas(settings.schema_files_module)
    return run_registry_ingest(
        settings,
        table_to_schema,
        build_resource=build_entity_resource,
        pipeline_name=JSONSCHEMA_PIPELINE_NAME,
        pipeline_run_kwargs={"loader_file_format": str(settings.loader_file_format)},
    )


cli = PipelineCli(JsonlJsonschemaIngestSettings, run_jsonlines_pipeline)


if __name__ == "__main__":
    cli()
