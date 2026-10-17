"""JSONL validate-and-load pipeline for the KBase CTS.

Reads JSONL files. Validates each record against a Pydantic model.
Writes valid and invalid records to separate tables.

The entity models come from an external module set by
`entity_models_module`. That module must define:

    ENTITY_MODELS: dict[str, type[BaseModel]]

Each key is a table name. Each value is the Pydantic model for that
table. Example:

    from datetime import date
    from pydantic import BaseModel

    class Sample(BaseModel):
        id: str
        collection_date: date

    ENTITY_MODELS: dict[str, type[BaseModel]] = {"sample": Sample}

`input_dir` must contain one subdirectory per table name in
`ENTITY_MODELS`. Each subdirectory holds that table's JSONL files.
Files may be plain text or gzip-compressed (`.jsonl` or `.jsonl.gz`).

Valid records go to a table named after the entity. Records that fail
JSON parsing or model validation go to a `<table>_rejected` table. Each
rejected record includes source file, line number, and error detail.
Write disposition is `append` or `replace` only. No merge is performed.

Each model also sets its table's schema: it is passed to dlt as the table's
column hints, so destination types come from the model rather than inference.
"""

import json
from importlib import import_module
from typing import Any

from dlt.common.pipeline import LoadInfo
from pydantic import BaseModel, ValidationError

from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.pipelines.jsonlines.base import (
    ParseErrorResolver,
    run_registry_ingest,
)
from cdm_data_loaders.pipelines.jsonlines.base import (
    build_entity_resource as build_shared_entity_resource,
)
from cdm_data_loaders.pipelines.jsonlines.settings import (
    PYDANTIC_PIPELINE_NAME,
    JsonlPydanticIngestSettings,
)


def load_entity_models(module_path: str) -> dict[str, type[BaseModel]]:
    """Import module_path and return its ENTITY_MODELS dict.

    :param module_path: dotted import path of the module defining ENTITY_MODELS
    :type module_path: str
    :return: mapping of table name to Pydantic model
    :rtype: dict[str, type[BaseModel]]
    :raises ValueError: if ENTITY_MODELS is missing, empty, or not a dict
    """
    module = import_module(module_path)
    registry = getattr(module, "ENTITY_MODELS", None)
    if not isinstance(registry, dict) or not registry:
        err_msg = f"{module_path} does not define a non-empty ENTITY_MODELS: dict[str, type[BaseModel]]"
        raise ValueError(err_msg)
    return registry


def _make_resolve_error(model: type[BaseModel]) -> ParseErrorResolver:
    """Build the per-record Pydantic validation closure.

    A record that fails `model.model_validate` yields its error list as JSON.
    """

    def _validate(record: dict[str, Any]) -> str | None:
        try:
            model.model_validate(record)
        except ValidationError as e:
            return json.dumps(e.errors(), default=str)
        return None

    return _validate


def build_entity_resource(table_name: str, model: type[BaseModel], settings: JsonlPydanticIngestSettings) -> Any:  # noqa: ANN401
    """Build the resource that reads, validates, and routes one entity's records.

    Reads files matching `settings.file_glob` from
    `settings.input_dir/table_name`. If that directory does not exist,
    uses an empty file list instead. A missing entity then produces zero
    rows, not an error. The model sets the entity table's dlt columns.

    :param table_name: entity name; also the input_dir subdirectory name
    :type table_name: str
    :param model: Pydantic model used to validate this entity's records
    :type model: type[BaseModel]
    :param settings: pipeline settings
    :type settings: JsonlPydanticIngestSettings
    :return: resource yielding validated and rejected records for this entity
    :rtype: Any
    """

    def _dump_validated(record: dict[str, Any]) -> dict[str, Any]:
        return model.model_validate(record).model_dump(mode="python")

    return build_shared_entity_resource(
        table_name,
        settings,
        resolve_error=_make_resolve_error(model),
        include_record=True,
        transform_valid=_dump_validated,
        table_columns=model,
    )


def run_jsonlines_with_validation_pipeline(settings: JsonlPydanticIngestSettings) -> LoadInfo | None:
    """Run the JSONL validate-and-load pipeline for every requested table.

    :param settings: pipeline configuration
    :type settings: JsonlPydanticIngestSettings
    :raises ValueError: if settings.table_names includes a name not in ENTITY_MODELS
    :return: load information for the pipeline
    :rtype: LoadInfo | None
    """
    entity_models = load_entity_models(settings.entity_models_module)
    return run_registry_ingest(
        settings,
        entity_models,
        build_resource=build_entity_resource,
        pipeline_name=PYDANTIC_PIPELINE_NAME,
        pipeline_run_kwargs={"loader_file_format": str(settings.loader_file_format)},
    )


cli = PipelineCli(JsonlPydanticIngestSettings, run_jsonlines_with_validation_pipeline)


if __name__ == "__main__":
    cli()
