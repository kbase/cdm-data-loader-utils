"""JSON Schema registry for individual assembly reports in the reference JSONL."""

import json
from importlib.resources import files

_response_schema = json.loads(
    files("cdm_data_loaders.parsers.ncbi.datasets_api")
    .joinpath("DatasetReportResponse.schema.json")
    .read_text(encoding="utf-8")
)

ENTITY_SCHEMAS = {
    "dataset": {**_response_schema, "$ref": "#/$defs/v2reportsAssemblyDataReport"},
}
