"""Fixed legacy outputs for the four converter routes before facade rewiring."""

import json
from pathlib import Path
from typing import Any

import pytest
from pyiceberg.schema import Schema

from cdm_data_loaders.converters.dlt_to_jsonschema import DltToJSONSchema
from cdm_data_loaders.converters.jsonschema_to_dlt import JSONSchemaToDlt
from cdm_data_loaders.converters.jsonschema_to_pyspark import JSONSchemaToPySpark
from cdm_data_loaders.converters.pyiceberg_to_jsonschema import table_to_json_schema
from cdm_data_loaders.converters.readers.dlt import DltReader
from cdm_data_loaders.converters.readers.iceberg import IcebergReader
from cdm_data_loaders.converters.readers.json_schema import JsonSchemaReader
from tests.cdm_data_loaders.converters.conftest import make_table

GOLDEN_PATH = Path(__file__).parents[3] / "data" / "converters" / "ir" / "legacy_outputs.json"
CASES: list[dict[str, Any]] = json.loads(GOLDEN_PATH.read_text(encoding="utf-8"))


def legacy_result(case: dict[str, Any]) -> dict[str, Any]:
    """Run a legacy converter over one fixed corpus input."""
    source, options = case["input"], case["options"]
    match case["route"]:
        case "json-dlt":
            return JSONSchemaToDlt(**options).convert(source)
        case "json-spark":
            return JSONSchemaToPySpark(**options).convert(source).jsonValue()
        case "dlt-json":
            return DltToJSONSchema(**options).convert(source)
        case "iceberg-json":
            identifier = tuple(source["identifier"])
            table = make_table(identifier, Schema.model_validate(source["schema"]), properties=source["properties"])
            return table_to_json_schema(table, identifier)
        case _:
            message = f"Unknown golden route: {case['route']}"
            raise ValueError(message)


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_legacy_result_pass_fixed_output(case: dict[str, Any]) -> None:
    """Existing converters retain exact structured output on the captured corpus."""
    assert legacy_result(case) == case["expected"]


@pytest.mark.parametrize("case", CASES, ids=[case["id"] for case in CASES])
def test_read_pass_golden_inputs(case: dict[str, Any]) -> None:
    """Each golden source is readable as a structural IR tree."""
    source = case["input"]
    if case["route"].startswith("json-"):
        document = JsonSchemaReader().read(source)
        assert tuple(prop.name for prop in document.root.properties) == tuple(source["properties"])
    elif case["route"] == "dlt-json":
        options = {key: value for key, value in case["options"].items() if key != "preserve_unknown_hints"}
        documents = DltReader(**options).read(source)
        assert tuple(documents) == tuple(case["expected"])
    else:
        schema = Schema.model_validate(source["schema"])
        document = IcebergReader().read(schema)
        assert tuple(prop.name for prop in document.root.properties) == tuple(field.name for field in schema.fields)
