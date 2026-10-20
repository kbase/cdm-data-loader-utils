# JSONL Ingestion Pipelines

Three pipelines under `pipelines/jsonlines/` load JSONL files through the
shared CTS runner. They share settings, file handling, and output options, and
differ only in how (or whether) records are validated before loading:

| Pipeline                                                                                                                       | CLI name           | Validation                 |
| ------------------------------------------------------------------------------------------------------------------------------ | ------------------ | -------------------------- |
| [extract_pipeline.py](../src/cdm_data_loaders/pipelines/jsonlines/extract_pipeline.py)                                         | `jsonl_extract`    | None beyond valid JSON     |
| [extract_jsonschema_validate_pipeline.py](../src/cdm_data_loaders/pipelines/jsonlines/extract_jsonschema_validate_pipeline.py) | `jsonl_jsonschema` | JSON Schema, per table     |
| [extract_pydantic_validate_pipeline.py](../src/cdm_data_loaders/pipelines/jsonlines/extract_pydantic_validate_pipeline.py)     | `jsonl_pydantic`   | Pydantic models, per table |

All three support plain and gzip-compressed inputs and Parquet or JSONL output.

## Extract-Only Pipeline

`jsonl_extract` loads every record from `input_dir` into a single destination
table named by `--table-name`, with no validation beyond parsing as JSON.
Files may live anywhere under `input_dir`; there is no per-table subdirectory
layout. Records that fail JSON parsing go to `<table_name>_invalid`, with
source file, line number, raw text, and parse error detail.

```sh
DESTINATION__LOCAL_FS__DESTINATION_TYPE=filesystem \
python -m cdm_data_loaders.pipelines.jsonlines.extract_pipeline \
  --input-dir ./input \
  --output-dir ./output \
  --dataset-name raw_widgets \
  --table-name widget \
  --use-destination local_fs \
  --use-output-dir-for-pipeline-metadata true \
  --loader-file-format parquet
```

## The Two Validated Pipelines

`jsonl_jsonschema` and `jsonl_pydantic` share the same input layout and
rejection behavior; only the validation mechanism and registry format differ.

### JSON Schema Registry

Provide an importable Python module, such as `my_schemas`, defining
`ENTITY_SCHEMAS`. Each key names a table and its input subdirectory. Each value
is a JSON Schema dictionary, not a path to a schema file:

```python
ENTITY_SCHEMAS = {
    "widget": {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "type": "object",
        "required": ["widget_id", "count"],
        "properties": {
            "widget_id": {"type": "string", "minLength": 1},
            "count": {"type": "integer", "minimum": 0},
            "collected": {"type": "string", "format": "date"},
        },
    },
}
```

Every schema must declare a supported `$schema` draft. The entire registry is
checked before extraction, including schemas for unselected tables. Missing or
unknown drafts fail explicitly; invalid schemas are reported together by table.
Validation uses each schema's declared draft and enables `jsonschema` format
checks. Unknown format names are ignored by that library.

```sh
DESTINATION__LOCAL_FS__DESTINATION_TYPE=filesystem \
python -m cdm_data_loaders.pipelines.jsonlines.extract_jsonschema_validate_pipeline \
  --schema-files-module my_schemas \
  --input-dir ./input \
  --output-dir ./output \
  --dataset-name validated_widgets \
  --use-destination local_fs \
  --use-output-dir-for-pipeline-metadata true \
  --log-config-file logger_config.json \
  --loader-file-format parquet
```

The `schema_files_module` option name is retained for compatibility, but it
references `ENTITY_SCHEMAS` dictionaries. The module must be on Python's import
path. The dataset name is required and must be nonempty.

### Pydantic Registry

Provide an importable Python module, such as `my_models`, defining
`ENTITY_MODELS`. Each key names a table and its input subdirectory. Each value
is the Pydantic model used to validate that table's records:

```python
from datetime import date
from pydantic import BaseModel


class Widget(BaseModel):
    widget_id: str
    count: int
    collected: date | None = None


ENTITY_MODELS = {"widget": Widget}
```

Valid records are loaded via `model.model_dump(mode="python")`, so type
coercion and field defaults from the model apply, unlike the JSON Schema
pipeline. Write disposition is `append` or `replace` only; no merge is
performed.

```sh
DESTINATION__LOCAL_FS__DESTINATION_TYPE=filesystem \
python -m cdm_data_loaders.pipelines.jsonlines.extract_pydantic_validate_pipeline \
  --entity-models-module my_models \
  --input-dir ./input \
  --output-dir ./output \
  --dataset-name validated_widgets \
  --use-destination local_fs \
  --use-output-dir-for-pipeline-metadata true \
  --log-config-file logger_config.json \
  --loader-file-format parquet
```

### Inputs And Outputs

```text
input/
  widget/
    first.jsonl
    second.jsonl.gz
```

- Each table has its own input subdirectory, matching `--file-glob`
  (default `*.jsonl*`).
- Valid records go to a table named after the entity. The JSON Schema pipeline
  does not coerce types, add defaults, or remove extra properties; the
  Pydantic pipeline applies its model's coercion and defaults.
- Parsing and validation failures go to `<table>_rejected`.
- Rejected rows contain `source_file`, `line_no`, `raw_record`, and `error_detail`.
- `source_file` is relative to the entity directory. `line_no` counts physical
  lines, including blank lines. `raw_record` preserves the reader's stripped
  source text, without reserializing parsed JSON.
- Parse errors use a plain error string. Validation errors use a JSON-encoded
  list containing all validation error messages (JSON Schema messages, or
  Pydantic's `error.errors()` dicts).
- Nested values stay in the parent table. Parquet represents them as JSON
  strings; JSONL retains the nested structures. dlt still applies its normal
  column naming, inference, and metadata rules.
- Blank lines are skipped. Missing entity directories and empty files produce
  no rows. Duplicate records are retained.
- `--table-names` selects a subset of tables; it defaults to every table in
  the registry. Unknown names raise `ValueError`.

## Shared Options

Standard CTS settings and `CDL_` environment variables apply to all three
pipelines.

| Option                         | Applies to           | Behavior                                                         |
| ------------------------------ | -------------------- | ---------------------------------------------------------------- |
| `--table-name`                 | extract              | Destination table for all records                                |
| `--table-names`, `-t`          | jsonschema, pydantic | Comma-separated selection; defaults to every registered table    |
| `--schema-files-module`, `-m`  | jsonschema           | Importable `ENTITY_SCHEMAS` registry module                      |
| `--entity-models-module`, `-m` | pydantic             | Importable `ENTITY_MODELS` registry module                       |
| `--file-glob`, `-g`            | all                  | File pattern within each input directory; defaults to `*.jsonl*` |
| `--buffer-size`                | all                  | Positive maximum page size                                       |
| `--loader-file-format`         | all                  | `parquet` (default) or `jsonl`                                   |

The registered pipeline names are `jsonlines_ingest` (extract and pydantic) and
`jsonschema_validator` (JSON Schema). Imports from the former benchmark module
are retained as compatibility re-exports; new callers should use the package
module.

## Tests

Run the unit and real filesystem end-to-end tests without external services:

```sh
uv run pytest \
  tests/cdm_data_loaders/pipelines/jsonlines/ \
  tests/integration/pipelines/jsonlines/
```
