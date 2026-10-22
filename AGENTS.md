# CLAUDE.md

Guidance for working in this repository. Rules with stricter, active requirements live in `.claude/rules/` (coding-style, testing, heresy) and apply to all changes.

## Project

`cdm-data-loaders`: data loaders and wranglers for the KBase CDM. Python >= 3.13, managed and executed with `uv`. Package source is under `src/cdm_data_loaders/` (~115 files, 8 top-level packages). Tests mirror the source layout under `tests/cdm_data_loaders/` (including pipeline end-to-end suites under `tests/cdm_data_loaders/pipelines/`), with multi-phase NCBI FTP integration tests in `tests/integration/ncbi_ftp/` and vcrpy cassettes in `tests/cassettes/`.

Autouse fixtures in `tests/conftest.py` (`_isolated_cli_and_env` and `_isolated_dlt_config` from `tests/dlt_config_isolation.py`) reset `sys.argv = ["pytest"]`, strip `CDL_*` environment variables, and snapshot/restore `dlt.config` between tests. Pipeline end-to-end tests use `PipelineTestRunner` (`tests/cdm_data_loaders/pipelines/framework.py`) to run pipelines across `ExecutionMode.CORE`, `ExecutionMode.CLI`, and `ExecutionMode.DOCKER` (Docker tests are marked `requires_docker` and require `PIPELINE_DOCKER_IMAGE` to be set).

Useful invocations:

- Standard local test run (skips Spark, Ceph, and unset Docker image tests): `uv run pytest -m "not requires_spark and not requires_ceph"`
- Full pipelines suite, excluding the bulk downloader: `uv run pytest tests/cdm_data_loaders/pipelines --ignore=tests/cdm_data_loaders/pipelines/all_the_bacteria`
- Unit vs. container test script: `./scripts/run_tests.sh [unit|container|all]`
- Lint: `uv run ruff check` (line length 120, target py313, full rule selection in `pyproject.toml`).

## Architecture

Domain logic lives in the top-level packages and is orchestration-free; orchestration lives in `pipelines/`. Data flows source -> reader -> intermediate representation (IR) -> emitter -> target:

- Converters: reader (source schema -> IR) + emitter (IR -> target format) + thin facade composing them. See `docs/converters_ir.md` for the full IR contract and package layout.
- IR is a frozen, slotted node tree (`converters/core/ir.py`) that preserves source facts without choosing target policy; emitters make all target-specific decisions.
- `converters/core/` imports no dlt, PyIceberg, or PySpark modules. The only cross-package import is the `x-xsv-config` contract reusing `X_XSV_CONFIG_SCHEMA` from `cdm_data_loaders.readers.jsonschema_xsv.xsv_validator.custom_metaschema`.

Settings hierarchy: pipeline settings in `pipelines/` extend the hierarchy in `core/settings.py`: `CdmDataLoadersBase` -> `LoggerSettings` -> `LocalInputOutputSettings` -> `CtsSettings` -> `FileFormatIngestSettings`. Settings use `env_prefix="CDL_"`, `cli_kebab_case=True`, and `default_settings_with_shortcuts(...)`. `CdmDataLoadersBase.__pydantic_init_subclass__` calls `check_aliases()` at class creation time to verify that CLI argument names, aliases, and `cli_shortcuts` do not collide or claim reserved flags (`h`, `help`). CLI parsing is invoked explicitly via `pydantic_settings.CliApp.run` in `pipelines/core.py` (`PipelineCli` / `run_cli`), not implicitly on model instantiation. Every pipeline registers a `cli` entry point in `pyproject.toml` under `[project.scripts]` (and in `scripts/entrypoint.sh` for Docker execution): `all_the_bacteria`, `jsonl_extract`, `jsonl_jsonschema`, `jsonl_pydantic`, `ncbi_ftp_promote`, `ncbi_ftp_sync`, `ncbi_rest_api`, `pdb_manifest`, `uniprot`, `uniref`, `xmltodict`, `xml2db`, `xsv`.

## Package Map

All paths below are relative to `src/cdm_data_loaders/`.

### `core/`

Framework-level foundations:

- `fields.py`: common annotated field types (`BufferSize`, `DatasetName`, `FileGlob`, `LoaderFileFormat`, `TableName`, `WriteDisposition`, etc.), `LoaderFileFormatEnum` (`parquet`, `jsonl`), and `DEFAULTS`.
- `path_fields.py`: `LocalInputDir` and `LocalOutputDir` path types (`file://` or plain local directory paths) and URL/path validation helpers (`normalize_local_path`, `validate_directory`, `validate_url_protocol`).
- `settings.py`: `pydantic-settings` base class hierarchy (`CdmDataLoadersBase` -> `LoggerSettings` -> `LocalInputOutputSettings` -> `CtsSettings` -> `FileFormatIngestSettings`), `CLI_SHORTCUTS`, `default_settings_with_shortcuts()`, and CLI flag/shortcut collision checking (`check_aliases()`).

### `converters/`

Schema-format translation. See `docs/converters_ir.md` for the full contract.

- `core/`: shared, target-independent foundation
  - `ir.py` (NodeType, Field, SchemaDocument, Reader/Emitter protocols), `ir_values.py` (`freeze_value` and friends), `errors.py` (ConversionError), `extensions.py` (ExtensionRegistry, DEFAULT_EXTENSIONS), `guards.py` (input guards), `inference.py` (type inference, decimal scale), `metadata.py` (draft-specific JSON Schema keyword selection), `io.py` (schema loading), `paths.py` (`set_nested`).
- `readers/`: source -> IR
  - `json_schema.py` (JsonSchemaReader), `dlt.py` (DltReader over dlt's normalized table forest), `iceberg.py` (IcebergReader).
- `emitters/`: IR -> target
  - `json_schema.py`, `dlt.py`, `pyspark.py`, `pyarrow.py`, `linkml.py`.
- Facades: `jsonschema_to_dlt.py`, `jsonschema_to_pyspark.py`, `jsonschema_to_pyarrow.py`, `dlt_to_jsonschema.py`, `pyiceberg_to_jsonschema.py`, `pyiceberg_to_linkml.py`, plus CLI dumps `pyiceberg_catalog_to_jsonschema.py` and `pyiceberg_catalog_to_linkml.py`.

### `readers/`

File-format readers feeding pipelines:

- `jsonlines.py`: JSONL line parsing (`parse_jsonl_line`), paged streaming (`stream_jsonl_lines`), and table routing (`route_jsonl_lines`).
- `xml.py`: shared XML streaming and `xmltodict` helpers (`parse_head_matter`, `stream_xml_file`, `process_xml_file_to_dict`, `process_xml_file`, `process_xml_file_batches`, `process_xml_file_items`, and `build_xml_file_resource` over dlt's `filesystem` source).
- `xml2db_doc.py`: `xml2db` model building (`build_xml2db_model`), bounded-memory chunk fragment streaming (`iter_xml2db_chunk_fragments`), primary/foreign key rewriting (content-addressed keys for reused tables, file/chunk-namespaced keys otherwise), document flattening (`flatten_xml2db_document`), and `process_xml_file_with_xml2db`.
- `xsd.py`: XSD schema loading (`load_schema`) and parent/child element cardinality classification (`find_list_and_single_child_paths`).
- `jsonschema_xsv/`: `xsv_reader.py` (casts validated CSV/TSV rows to JSON Schema property types and resolves `x-xsv-config` parsing parameters) plus `xsv_validator/` (`qsv.py` wraps the external `qsv` binary for two-pass validation, normalization, and null replacement; `custom_metaschema.py` defines the `https://json-schema.kbase.us/draft/2020-12-xsv/schema` metaschema and `X_XSV_CONFIG_SCHEMA`, also imported by `converters`; `schema_utils.py` and `helpers.py`).

### `parsers/`

Parsers for third-party outputs and formats:

- `uniprot/`: `uniprot_kb.py` (UniProt XML), `uniref.py` (UniRef XML), `metalink.py`, `relnotes.py`, `idmapping.schema.json`; shared XML helpers in `utils/xml_utils.py`.
- `ncbi/datasets_api/openapi_pydantic.py`: pydantic models for NCBI Datasets API responses (bundled `AnnotationReportResponse.schema.json` and `DatasetReportResponse.schema.json` alongside).
- `ncbi/ftp_mod_times/`: `ftp_mod_times.schema.json` for NCBI FTP modification-time files.
- `bbmap_stats.py` (bbmap `stats.sh` output), `checkm2.py` (quality reports), `pdb_holdings_file.py`.
- `all_the_bacteria/schema/` and `gtdb/schema/`: JSON Schemas for ATB and GTDB tabular files (no Python modules).

### `ncbi_ftp/`

NCBI FTP assembly domain logic, phased: `manifest.py` (phase 1: assembly summary diffing and manifest generation), download is the phase 2 pipeline, `promote.py` (phase 3: promote staged files to final Lakehouse paths in S3), plus `assembly.py` (assembly-specific logic), `metadata.py` (frictionless data package descriptors), `constants.py`.

### `pdb/`

PDB pipeline constants, enums, and globals (`constants.py`); the pipeline itself is `pipelines/pdb_manifest.py`.

### `pipelines/`

dlt-based orchestrators; each exposes a `cli = PipelineCli(...)` registered in `pyproject.toml`:

- `core.py`: `PipelineCli`, `run_cli`, `run_pipeline`, `setup_dlt_Pipeline`, `execute_dlt_pipeline`, and `configure_pipeline_logging`.
- `all_the_bacteria.py`: All The Bacteria bulk data downloader (excluded from fast test runs).
- `jsonlines/`: JSONL ingest pipelines (`settings.py`, shared machinery in `base.py`) - `extract_pipeline.py` (`jsonl_extract`, plain load into `<table_name>` / `<table_name>_invalid`), `extract_jsonschema_validate_pipeline.py` (`jsonl_jsonschema`, validates against `ENTITY_SCHEMAS`), `extract_pydantic_validate_pipeline.py` (`jsonl_pydantic`, validates against `ENTITY_MODELS`). See `docs/jsonschema_ingest.md`.
- `ncbi_ftp_download.py` / `ncbi_ftp_promote.py`: phases 2 (`ncbi_ftp_sync`) and 3 (`ncbi_ftp_promote`) of the NCBI FTP flow.
- `ncbi_rest_api.py`: NCBI REST API import.
- `pdb_manifest.py`: PDB holdings transfer manifest.
- `uniprot_kb.py` / `uniref.py`: UniProt (`uniprot`) and UniRef (`uniref`) dlt pipelines.
- `xml2db/` / `xmltodict/` / `xsv/`: XML-to-relational, XML-to-dictionary, and XSV (CSV/TSV) ingestion, each with `pipeline.py` and `settings.py`; `xml2db/compaction.py` runs DuckDB-based post-load deduplication of "reused" tables on filesystem destinations.

### `utils/`

Cross-cutting helpers: `batcher.py` (`NumericFileSequenceBatcher`, `get_file_batches`), `buffer.py` (per-table row buffers for dlt: `DictBuffer`, `ListBuffer`, `ArrowDictBuffer`, `ArrowListBuffer`), `cdm_logger.py` (structured JSON logging), `checksums.py`, `gz.py` (`open_maybe_gzip`), `helpers.py`, `progress.py` (tqdm helpers), `spark.py`, `xml_utils.py`, `jsonschema/dereferencer.py` (`$ref` resolution), `dlt/dlt_normalization.py` (lossless parent-linked normalization of dlt table definitions, used by the converters' `DltReader`), and `file_transfer/` - `core.py` (tenacity-retried downloads), `sync_client.py` / `async_client.py` / `ftp_client.py`, and `s3/` (`client.py`, `object_utils.py`, HTTP-to-S3 streaming in `streamer.py`, `transfer_config.py`, `cli.py`).

### Other Top-Level Directories

- `scripts/`: `entrypoint.sh` (container entrypoint dispatching to registered pipeline CLIs), `run_tests.sh` (`unit`, `container`, or `all` test suites), `filter_ncbi_updates.py` (CLI filtering NCBI assembly IDs by modification date), `s3_local.py` (thin S3 CLI for local Ceph testing), `dtd2xsd.pl` (DTD-to-XSD converter), `run_tools.sh`, `parse_index.py`.
- `docs/`: architecture and pipeline guides (`converters_ir.md`, `jsonschema_ingest.md`, `ncbi_ftp_e2e_walkthrough.md`, `pdb_e2e_walkthrough.md`, `rclone_instructions.md`, `uniprot.md`).
- `cdl_skill/`: skill definition (`SKILL.md`) and generator (`scripts/skill_gen.sh`) for packaging pipeline documentation and schemas.

## Conventions

- Frozen, slotted dataclasses for models; constructors validate and recursively copy mappings.
- CLI argument parsing by pydantic-settings
- Fully type-annotated code (ruff `ANN`); datetime values carry an explicit timezone (`DTZ`); use `Decimal` for exact arithmetic, never constructed from floats.
- External HTTP calls in tests are recorded and replayed with vcrpy; AWS is mocked with moto; Spark and Ceph tests run in containers via `docker-compose.yml` and are skipped by default (`-m "not requires_spark and not requires_ceph"`).
- Tests are module-based `def test_*() -> None:` functions with pytest fixtures in `conftest.py`, named `test_<function>_<pass|fail>_<conditions>`.

## Gotchas

- dlt 1.30 DuckDB destination fails with a Binder Error "Ambiguous reference to catalog or schema" when the db filename stem equals the dataset name; keep the db filename distinct from both pipeline and dataset names.
- dlt 1.30: `AnyDestination` is a quoted TypeAlias (a str at runtime), so `AnyDestination | None` in an evaluated annotation raises TypeError on import. Use `TDestinationReferenceArg` from `dlt.common.destination.reference`.
- pyiceberg SqlCatalog engines leak sqlite connections; use the context-manager protocol (`with load_catalog(...) as catalog`) so `close()` disposes the engine.
- pipelines are run on the KBase CTS (Task Service), which only saves the contents of the output directory (or dlt destination) if the container exits with exit code 0. For any other error code, the container STDERR and STDOUT are saved, but nothing else. As a result, even if pipelines experience failure (e.g. `load_info.has_failed_packages` is true), the container must exit 0 to ensure that logs, parsed data, and other artefacts are not deleted. Errors that occur before a pipeline is run (e.g. argument parsing errors) are not caught, so they should be picked up by the container STDERR/STDOUT.
