# RefSeq Pipeline – Architecture & Developer Guide

The RefSeq Pipeline is a Spark- and Delta Lake–based data ingestion and update system designed to
efficiently track and process genome assembly data from NCBI RefSeq.

Key capabilities include:

- Fetching genome assembly metadata from the NCBI Datasets API
- Tracking content-level changes using hash-based snapshots
- Supporting incremental updates instead of full reprocessing
- Normalizing heterogeneous NCBI responses into a stable Common Data Model (CDM)

This document focuses on the internal architecture, module responsibilities, and execution flow,
rather than end-user CLI usage.


## Design Principles

The RefSeq Pipeline follows a set of explicit design principles to ensure
scalability, reproducibility, and maintainability:

- **Deterministic IDs**  
  CDM identifiers are UUIDv5-based and remain stable across runs given the same input.

- **Incremental by default**  
  Hash-based snapshots are used to detect content-level changes, avoiding unnecessary reprocessing.

- **Pure Spark execution**  
  The core pipeline avoids Pandas and relies exclusively on Spark for scalability.

- **Schema-first design**  
  All outputs strictly conform to the predefined `CDM_SCHEMA`.

- **Separation of concerns**  
  API access, parsing, hashing, and storage are implemented as independent modules.

## Core Modules (Execution Order)
> Note: This order reflects logical dependencies, not standalone execution steps.

### Step 1. config.py (Global Configuration & Schema)

**Responsibility**  
Central configuration and schema definition for the entire pipeline.

**Defines**
- `CDM_NAMESPACE`: UUID namespace for deterministic CDM identifiers
- `NCBI_BASE_V2`: Base URL for the NCBI Datasets API
- `EXPECTED_COLS`: Required output columns for CDM normalization
- `CDM_SCHEMA`: Spark `StructType` defining the CDM data contract

This module is imported by nearly all other core components and establishes the schema and identifier semantics used throughout the pipeline.


### Step 2. datasets_api.py (NCBI Datasets API Client)

**Responsibility**
- Fetch genome assembly dataset reports from the NCBI Datasets API
- Serve as the sole external API ingress for RefSeq metadata
- Handle pagination, retries, and transient API failures
- Stream raw assembly reports for downstream Spark-based processing

**Notes**
The RefSeq pipeline retrieves genome assembly metadata directly from the NCBI Datasets API, which serves as the authoritative and up-to-date source
for RefSeq assembly reports.

This module implements a retry-enabled, session-based API client that iteratively streams assembly reports for a given NCBI Taxonomy ID.
Pagination is handled via API-provided page tokens, and transient failures are mitigated through bounded retries to avoid infinite loops or API abuse.

**API Endpoint**
All requests are made against the NCBI Datasets V2 API:
- Base URL: `https://api.ncbi.nlm.nih.gov/datasets/v2`
- Genome assembly reports:
  - `/genome/taxon/{taxon}/dataset_report`



### Step 3. refseq_io.py (RefSeq FTP & Assembly Index Utilities)

**Responsibility**
- Load and parse the RefSeq assembly index (assembly_summary_refseq.txt)
- Resolve stable mappings:
  - accession → FTP path
  - accession → taxonomic identifiers (taxid, species_taxid)
- Fetch remote content from RefSeq FTP:
  - annotation hash files
  - MD5 checksum files
- Provide normalized, cached access to FTP-based metadata and file content

**Notes**
This module acts as the boundary layer between NCBI’s structured metadata (API responses and assembly summaries) and the RefSeq FTP filesystem.

While upstream modules operate on JSON-based assembly reports, `refseq_io.py` resolves those records into concrete FTP locations and retrieves content used for downstream change detection (hash snapshots).

Network access is centralized and stabilized via shared HTTP sessions, retry logic, and optional caching to minimize redundant downloads.

	
### Step 4. cdm_parse.py (Normalize NCBI Reports into CDM)

**Responsibility**
- Normalize heterogeneous NCBI assembly reports into a stable, schema-aligned CDM representation
- Generate deterministic CDM entity IDs (UUIDv5-based) to ensure cross-run and cross-release stability
- Perform safe and defensive type conversions on numeric and percentage fields
- Bridge raw NCBI JSON structures into Spark-native rows conforming to `CDM_SCHEMA`

**Notes**
NCBI assembly reports contain heterogeneous field naming conventions (e.g. snake_case vs camelCase), optional sections, and loosely typed values.
This module isolates all normalization logic, ensuring that downstream Spark and Delta Lake operations operate on a clean, predictable schema.
By centralizing ID generation, field selection, and type coercion, `cdm_parse.py` guarantees that identical biological entities are consistently mapped to the same CDM identifiers across pipeline runs.


### Step 5. spark_delta.py (Spark & Delta Lake I/O Layer)
**Responsibility**
- Initialize SparkSession with Delta Lake support and metastore integration
- Persist Spark DataFrames into Delta Lake as either:
  - managed tables (metastore-managed)
  - external Delta tables (path-based)
- Enforce schema consistency and controlled schema evolution
- Support append and overwrite semantics with safety checks
- Perform post-write cleanup, deduplication, and optional optimization
- Register external Delta paths into the Spark metastore for SQL access

**Notes**

This module serves as the infrastructure boundary between Spark-based computation and persistent storage.

All Delta Lake–specific behaviors — including schema evolution, overwrite safeguards, deduplication rules, and table lifecycle management — are centralized here to prevent leakage of storage logic into parsing or business logic layers.

By isolating write semantics and cleanup policies, `spark_delta.py` ensures that upstream modules can focus solely on data correctness, while downstream consumers interact with stable, queryable Delta tables.


### Step 6. hashes_snapshot.py (Content Hash Snapshot Generation)

**Responsibility**
- Generate deterministic content fingerprints for RefSeq assemblies
- Fetch remote assembly content from NCBI FTP, including:
  - annotation hash files (preferred)
  - MD5 checksum files (fallback)
- Normalize raw content and compute SHA256 digests
- Materialize hash snapshots as Spark DataFrames suitable for Delta persistence

**Purpose**
Enable content-based change detection at the assembly level.

Rather than relying on timestamps or metadata fields, this module fingerprints the actual biological deliverables (annotation files and checksums) associated with each assembly.
These hash snapshots form the foundation for incremental updates, ensuring that downstream processing is triggered only when the underlying biological content has genuinely changed.



### Step 7. hashes_diff.py (Incremental Change Detection)

**Responsibility**
- Compare two hash snapshots (old vs new) stored in Delta Lake
- Identify assembly-level changes, including:
  - newly introduced assemblies
  - updated assemblies with content changes
  - removed or missing assemblies
- Resolve changed accessions to affected taxonomy IDs

**Purpose**

This module is the decision engine of the incremental pipeline.

By diffing content-based hash snapshots rather than metadata, it determines the minimal set of assemblies and taxa that require reprocessing.
The output of this step directly drives downstream execution, ensuring that only biologically meaningful changes propagate through the system.


### Step 8. snapshot_utils.py (Delta Snapshot Diff Utilities)

**Responsibility**

- Provide lightweight, reusable helpers for comparing two Delta snapshots
- Identify:
  - newly added accessions
  - removed accessions
  - accessions with content changes
- Operate directly on Delta table paths rather than Spark metastore tables

**Purpose**

This module exposes low-level snapshot diff primitives that can be reused by CLI commands, orchestration layers, and ad-hoc workflows.
It deliberately avoids any domain-specific logic (e.g. taxonomic resolution), serving as a thin abstraction over Delta Lake snapshot comparisons.


### Step 9. debug_snapshot.py (System Sanity Check & Debug Harness)

**Responsibility**

Provide a minimal, end-to-end runnable workflow to validate that the core
RefSeq pipeline infrastructure is functioning correctly.

Specifically, this script verifies:

1. Spark + Delta Lake initialization
2. RefSeq assembly index download and parsing
3. FTP-based hash retrieval (annotation / MD5)
4. Hash snapshot DataFrame construction
5. Delta write path correctness and SQL-level readability

**Usage**
python -m refseq_pipeline.core.debug_snapshot


### Step 10. driver.py (Pipeline Orchestration Layer)
**Responsibility**
Provide a high-level orchestration layer that wires together the core RefSeq pipeline modules into a coherent execution flow.

Specifically, this module coordinates:

1. Fetching genome assembly metadata
2. Generating or loading hash snapshots
3. Detecting incremental changes
4. Parsing selected reports into CDM format
5. Writing normalized outputs into Delta Lake
6. Performing post-write cleanup and optimization

**Role in the Architecture**

This module acts as the **glue layer** of the pipeline.

- It does not implement domain logic itself
- It does not fetch data directly from APIs or FTP
- It does not define schemas or parsing rules

Instead, it composes and sequences lower-level modules such as:

- `datasets_api.py`
- `refseq_io.py`
- `hashes_snapshot.py`
- `hashes_diff.py`
- `cdm_parse.py`
- `spark_delta.py`



## Execution Entry Points

The following modules are intended to be executed directly:
1. **driver.py**  
   Primary pipeline entry point. Orchestrates metadata fetch, hash snapshot generation, incremental diffing, CDM parsing, and Delta Lake writes.

2. **debug_snapshot.py**  
   Diagnostic and validation script. Verifies Spark + Delta setup, RefSeq index resolution, FTP hash fetching, snapshot creation, and Delta writes.
   Intended for one-time setup validation and troubleshooting.

All other modules are designed to be imported and composed, not executed directly.

## Incremental Update Workflow

The RefSeq pipeline is incremental by design and avoids full re-ingestion whenever possible.

Typical update flow:

1. Create a new hash snapshot for all (or selected) assemblies
2. Compare the new snapshot against a previous snapshot
3. Identify accessions and taxonomy IDs whose content has changed
4. Fetch metadata only for affected taxa
5. Re-parse selected reports into CDM format
6. Overwrite or merge Delta tables with deduplication

This workflow ensures that only biologically meaningful changes trigger downstream recomputation.

## Results and Benefits

- Orders-of-magnitude faster than full RefSeq re-ingestion
- Deterministic and reproducible updates
- Scales to large taxonomic scopes
- Minimizes unnecessary Spark recomputation

## Command-Line Interface (CLI) Modules

The RefSeq pipeline exposes a small set of CLI-oriented entry points designed for:

- Incremental updates
- Snapshot comparison
- Operational debugging
- Automation (cron, Airflow, CI/CD jobs)

CLI modules act as orchestration layers.

They do **not** implement business logic themselves.
Instead, they coordinate functionality from the `core/` modules, including:
- API access
- Hash snapshot generation
- Snapshot diffing
- CDM parsing
- Delta Lake I/O

This design allows CLI interfaces to remain thin, stable, and easy to evolve independently of the underlying data processing logic.


# RefSeq Pipeline  
## Incremental Update Architecture & Execution Guide

The RefSeq Pipeline is a **Spark + Delta Lake–based incremental ingestion system** for RefSeq genome assemblies.
Instead of reprocessing the full RefSeq dataset on every release, the pipeline:

- Tracks **content-level changes** using hash snapshots
- Detects **incremental differences** between releases
- Reprocesses **only affected genomes**
- Normalizes heterogeneous NCBI responses into a stable **CDM (Common Data Model)**

This document focuses on **architecture, execution order, and operational usage**, rather than low-level implementation details.

---

## Recommended Execution Order (CLI Modules)

The pipeline provides two operational modes:

- **Option A**: Single-command production entry point  
- **Option B**: Step-by-step execution for debugging or fine-grained control  

---

## Option A: Production Entry Point (Recommended)

For routine RefSeq updates, use the unified entry point:
python -m refseq_pipeline.cli.refseq_update_manager


## Option B: Step-by-step Execution

The following steps reflect the **logical dependency chain** between CLI modules.
This mode is recommended for debugging, validation, or fine-grained control.

---

### Step 1. Save RefSeq Assembly Index

**Module:** `save_index_tsv.py`

**Purpose:**  
Download and version the RefSeq assembly index (`assembly_summary_refseq.tsv`).

This index acts as the authoritative source for:
- Assembly accessions
- FTP paths
- Taxonomic relationships

**When to run this step:**
- On **first-time setup**
- When **manually versioning** RefSeq releases
- When validating index-level changes before snapshot creation

### Step 2. Create Hash Snapshot

**Module:** `snapshot_hashes.py`

**Purpose:**  
Generate **content-level hash snapshots** for all RefSeq assemblies.

Each snapshot represents the biological content state of assemblies at a specific point in time, independent of metadata or release labels.

**Output:**
- Delta table: `refseq_api.assembly_hashes`
- Key columns:
  - `accession`
  - `content_sha256`
  - `tag` (snapshot identifier)
  - `retrieved_at`

**Notes:**
- Hash snapshots are used as the **single source of truth** for change detection.
- This step **must be executed before** any snapshot diff or incremental update.

### Step 3. Detect Updated or New Assemblies

**Module:** `detect_updates.py`

**Purpose:**  
Compare two hash snapshot tags and identify assemblies whose **content state has changed**.

This step operates purely on hash snapshots and does **not** fetch any external metadata.

**Output:**
- CSV file containing accession-level differences
- Each record includes:
  - `accession`
  - `change_type` (`new`, `updated`)

### Step 4. Map Accession-Level Changes to TaxIDs

**Module:** `diff_changed_taxids.py`

**Purpose:**  
Translate accession-level content changes into the corresponding **affected Taxonomy IDs (TaxIDs)**.

This step resolves biological impact boundaries by grouping changed assemblies under their parent taxa.

**Output:**
- `changed_taxids.json` — a deduplicated list of TaxIDs whose assemblies were:
  - newly added
  - updated
  - removed

### Step 5. Fetch and Rebuild CDM Tables for Affected TaxIDs

**Module:** `fetch_taxon_reports.py`

**Purpose:**  
Fetch genome assembly metadata from the NCBI Datasets API **only for TaxIDs identified as changed**,  
normalize the responses into the Common Data Model (CDM),  
and write the results into Delta Lake tables.

This step performs the **actual incremental reprocessing**, ensuring that only biologically affected taxa are rebuilt.

**Output:**
- Delta tables containing CDM-normalized assembly records
- Records are:
  - schema-aligned (CDM_SCHEMA)
  - deduplicated by deterministic CDM IDs
  - safe for overwrite or append operations

## Debug and Utility Scripts (Non-Production)

The following scripts are **not part of the production execution path**,  
but are provided to support validation, inspection, and debugging during development and operations:

- **debug_snapshot.py**  
  Validate Spark and Delta Lake setup, RefSeq index loading, FTP hash fetching, and snapshot writes.

- **compare_snapshots.py**  
  Inspect and explore differences between two snapshot tags interactively.

- **debug_register.py**  
  Manually register an existing Delta table into the Spark metastore for inspection or recovery.

These utilities are intentionally separated from the main pipeline to keep production execution deterministic and minimal.

## Summary

The RefSeq Pipeline implements a **deterministic, incremental update architecture** that enables scalable and efficient processing of large genome datasets by:

- Tracking **real content-level changes** via hash snapshots  
- Avoiding unnecessary full re-ingestion and recomputation  
- Isolating responsibilities across clearly defined core and CLI modules  
- Enforcing schema consistency through a Common Data Model (CDM)

### Recommended Operational Entry Point

For routine updates and scheduled runs, use:

```bash
python -m refseq_pipeline.cli.refseq_update_manager
