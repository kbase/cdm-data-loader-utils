# PDB Pipelines - End-to-End Transfer Walkthrough

## Workflow Overview

Semi-automated PDB transfers involve three distinct pipelines run sequentially:
* `manifest-generation`: Triggered manually in the Lakehouse.
    - Downloads current PDB holdings files
    - Optionally bootstraps current Lakehouse PDB snapshot (if previous snapshot was not saved)
    - Compares current PDB holdings to Lakehouse PDB snapshot to generate:
        * `transfer_manifest.txt`: New and updated PDB records to download
        * `updated_manifest.txt`: Lakehouse records to archive before promoting new downloads
        * `removed_manifest.txt`: Lakehouse records to archive (will not be replaced with new download)
        * `missing_dates.txt`: Current PDB holdings for which a last-modified date was not included. (Manual review)
* `download`: cdm-task-service containerized pipeline.
    - Downloads records in `transfer_manifest.txt` to staging folder
    - (possibly verifies file contents against JSONSchema - TBD)
* `promote`: Triggered manually in the Lakehouse.
    - Archives PDB records in `updated_manifest.txt` and `removed_manifest.txt`
    - Promotes staged new downloads to PDB Lakehouse destination
    - Writes frictionless file descriptors for newly promoted PDB records

### Note for local testing

It is possible to test the end-to-end workflow with a local containerized CEPH test store.
To start the container:
```sh
docker run -d \
  --name ceph \
  -p 9000:8080 \
  -p 9001:8443 \
  -e RGW_PORT=8080 \
  -e RGW_ACCESS_KEY=test_access_key \
  -e RGW_SECRET_KEY=test_access_secret \
  ghcr.io/kbasetest/ceph-rgw-test-image:0.1.5
```
Then, set CEPH credentials as environment variables:

```sh
export AWS_ENDPOINT_URL=http://localhost:9000
export AWS_ACCESS_KEY_ID=test_access_key
export AWS_SECRET_ACCESS_KEY=test_access_secret
```

Finally, create the standard staging and destination buckets in the test store:

```sh
uv run python scripts/s3_local.py mb s3://cdm-lake
uv run python scripts/s3_local.py mb s3://cts
```

### Package Build

The `cdm-data-loaders` package must be present locally. If you haven't already build it, follow these steps:

```sh
cd cdm-data-loaders
uv sync
source .venv/bin/activate
uv pip install -e .
```

## Manifest Generation

Manifest generation is triggered by calling the `pdb_manifest` cli tool. At minimum, you must provide a path
to the folder the manifests should be written to. If that folder doesn't exist, it will be created:

```sh
uv run pdb_manifest --output-path out
```

This will trigger the manifest generation and save the manifest files to `./out/`.

### `pdb_manifest` CLI Tool options

* `--help`: Show full usage (includes some options that are ignored by `pdb_manifest`).
* `--output-path {path}`: **(Required)** Path to the local folder where manifest files are saved. Created automatically if it does not exist.
* `--destination-bucket {str}`: Specify a non-standard S3 bucket for the current Lakehouse (or test store) records.
* `--destination-prefix {str}`: Specify a non-standard S3 key prefix for the current Lakehouse (or test store) records. Must contain `raw_data/`.
* `--snapshot {str}`: Specify a non-standard path (relative to the `destination-prefix`) for the current snapshot file.
* `--bootstrap {date}`: Bootstrap a snapshot file using the current Lakehouse (or test store) state, assigning the specified ISO 8601 date as the last-modified date for all records (e.g. `2024-01-01`).
* `--skip-diff True`: Disable reading of the snapshot file. All current records in the PDB holdings file will be included in the `transfer_manifest.txt`.
* `--regex {str}` / `-r {str}`: Optionally provide a regex filter expression. Only PDB record ids that pass the filter will be included in manifest files.

## Download

_Coming Soon!_

## Promote

_Coming Soon!_
