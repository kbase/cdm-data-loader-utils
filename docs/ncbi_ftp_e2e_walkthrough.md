# NCBI FTP Pipeline — Local End-to-End Walkthrough

Step-by-step instructions for running a full production transfer in the Lakehouse,
or running a small (≤ 10 assembly) local end-to-end sync
of NCBI RefSeq records against a local CEPH container.  The walkthrough uses
a Jupyter notebook for manifest generation, a local or containerized download
and staging CLI tool, and a CLI tool for promotion of the staged records to an S3 store.

> **Prerequisites:**
> - Docker or Podman
> - [uv](https://docs.astral.sh/uv/) (for running notebooks locally)
> - Network access to `ftp.ncbi.nlm.nih.gov`

---

## Architecture overview

```
 Manifest (notebook)
 (to be replaced with CLI)  Download (container CLI)      Promote (local CLI)
┌────────────────────┐     ┌───────────────────────┐     ┌──────────────────────┐
│ Manifest notebook  │     │ ncbi_ftp_sync CLI     │     │ Promote CLI tool     │
│ ─ download FTP     │────▶│ ─ read manifest       │────▶│ ─ promote staged     │
│   assembly summary │     │ ─ parallel FTP DL     │     │   files to Lakehouse │
│ ─ diff against     │     │ ─ MD5 verify          │     │ ─ archive old ver.   │
│   previous         │     │ ─ write .md5 sidecars │     │ ─ trim manifest      │
│ ─ write manifests  │     └──────────┬────────────┘     └──────────────────────┘
└────────────────────┘                │
                                 local volume
                                 mounted into
                                 the container
```

---

### Path formats used

| Format | Example | Description |
|--------|---------|-------------|
| **s3:// URI** | `s3://cdm-lake/staging/run1/` | Full URI with scheme + bucket + key |
| **bucket name** | `cdm-lake` | Just the bucket, no scheme |
| **S3 key prefix** | `tenant-general-warehouse/kbase/datasets/ncbi/` | Path within a bucket (no scheme, no bucket) |
| **S3 object key** | `staging/transfer_manifest.txt` | Single object key within a bucket |
| **local path** | `output/removed_manifest.txt` | Filesystem path on the host |

### Lakehouse object (final location)

```
s3://{LAKEHOUSE_BUCKET}/{LAKEHOUSE_KEY_PREFIX}/raw_data/{GCF|GCA}/{nnn}/{nnn}/{nnn}/{assembly_dir}/{filename}
     └── bucket ─────┘ └── key prefix ───────┘└── build_accession_path() ────────────────────────┘
```

Example:
```
s3://cdm-lake/tenant-general-warehouse/kbase/datasets/ncbi/raw_data/GCF/900/000/615/GCF_900000615.1_PRJEB7657_assembly/GCF_900000615.1_PRJEB7657_assembly_genomic.fna.gz
```

### Staging object (Phase 2 output)

```
s3://{STAGING_BUCKET}/{STAGING_KEY_PREFIX}/raw_data/{GCF|GCA}/{nnn}/{nnn}/{nnn}/{assembly_dir}/{filename}
     └── bucket ─────┘ └── key prefix ────┘└── build_accession_path() ────────────────────────┘
```

### Local output (Phase 1)

```
{OUTPUT_DIR}/transfer_manifest.txt
{OUTPUT_DIR}/removed_manifest.txt
{OUTPUT_DIR}/updated_manifest.txt
{OUTPUT_DIR}/diff_summary.json
```

---

## 1. Setup

### Local testing

### Start CEPH

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
Set CEPH credentials as environment variables:

```sh
export AWS_ENDPOINT_URL=http://localhost:9000
export AWS_ACCESS_KEY_ID=test_access_key
export AWS_SECRET_ACCESS_KEY=test_access_secret
```

(Note that a similar service is included in the `docker-compose` configuration file at the root of
this repository that is used in CI test workflows.)

Create test buckets from the command line using the included `scripts/s3_local.py` helper (requires no extra installs — only `boto3` which is already a project dependency):

```sh
uv run python scripts/s3_local.py mb s3://cdm-lake
uv run python scripts/s3_local.py mb s3://cts
```

### Lakehouse

#### Build `cdm-data-loaders`

First, clone the `cdm-data-loaders` repo in your Lakehouse user space. Then, build the package
in a virtual environment and register it as a Jupyter kernel:
```bash
cd cdm-data-loaders
uv sync
source .venv/bin/activate
uv pip install -e .
uv pip install ipykernel
uv run python -m ipykernel install --user --name cdm-data-loaders --display-name "cdm-data-loaders"
```
Then, when you open the manifest or promote notebooks, choose the `cdm-data-loaders` kernel.


---

## 2. Phase 1 — Generate manifests (notebook)

Open `notebooks/ncbi_ftp_manifest.ipynb` in JupyterLab or VS Code.

### Constants to change (Cell 3)

| Constant              | Walkthrough value                | Format | Why                                                     |
|-----------------------|----------------------------------|--------|---------------------------------------------------------|
| `DATABASE`            | `"refseq"`                       | string | keep as-is                                              |
| `PREFIX_FROM`         | `"900"`                          | string | high-numbered prefix → few assemblies, fast diffing     |
| `PREFIX_TO`           | `"900"`                          | string | single prefix bucket                                    |
| `LIMIT`               | `10`                             | int    | cap to 10 assemblies                                    |
| `PREVIOUS_SUMMARY_URI` | `None`                          | s3:// URI | first run — everything is "new"                       |
| `SNAPSHOT_UPLOAD_URI`  | `None`                          | s3:// URI | skip S3 upload for local testing                      |
| `LAKEHOUSE_BUCKET`    | `"cdm-lake"` (or `None`)         | bucket name | set to prune assemblies already in the Lakehouse   |
| `STORE_KEY_PREFIX`    | `"tenant-general-warehouse/kbase/datasets/ncbi/"` | S3 key prefix | default Lakehouse path prefix    |
| `OUTPUT_DIR`          | `Path("output")`                 | local path | keep as-is (local directory)                        |

### Initialise the S3 client for CEPH

If you set `PREVIOUS_SUMMARY_URI`, `SNAPSHOT_UPLOAD_URI`, `LAKEHOUSE_BUCKET`,
or `STAGING_URI` to point at your local CEPH, you must initialise
the S3 client **before** running the cells that use them.  Insert a new cell
after Cell 1 (Imports) with:

```python
from cdm_data_loaders.utils.file_transfer.s3.client import get_s3_client, reset_s3_client

reset_s3_client()
get_s3_client(
    {
        "endpoint_url": "http://localhost:9000",
        "aws_access_key_id": "test_access_key",
        "aws_secret_access_key": "test_access_secret",
    }
)
```

If all three S3 variables are `None` (purely local testing), this cell can
be skipped — though on repeat runs you should set `LAKEHOUSE_BUCKET` so
assemblies already promoted to the Lakehouse are pruned from the transfer
manifest.

### Optional: Bootstrap from existing store (Cell 5)

If you have a pre-populated S3 store but lack a baseline assembly summary,
you can scan the store to generate a synthetic baseline. This is especially
useful for large stores (100K+ assemblies) where verifying against FTP
checksums would take days.

**When to use this:**
- First run against an existing, pre-populated store
- You want to start diffing without waiting for checksum verification
- You don't have a previous assembly summary snapshot to compare against

**How it works:**
1. Set `SCAN_STORE = True` in Cell 5
2. The notebook scans all objects under `s3://{LAKEHOUSE_BUCKET}/{STORE_KEY_PREFIX}`
3. For each unique assembly found, it extracts the accession and uses the
   earliest object `LastModified` as a conservative `seq_rel_date`
4. It saves the synthetic summary to `LOCAL_SYNTHETIC_SUMMARY` (default:
   `output/synthetic_summary_from_store.txt`)
5. This becomes the baseline for diffing; subsequent runs can load this
   file as `PREVIOUS_SUMMARY_URI`

**Example (for a 500K-assembly store):**
```python
SCAN_STORE = True
LAKEHOUSE_BUCKET = "cdm-lake"
STORE_KEY_PREFIX = "tenant-general-warehouse/kbase/datasets/ncbi/"
LOCAL_SYNTHETIC_SUMMARY = Path("output/synthetic_summary_from_store.txt")

# After running Cell 5, upload the result to S3 for future runs:
# s3 cp output/synthetic_summary_from_store.txt s3://cdm-lake/assembly_summaries/synthetic_base.txt
# Then in future runs, set:
# PREVIOUS_SUMMARY_URI = "s3://cdm-lake/assembly_summaries/synthetic_base.txt"
```

**Performance:** Scanning typically takes 5–10 minutes for 500K assemblies
(vs. ~6 days of checksum verification).

### Run the notebook

Execute all cells in order.  After Cell 7 finishes you should see files in
`output/`:

```
output/
  transfer_manifest.txt   # ≤ 10 FTP directory paths
  removed_manifest.txt    # empty on first run
  updated_manifest.txt    # empty on first run
  diff_summary.json       # counts of new/updated/replaced/suppressed
```

Inspect `transfer_manifest.txt` — each line is an FTP directory path like:

```
/genomes/all/GCF/900/000/615/GCF_900000615.1_PRJEB7657_assembly
```

### Optional: upload manifests to S3 for CTS

Cell 7 optionally uploads the manifests to an S3 staging prefix so that CTS
can stage them into the container.  For local testing, set
`STAGING_URI = None` (the default) and copy the manifest manually in
Step 3b below.

If you are testing against CEPH and want to exercise the S3 upload path:

```python
STAGING_URI = "s3://cts/staging/run1/"
```

> **Tip:** If you re-run later with `PREVIOUS_SUMMARY_URI` pointing at a
> snapshot from a prior run you will see `updated`, `replaced`, and
> `suppressed` entries in the diff.

---

## 3. Phase 2 — Download assemblies (container)

Phase 2 uses the `ncbi_ftp_sync` CLI, which is the container's built-in entry
point for parallel FTP downloads.

> **CTS (CDM Task Service):** In production, Phase 2 runs as a CTS job.
> CTS stages input files from S3 into the container's filesystem mount
> (`/input_dir`) and copies container output back to S3 (`/output_dir`).
> The container itself never receives S3 credentials.
> See [cdm-task-service](https://github.com/kbase/cdm-task-service) for details.

For local testing without a CTS instance we run the container directly with
Docker (or Podman), mounting the manifest produced in Phase 1 as input and a
local staging directory as output.

### 3a. Build the container image

```sh
# From the repository root
docker build -t cdm-data-loaders .
```

### 3b. Prepare local directories

```sh
mkdir -p notebooks/staging/input
cp notebooks/output/* notebooks/staging/input
```

### 3c. Run the download

```sh
docker run --rm \
  --userns=keep-id \
  -v "$(pwd)/notebooks/staging/input:/input:ro" \
  -v "$(pwd)/notebooks/staging:/output" \
  cdm-data-loaders ncbi_ftp_sync \
    --manifest /input/transfer_manifest.txt \
    --output-dir /output \
    --threads 2 \
    --limit 10
```

> **Note:** `--userns=keep-id` maps your host UID into the container so
> bind-mount writes work with Podman's rootless mode.  If you use Docker
> instead, replace it with `--user "$(id -u):$(id -g)"`.

| Flag            | Purpose                                                   |
|-----------------|-----------------------------------------------------------|
| `--manifest`    | Path to the transfer manifest inside the container        |
| `--output-dir`  | Where downloads land (mounted from host `staging/`)       |
| `--threads`     | Parallel FTP connections (2 is polite for testing)        |
| `--limit`       | Redundant safety cap (already limited in Phase 1)         |

After the container exits, `notebooks/staging/` will contain something like:

```
staging/
  raw_data/GCF/900/000/615/GCF_900000615.1_PRJEB7657_assembly/
    GCF_900000615.1_PRJEB7657_assembly_genomic.fna.gz
    GCF_900000615.1_PRJEB7657_assembly_genomic.fna.gz.md5
    GCF_900000615.1_PRJEB7657_assembly_protein.faa.gz
    GCF_900000615.1_PRJEB7657_assembly_protein.faa.gz.md5
    ...
  download_report.json
```

Each data file has a `.md5` sidecar containing the hex digest verified against
the FTP server's `md5checksums.txt`.

> **Without Docker:** You can also run the CLI directly if you have the project
> installed locally:
>
> ```sh
> uv run ncbi_ftp_sync \
>   --manifest notebooks/output/transfer_manifest.txt \
>   --output-dir staging \
>   --threads 2 --limit 10
> ```

### 3d. Local Testing: Upload staged files to CEPH

The download step writes to the local filesystem.  To feed Phase 3 we need
to upload the staged files into CEPH under a staging prefix:

```sh
uv run python scripts/s3_local.py cp notebooks/staging/raw_data/ s3://cts/staging/run1/raw_data/
```

Verify the upload:

```sh
uv run python scripts/s3_local.py ls s3://cts/staging/run1/
```

---

## 4. Phase 3 — Promote & archive (CLI tool)

Phase 3 uses the `ncbi_ftp_promote` CLI tool to promote staged assemblies from
the S3 staging prefix to their final Lakehouse paths, archive replaced or
suppressed assemblies, and trim the transfer manifest for resumability.

### Arguments reference

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--staging-path` | `-s` | *(required)* | S3 key prefix where Phase 2 wrote its output (must contain a `raw_data/` folder) |
| `--destination-path` | | `tenant-general-warehouse/kbase/datasets/ncbi` | S3 key prefix in the destination bucket to promote files into |
| `--staging-bucket` | | `cts` | S3 bucket containing the staged files |
| `--destination-bucket` | | `cdm-lake` | S3 bucket to promote files into (Lakehouse) |
| `--removed-manifest` | `-r` | *(none)* | Local path to the removed manifest from Phase 1; omit to skip archiving removed assemblies |
| `--updated-manifest` | `-u` | *(none)* | Local path to the updated manifest from Phase 1; omit to skip archiving updated assemblies |
| `--transfer-manifest` | `-t` | `{staging-path}/transfer_manifest.txt` | S3 key of the transfer manifest to trim after a successful promote |
| `--dry-run` | | `False` | Log what would happen without making any changes |



First, do a dry run to make sure everything looks as expected:

```sh
uv run ncbi_ftp_promote \
  --staging-path staging/run1 \
  --removed-manifest notebooks/output/removed.txt \
  --updated-manifest notebooks/output/updated.txt \
  --dry-run
```

Once everything looks good, run the actual promotion:
```sh
uv run ncbi_ftp_promote \
  --staging-path staging/run1 \
  --removed-manifest notebooks/output/removed.txt \
  --updated-manifest notebooks/output/updated.txt
```

If you exclude the removed and updated manifests, no archiving will occur, just promotion of staged records.

The CLI prints a promote summary on completion:

```
PROMOTE SUMMARY: 10 promoted, 0 archived, 0 failed
```

After promotion the final Lakehouse layout in CEPH will look like:

```
cdm-lake/
  tenant-general-warehouse/kbase/datasets/ncbi/
    raw_data/GCF/900/000/615/GCF_900000615.1_.../
      GCF_900000615.1_..._genomic.fna.gz          (with md5 in user metadata)
      GCF_900000615.1_..._protein.faa.gz
      ...
```

---

## 5. Inspect results in CEPH

Use the CLI to inspect the final state of the store:

```sh
# List final Lakehouse objects
uv run python scripts/s3_local.py ls \
  s3://cdm-lake/tenant-general-warehouse/kbase/datasets/ncbi/raw_data/

# Check user metadata (md5) on a specific object
uv run python scripts/s3_local.py head \
  s3://cdm-lake/tenant-general-warehouse/kbase/datasets/ncbi/raw_data/GCF/900/000/615/GCF_900000615.1_PRJEB7657_assembly/GCF_900000615.1_PRJEB7657_assembly_genomic.fna.gz
```

### Frictionless metadata descriptors

Each promoted assembly gets a [frictionless](https://framework.frictionlessdata.io/) data package descriptor stored at:

```
s3://{LAKEHOUSE_BUCKET}/{LAKEHOUSE_KEY_PREFIX}metadata/{assembly_dir}_datapackage.json
```

For example:

```
s3://cdm-lake/tenant-general-warehouse/kbase/datasets/ncbi/metadata/GCF_900000615.1_PRJEB7657_assembly_datapackage.json
```

The descriptor follows the KBase credit metadata schema (v1.0) and records:

- **identifier** — `NCBI:{accession}`, e.g. `NCBI:GCF_900000615.1`
- **resource_type** — always `"dataset"`
- **resources** — list of promoted files with their final S3 key, byte size,
  file format, and MD5 hash (when available)
- **contributors / publisher** — NCBI organizational metadata
- **meta.saved_by** — `"cdm-data-loaders-ncbi-ftp"`

When an assembly is archived (updated or removed), its live descriptor is
copied to:

```
s3://{LAKEHOUSE_BUCKET}/{LAKEHOUSE_KEY_PREFIX}archive/{release_tag}/metadata/{assembly_dir}_datapackage.json
```

Use `scripts/s3_local.py ls` to list all descriptors written in a promote run:

```sh
uv run python scripts/s3_local.py ls \
  s3://cdm-lake/tenant-general-warehouse/kbase/datasets/ncbi/metadata/
```

To inspect a descriptor directly:

```sh
uv run python scripts/s3_local.py cat \
  s3://cdm-lake/tenant-general-warehouse/kbase/datasets/ncbi/metadata/GCF_900000615.1_PRJEB7657_assembly_datapackage.json
```

---

## 6. Incremental run (second sync)

To exercise the diff/update/archive logic, repeat the pipeline with a
previous snapshot:

2. **Phase 1:** Set `PREVIOUS_SUMMARY_URI` to an S3 path where you upload the
   raw summary from the first run, or save the `raw_summary` string from Cell 4
   to a local file and pass it via `parse_assembly_summary(Path("prev.txt"))`.
2. **Phase 1:** The diff will now show `updated`, `replaced`, and
   `suppressed` entries (if any changed between runs).
3. **Phase 2:** Download the new manifest.
4. **Phase 3:** Pass `--removed-manifest` and `--updated-manifest` pointing at
   the files produced in Phase 1.  Updated assemblies are archived before
   overwrite; removed assemblies are archived and deleted.

---

## 7. Cleanup

```sh
# Stop and remove CEPH
docker stop ceph && docker rm ceph

# Remove local staging data
rm -rf staging/ output/
```

---

## Troubleshooting

| Symptom | Cause | Fix |
|---------|-------|-----|
| `connect_ftp() timeout` | NCBI FTP may be slow or rate-limited | Retry; reduce `--threads` to 1 |
| Phase 3 shows 0 promoted | Staging prefix doesn't match or bucket is wrong | Verify `--staging-path` matches the S3 upload path from Step 3d |
| Phase 3 S3 auth error | Missing credentials for CEPH | Export `AWS_ENDPOINT_URL`, `AWS_ACCESS_KEY_ID`, and `AWS_SECRET_ACCESS_KEY` before running |
| Container can't reach FTP | Docker network isolation | Use `--network host` or ensure DNS resolution works inside the container |

---

## Reference: file filters

Phase 2 downloads only files matching these suffixes (defined in
`cdm_data_loaders.ncbi_ftp.assembly.FILE_FILTERS`):

| Suffix | Content |
|--------|---------|
| `_genomic.fna.gz` | Genome nucleotide sequences |
| `_genomic.gff.gz` | Genome annotations (GFF3) |
| `_protein.faa.gz` | Protein sequences |
| `_gene_ontology.gaf.gz` | GO annotations |
| `_assembly_report.txt` | Assembly metadata |
| `_assembly_stats.txt` | Assembly statistics |
| `_assembly_regions.txt` | Assembly regions |
| `_ani_contam_ranges.tsv` | ANI contamination ranges |
| `_gene_expression_counts.txt.gz` | Gene expression counts |
| `_normalized_gene_expression_counts.txt.gz` | Normalised expression counts |

Plus the per-assembly `md5checksums.txt` which is always downloaded for
integrity verification.
