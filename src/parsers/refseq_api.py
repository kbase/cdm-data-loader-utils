import os
import uuid
from datetime import datetime, date
import re
import click
import pandas as pd
import requests
from typing import Optional
from typing import Any, Literal
from pyspark.sql.types import StructType, StructField, StringType

from pyspark.sql import SparkSession
from delta import configure_spark_with_delta_pip


"""

python refseq_api.py \
  --taxid "224325, 2741724, 193567" \
  --database refseq_api \
  --mode overwrite \
  --debug \
  --unique-per-taxon

"""


# ---------------- Spark + Delta ----------------

def build_spark(database: str) -> SparkSession:
    """
    Initialize a Spark session with Delta Lake support and create the specified database if it doesn't exist.
    """
    builder = (
        SparkSession.builder.appName("NCBI Datasets -> CDM")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )
    spark = configure_spark_with_delta_pip(builder).getOrCreate()

    # Create the database namespace if it doesn't already exist
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {database}")
    return spark


def write_delta(
        spark: SparkSession,
        pandas_df: pd.DataFrame,
        database: str,
        table: str,
        mode: str = "append") -> None:
    """
    Write Pandas DataFrame to a Delta Lake table using Spark. 
    Supports append or overwrite.
    Special handling for 'contig_collection' schema to avoid inference error.
    """

    if pandas_df is None or pandas_df.empty:
        print(f"No data to write to {database}.{table}")
        return

    print(f"Writing {table} with {len(pandas_df)} rows")
    print(pandas_df.dtypes)
    print(pandas_df.head(10))

    # -------- Special schema for problematic table --------
    if table == "contig_collection":
        schema = StructType([
            StructField("collection_id", StringType(), True),
            StructField("contig_collection_type", StringType(), True),
            StructField("ncbi_taxon_id", StringType(), True),
            StructField("gtdb_taxon_id", StringType(), True),
        ])
        pandas_df = pandas_df.astype(str).where(pandas_df.notnull(), None)
        spark_df = spark.createDataFrame(pandas_df, schema=schema)
    else:
        spark_df = spark.createDataFrame(pandas_df)

    # -------- Writer with proper options --------
    writer = spark_df.write.format("delta").mode(mode)
    if mode == "append":
        writer = writer.option("mergeSchema", "true")
    elif mode == "overwrite":
        writer = writer.option("overwriteSchema", "true")

    writer.saveAsTable(f"{database}.{table}")
    print(f"Saved {len(pandas_df)} rows to {database}.{table} (mode={mode})")


def preview_or_skip(
        spark: SparkSession,
        database: str,
        table: str,
        limit: int = 20) -> None:
    
    """
    Preview the first N rows of a Delta table if it exists.
    """

    full_table = f"{database}.{table}"
    if spark.catalog.tableExists(full_table):
        print(f"Showing the first {limit} row from {full_table}:")
        spark.sql(f"SELECT * FROM {full_table} LIMIT {limit}").show(truncate=False)
    else:
        print(f"Table {full_table} not found. Skipping preview.")


# ---------------- Helpers ----------------

def parse_taxid_args(taxid_arg: Optional[str], taxid_file: Optional[str]) -> list[str]:
    """
    Parse and collect valid numeric TaxIDs from command-line arguments and file.
    """

    ## empty list to collect taxids，avoid the duplicate TaxIDs 
    taxids: list[str] = []

    # Parse taxid argument: --taxid "224325, 2741724" 
    if taxid_arg:
        id_list = taxid_arg.split(",") ## separate them into a list using commas
        for num in id_list:
            # Keep only digits
            id = re.sub(r"\D+", "", num.strip()) ## Remove all non-numeric characters
            if id:
                taxids.append(id)

    # Parse --taxid-file 
    if taxid_file:
        if not os.path.exists(taxid_file):
            raise click.BadParameter(f"Path '{taxid_file}' does not exist.", param_hint="--taxid-file")
        with open(taxid_file, "r", encoding="utf-8") as f:
            for line in f:
                id = re.sub(r"\D+", "", line.strip())
                if id:
                    taxids.append(id)

    # Deduplicate while preserving order
    seen = set()
    unique_taxids = []
    for id in taxids:
        if id not in seen:
            seen.add(id)
            unique_taxids.append(id)

    return unique_taxids


# ---------------- NCBI Datasets v2 ----------------

def fetch_reports_by_taxon(
    taxon: str,
    api_key: str | None = None,
    page_size: int = 500,
    refseq_only: bool = True,
    current_only: bool = True,
    debug: bool = False,
):
    """
    Generator to iterate through genome dataset reports from NCBI Datasets v2 API by TaxID.
    
    Features:
    - Calls the NCBI Datasets v2 endpoint for genome reports.
    - Applies filters: RefSeq only / current assemblies only.
    - Handles pagination via `next_page_token`.
    - Yields report dicts for each assembly.

    """

    # ---------------- API endpoint ----------------
    # Base URL for NCBI Datasets REST API v2
    base = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    url = f"{base}/genome/taxon/{taxon}/dataset_report"

    # ---------------- Request params ----------------
    # metadata + assembly report text
    params = {
        "page_size": page_size,
        "returned_content": "COMPLETE",
        "filters.report_type": "assembly_report"
    }
   
    if current_only:
        params["filters.assembly_version"] = "current"
    if refseq_only:
        params["filters.assembly_source"] = "refseq"

    # ---------------- Headers ----------------
    headers = {"Accept": "application/json"}
    if api_key:
        headers["api-key"] = api_key

    # ---------------- Pagination loop ----------------
    token = None
    while True:
        if token:
            params["page_token"] = token

        # ---- request ----
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=60)
            resp.raise_for_status()
            payload = resp.json()
        except (requests.RequestException, ValueError) as e:
            print(f"Request failed for taxon {taxon}: {e}")
            break

        # ---- Extract reports ----
        reports = payload.get("reports", [])
        if not reports:
            print(f"No reports returned for taxon {taxon}")
            break

        # ---------------- Filter loop ----------------
        for rep in reports:
            info = rep.get("assemblyInfo") or rep.get("assembly_info") or {}
            src_db = info.get("sourceDatabase")

            # Skip if explicitly marked as GenBank 
            if src_db and src_db != "SOURCE_DATABASE_REFSEQ":
                continue

            # Print source info if debugging
            if debug:
                if src_db is None:
                    print(f"[DEBUG] accession={rep.get('accession')} has no sourceDatabase field")
                else:
                    print(f"[DEBUG] accession={rep.get('accession')} sourceDatabase={src_db}")

                # Print the first 200 chars of assemblyReport for inspection
                if info.get("assemblyReport"):
                    snippet = info["assemblyReport"][:200].replace("\n", " ")
                    print(f"[DEBUG] {snippet}")

            # Yield one assembly report to caller
            yield rep

        # ---------------- Handle pagination ----------------
        token = payload.get("next_page_token")
        if not token:
            break



# ---------------- Robust extractors set up ----------------

# regex patterns
PAT_BIOSAMPLE = re.compile(r"\bSAMN\d+\b")
PAT_BIOPROJECT = re.compile(r"\bPRJNA\d+\b")
PAT_GCF = re.compile(r"\bGCF_\d{9}\.\d+\b")
PAT_GCA = re.compile(r"\bGCA_\d{9}\.\d+\b")


def _coalesce(*vals: Any) -> str | None:
    """
    Return the first non-empty, non-whitespace string from a list of inputs.
    """
    for v in vals:
        if isinstance(v, str):
            trimmed = v.strip()
            if trimmed:
                return trimmed
    return None


def _deep_find_str(obj: Any, target_keys: set[str]) -> str | None:
    """
    Recursively search a nested dict/list structure for the first non-empty string value under any of the target_keys.
    _deep_find_str({"assemblyDate": "2000-12-01"}, {"assemblyDate"}) -> "2000-12-01"

    """
    
    ## If the object is a dictionary
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str) and k in target_keys and isinstance(v, str) and v.strip():
                return v.strip()
            
            ## recursively search the v
            found = _deep_find_str(v, target_keys)
            if found:
                return found
            
    ## If the object is a list, search each element
    elif isinstance(obj, list):
        for it in obj:
            found = _deep_find_str(it, target_keys)
            if found:
                return found
    return None


def _deep_collect_regex(obj: Any, pattern: re.Pattern) -> list[str]:
    """
    Recursively collect all unique regex matches from a nested structure

    Args:
        obj: The nested object to search (can be dict, list, or string).
        pattern: Compiled regex pattern to search for within string values.

    Returns:
        A sorted list of unique regex matches found anywhere inside the object.

    """

    results = set()  # Use a set to avoid duplicate matches

    def _walk(x):
        if isinstance(x, dict):
            # Recursively process each value in the dictionary
            for v in x.values():
                _walk(v)
        elif isinstance(x, list):
            # Recursively process each element in the list
            for v in x:
                _walk(v)
        elif isinstance(x, str):
            # Apply regex to the string and add matches to results
            for m in pattern.findall(x):
                results.add(m)

    # Start recursion
    _walk(obj)

    # Convert set to sorted list for consistent ordering
    return sorted(results)


# ---------------- Robust extractors ----------------

def extract_created_date(rep: dict[str, Any], allow_genbank_date: bool = False, debug: bool = False) -> str | None:
    """

    Extract creation/release date for a genome assembly.

    Priority:
    - For RefSeq: release_date > assembly_date > submission_date
    - For GenBank (if allowed): submission_date
    Returns None if no valid date is found.

    """

    # Normalize assembly info
    assem_data = rep.get("assembly_info") or rep.get("assemblyInfo") or {}
    src_db = rep.get("source_database") or assem_data.get("sourceDatabase")

    # Collect candidate dates
    candidates: dict[str, str] = {}
    for src in (assem_data, rep.get("assembly") or {}, rep):
        for key in ["releaseDate", "assemblyDate", "submissionDate",
                    "release_date", "assembly_date", "submission_date"]:
            v = src.get(key) # safely fetch the value (None if key not found)

            # Only accept non-empty string values
            if isinstance(v, str) and v.strip():
                # Normalize the key: "release_date" -> "releasedate"
                norm_key = key.lower().replace("_", "")

                # Store the cleaned value in candidates under the normalized key
                candidates[norm_key] = v.strip()

    if debug and candidates:
        print(f"[DEBUG] found candidates={candidates}, source={src_db}")

    # RefSeq: prioritize release > assembly > submission
    if src_db == "SOURCE_DATABASE_REFSEQ":
        for pref in ("releasedate", "assemblydate", "submissiondate"):
            if pref in candidates:
                return candidates[pref]

    # GenBank: fallback only submission date
    if allow_genbank_date and src_db == "SOURCE_DATABASE_GENBANK":
        if "submissiondate" in candidates:
            return candidates["submissiondate"]

    return None



def extract_assembly_name(rep: dict[str, Any]) -> str | None:
    """
    Extract the assembly name from a genome report record.
    This function normalizes extraction because the assembly name may appear in different parts of the JSON depending on API. 

    """

    assembly_info = rep.get("assemblyInfo") or {}
    a = rep.get("assembly") or {}

    # Try to extract directly from the most common locations
    v = _coalesce(
        assembly_info.get("assemblyName"),
        a.get("assemblyName"),
        rep.get("assemblyName"),
    )
    if v:
        return v

    # Fallback: recursively search the nested structure
    return _deep_find_str(rep, {"assemblyName", "assembly_name"})



def extract_organism_name(rep: dict[str, Any]) -> str | None:
    assembly_info = rep.get("assemblyInfo") or {}
    a  = rep.get("assembly") or {}
    org_top = rep.get("organism") or {}
    
    # Candidate locations where the organism name might appear
    candidates = [
        org_top.get("organismName"),
        org_top.get("scientificName"),
        org_top.get("taxName"),
        (assembly_info.get("organism") or {}).get("organismName") if isinstance(assembly_info.get("organism"), dict) else None,
        (a.get("organism")  or {}).get("organismName") if isinstance(a.get("organism"), dict) else None
        ]
    
    # Return the first non-empty candidate string
    for v in candidates:
        if isinstance(v, str) and v.strip():
            return v.strip()
    
    # Fallback deep search across nested structure
    return _deep_find_str(rep, 
                          {"organismName", "scientificName", "sciName", "taxName", 
                                "displayName", "organism_name"})



def extract_taxid(rep: dict[str, Any]) -> str | None:
    """
    Extract the NCBI Taxonomy ID (taxid) from a genome assembly report.

    """

    def is_numeric_id(value) -> bool:
        """Check if a value is a valid numeric taxid."""
        return isinstance(value, (int, float)) or (isinstance(value, str) and value.strip().isdigit())

    # --- Try top-level organism block ---
    org_top = rep.get("organism") or {}
    v = org_top.get("taxId") or org_top.get("taxid") or org_top.get("taxID")
    if is_numeric_id(v):
        return str(int(v))

    # --- Recursive search ---
    def _deep_find_taxid(x):
        if isinstance(x, dict):
            for k, v in x.items():
                # normalize key (lowercase) and check if it contains 'taxid'
                if isinstance(k, str) and k.lower().replace("_", "") in {"taxid", "taxidvalue"}:
                    if is_numeric_id(v):
                        return str(int(v))
                found = _deep_find_taxid(v)
                if found:
                    return found
        elif isinstance(x, list):
            for it in x:
                found = _deep_find_taxid(it)
                if found:
                    return found
        return None

    return _deep_find_taxid(rep)



def extract_biosample_ids(rep: dict[str, Any]) -> list[str]:
    """
    Extract BioSample IDs from a genome assembly report.
    """
    accs = set()

    # --- Check standard paths ---
    for path in [
        rep.get("assemblyInfo", {}).get("biosample"),
        rep.get("assembly", {}).get("biosample"),
        rep.get("biosample")
    ]:
        if isinstance(path, dict):
            v = path.get("accession") or path.get("biosampleAccession")
            if isinstance(v, str) and v.strip():
                accs.add(v.strip())
        elif isinstance(path, list):
            for it in path:
                if isinstance(it, dict):
                    v = it.get("accession") or it.get("biosampleAccession")
                    if isinstance(v, str) and v.strip():
                        accs.add(v.strip())

    # --- Regex fallback if no BioSamples found ---
    if not accs:
        accs.update(_deep_collect_regex(rep, PAT_BIOSAMPLE))

    # --- Return unique + sorted IDs ---
    return sorted(accs)



def extract_bioproject_ids(rep: dict[str, Any]) -> list[str]:
    accs = set()
    for path in [
        rep.get("assemblyInfo", {}).get("bioproject"),
        rep.get("assembly", {}).get("bioproject"),
        rep.get("bioproject"),
    ]:
        if isinstance(path, dict):
            v = path.get("accession") or path.get("bioprojectAccession")
            if isinstance(v, str) and v.strip():
                accs.add(v.strip())
        elif isinstance(path, list):
            for it in path:
                if isinstance(it, dict):
                    v = it.get("accession") or it.get("bioprojectAccession")
                    if isinstance(v, str) and v.strip():
                        accs.add(v.strip())
    if not accs:
        accs.update(_deep_collect_regex(rep, PAT_BIOPROJECT))
    return sorted(accs)


def extract_assembly_accessions(rep: dict[str, Any]) -> tuple[list[str], list[str]]:
    """
    Extract RefSeq (GCF_) and GenBank (GCA_) accession IDs from an assembly report.
    - Consider only the top-level accession and assembly_info.paired_assembly.
    - Classify accessions into GCF (RefSeq) and GCA (GenBank).
    - Return sorted unique lists for both GCF and GCA.

    Returns:
        A tuple of two lists:
            - List of GCF accessions (RefSeq)
            - List of GCA accessions (GenBank)
    """

    gcf, gca = set(), set()

    def _add_if_valid(acc: str):
        """Helper: classify accession into GCF or GCA bucket."""
        if not isinstance(acc, str) or not acc.strip():
            return
        acc = acc.strip()
        if acc.startswith("GCF_"):
            gcf.add(acc)
        elif acc.startswith("GCA_"):
            gca.add(acc)

    # --- Top-level accession ---
    _add_if_valid(rep.get("accession"))

    # --- Paired assembly accession (from assembly_info) ---
    ai = rep.get("assembly_info") or rep.get("assemblyInfo") or {}
    paired = ai.get("paired_assembly", {})
    if isinstance(paired, dict):
        _add_if_valid(paired.get("accession"))

    return sorted(gcf), sorted(gca)



# ---------------- CDM Builders ----------------
def build_cdm_datasource() -> pd.DataFrame:
    """Generate a single CDM datasource record for NCBI RefSeq (via API)."""
    record = {
        "name": "RefSeq",
        "source": "NCBI RefSeq",
        "url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/",  
        "accessed": date.today().isoformat(),
        "version": "231", 
    }
    return pd.DataFrame([record])


## entity 
CDM_NAMESPACE = uuid.UUID("11111111-2222-3333-4444-555555555555")

def build_entity_id(key: str) -> str:
    """
    Generate a deterministic CDM ID from a given key using UUIDv5.
    """
    k = (key or "").strip()
    return f"CDM:{uuid.uuid5(CDM_NAMESPACE, k)}"


def build_cdm_entity(
    key_for_uuid: str,
    created_date: str | None,
    *,
    entity_type: Literal["contig_collection", "genome", "protein", "gene"] = "contig_collection",
    data_source: str = "RefSeq") -> tuple[pd.DataFrame, str]:

    """
    Build a CDM 'entity' row with CDM-style entity_id, creation, and update timestamps.
    
    Parameters:
    - key_for_uuid: a unique identifier string to derive the CDM UUID
    - created_date: original creation date (from source), or None to use today
    - entity_type: CDM entity type (default: 'contig_collection')
    - data_source: Source name for traceability (default: 'RefSeq')
    
    Returns:
    - A tuple of (DataFrame with 1 row, entity_id string)
    """

    entity_id = build_entity_id(key_for_uuid)
    record = {
        "entity_id": entity_id,
        "entity_type": entity_type,
        "data_source": data_source,
        "created": created_date or date.today().isoformat(),
        "updated": datetime.now().isoformat(timespec="seconds"),
    }
    return pd.DataFrame([record]), entity_id


## contig_collection
def build_cdm_contig_collection(entity_id: str, taxid: str | None = None, collection_type: str = "isolate") -> pd.DataFrame:
    """
    Build the contig_collection CDM table with taxon ID and collection type.
    """
    return pd.DataFrame([{
        "collection_id": str(entity_id),
        "contig_collection_type": str(collection_type),
        "ncbi_taxon_id": f"NCBITaxon:{taxid}" if taxid else None,
        "gtdb_taxon_id": None,  # reserved for future GTDB support
    }])


def build_cdm_name_rows(entity_id: str, rep: dict[str, Any]) -> pd.DataFrame:
    """
    Build two name rows: organism name and assembly name.
    - Organism name comes from extract_organism_name()
    - Assembly name comes from extract_assembly_name()
    """
    rows = []

    # ==== organism name ====
    full_org_name = extract_organism_name(rep)
    if isinstance(full_org_name, str) and full_org_name.strip():
        rows.append({
            "entity_id": str(entity_id),
            "name": full_org_name.strip(),
            "description": "RefSeq organism name",
            "source": "RefSeq"
        })

    # ==== assembly name ====
    asm_name = extract_assembly_name(rep)
    if isinstance(asm_name, str) and asm_name.strip():
        rows.append({
            "entity_id": str(entity_id),
            "name": asm_name.strip(),
            "description": "RefSeq assembly name",
            "source": "RefSeq"
        })

    return pd.DataFrame(rows)


# ---------------- Identifier Constants ----------------
IDENTIFIER_PREFIXES = {
    "biosample": ("Biosample", "BioSample ID"),
    "bioproject": ("BioProject", "BioProject ID"),
    "taxon": ("NCBITaxon", "NCBI Taxon ID"),
    "gcf": ("ncbi.assembly", "NCBI Assembly ID"),
    "gca": ("insdc.gca", "GenBank Assembly ID"),
    "gcf_unit": ("insdc.gca", "RefSeq Unit Assembly ID"),
    "gca_unit": ("insdc.gca", "GenBank Unit Assembly ID"),
}



def build_cdm_identifier_rows(entity_id: str, rep: dict, request_taxid: str | None) -> list[dict]:
    """
    Build CDM 'identifier' rows from parsed NCBI RefSeq metadata.
    """

    rows = []

    # ---- BioSample IDs ----
    for bs in extract_biosample_ids(rep):
        rows.append({
            "entity_id": entity_id,
            "identifier": f"{IDENTIFIER_PREFIXES['biosample'][0]}:{bs.strip()}",
            "source": "RefSeq",
            "description": IDENTIFIER_PREFIXES['biosample'][1],
        })

    # ---- BioProject IDs ----
    for bp in extract_bioproject_ids(rep):
        rows.append({
            "entity_id": entity_id,
            "identifier": f"{IDENTIFIER_PREFIXES['bioproject'][0]}:{bp.strip()}",
            "source": "RefSeq",
            "description": IDENTIFIER_PREFIXES['bioproject'][1],
        })

    # ---- Taxon ID ----
    tx = extract_taxid(rep) or (str(request_taxid).strip() if request_taxid else None)
    if tx and tx.isdigit():
        rows.append({
            "entity_id": entity_id,
            "identifier": f"{IDENTIFIER_PREFIXES['taxon'][0]}:{tx}",
            "source": "RefSeq",
            "description": IDENTIFIER_PREFIXES['taxon'][1],
        })

    # ---- Assembly Accessions (GCF / GCA) ----
    gcf_list, gca_list = extract_assembly_accessions(rep)
    for gcf in gcf_list:
        rows.append({
            "entity_id": entity_id,
            "identifier": f"{IDENTIFIER_PREFIXES['gcf'][0]}:{gcf.strip()}",
            "source": "RefSeq",
            "description": IDENTIFIER_PREFIXES['gcf'][1],
        })
    for gca in gca_list:
        rows.append({
            "entity_id": entity_id,
            "identifier": f"{IDENTIFIER_PREFIXES['gca'][0]}:{gca.strip()}",
            "source": "RefSeq",
            "description": IDENTIFIER_PREFIXES['gca'][1],
        })

    # ---- Deduplicate ----
    seen = set()
    uniq = []
    for r in rows:
        ident = r["identifier"]

        k = (r["entity_id"], ident)
        if k not in seen:
            seen.add(k)
            uniq.append(r)

    return uniq


# ---------------- CLI ----------------
@click.command()
@click.option("--taxid", required=True,
              help="Comma-separated NCBI TaxIDs, e.g. '224325,2741724'.")
@click.option("--api-key", default=None,
              help="Optional NCBI API key (increases rate limits).")
@click.option("--database", default="refseq_api", show_default=True,
              help="Delta schema/database.")
@click.option("--mode", default="overwrite", type=click.Choice(["overwrite", "append"]), show_default=True,
              help="Write mode for Delta tables.")
@click.option("--debug/--no-debug", default=False, show_default=True,
              help="Print per-record parsed fields for troubleshooting.")
@click.option("--allow-genbank-date/--no-allow-genbank-date", default=False, show_default=True,
              help="Allow using GenBank submissionDate as fallback for RefSeq created date.")
@click.option("--unique-per-taxon/--all-assemblies", default=False, show_default=True,
              help="Keep only one assembly per taxon (latest by release_date).")


def cli(taxid, api_key, database, mode, debug, allow_genbank_date, unique_per_taxon):
    main(
        taxid=taxid,
        api_key=api_key,
        database=database,
        mode=mode,
        debug=debug,
        allow_genbank_date=allow_genbank_date,
        unique_per_taxon=unique_per_taxon
    )


def process_report(rep: dict, tx: str, seen: set, debug: bool, allow_genbank_date: bool):
    """
    Process a single assembly report and return partial CDM tables (entity, collection, names, identifiers).
    Skip duplicates using 'seen'.
    """
    entities, collections, names, identifiers = [], [], [], []

    # === accession (prefer GCF, fallback to GCA) ===
    gcf_list, gca_list = extract_assembly_accessions(rep)
    acc = gcf_list[0] if gcf_list else (gca_list[0] if gca_list else None)

    # === assembly + organism names ===
    asm_name = extract_assembly_name(rep)
    org_name = extract_organism_name(rep)

    # === creation date ===
    created = extract_created_date(rep, allow_genbank_date=allow_genbank_date)
    if not created:
        if debug:
            print(f"[WARN] No RefSeq date for accession={acc}")
        created = date.today().isoformat()

    # === unique key for entity ===
    key = _coalesce(acc, asm_name, org_name, tx)
    if not key or key in seen:
        return entities, collections, names, identifiers
    seen.add(key)

    # ---------------- CDM tables ----------------
    df_entity, entity_id = build_cdm_entity(key, created)
    entities.append(df_entity)

    collections.append(build_cdm_contig_collection(entity_id, taxid=tx))

    rows_name = build_cdm_name_rows(entity_id, rep)
    if not rows_name.empty:
        names.extend(rows_name.to_dict(orient="records"))

    rows_id = build_cdm_identifier_rows(entity_id, rep, tx)
    identifiers.extend(rows_id)

    return entities, collections, names, identifiers


def process_taxon(tx: str, api_key: str, debug: bool, allow_genbank_date: bool, unique_per_taxon: bool, seen: set):
    """
    Process all reports for a given TaxID and return combined partial tables.
    """
    entities, collections, names, identifiers = [], [], [], []

    reports = list(fetch_reports_by_taxon(taxon=tx, api_key=api_key))

    # If only unique assembly per taxon → keep latest
    if unique_per_taxon and reports:
        reports.sort(key=lambda r: extract_created_date(r, allow_genbank_date) or "0000-00-00", reverse=True)
        reports = [reports[0]]

    for rep in reports:
        e, c, n, i = process_report(rep, tx, seen, debug, allow_genbank_date)
        entities.extend(e)
        collections.extend(c)
        names.extend(n)
        identifiers.extend(i)

    return entities, collections, names, identifiers



def finalize_tables(entities, collections, names, identifiers):
    """
    Concatenate and deduplicate CDM tables.
    """
    def _concat(frames): return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    def _dedup(pdf, cols): return pdf.drop_duplicates(subset=cols) if not pdf.empty else pdf

    pdf_entity = _dedup(_concat(entities), ["entity_id"])
    pdf_coll   = _dedup(_concat(collections), ["collection_id"])
    pdf_name   = _dedup(pd.DataFrame(names), ["entity_id", "name"]) if names else pd.DataFrame()
    pdf_ident  = _dedup(pd.DataFrame(identifiers), ["entity_id", "identifier"]) if identifiers else pd.DataFrame()

    return pdf_entity, pdf_coll, pdf_name, pdf_ident



def write_and_preview(spark, database, mode, pdf_entity, pdf_coll, pdf_name, pdf_ident):
    """
    Write tables to Delta and preview results.
    """
    write_delta(spark, pdf_entity, database, "entity", mode)
    write_delta(spark, pdf_coll, database, "contig_collection", mode)
    write_delta(spark, pdf_name, database, "name", mode)
    write_delta(spark, pdf_ident, database, "identifier", mode)

    print("\nDelta tables written:")
    for tbl in ["datasource", "entity", "contig_collection", "name", "identifier"]:
        preview_or_skip(spark, database, tbl)



def main(taxid, api_key, database, mode, debug, allow_genbank_date=False, unique_per_taxon=False):
    spark = build_spark(database)

    # datasource table (fixed record)
    df_ds = build_cdm_datasource()
    write_delta(spark, df_ds, database, "datasource", mode)

    entities, collections, names, identifiers = [], [], [], []
    seen = set()

    taxids = [t.strip() for t in taxid.split(",") if t.strip()]
    print(f"Using TaxIDs: {taxids}")

    for tx in taxids:
        print(f"Fetching taxon={tx}")
        e, c, n, i = process_taxon(tx, api_key, debug, allow_genbank_date, unique_per_taxon, seen)

        # extend results into global containers
        entities.extend(e)
        collections.extend(c)
        names.extend(n)
        identifiers.extend(i)

    pdf_entity, pdf_coll, pdf_name, pdf_ident = finalize_tables(entities, collections, names, identifiers)
    write_and_preview(spark, database, mode, pdf_entity, pdf_coll, pdf_name, pdf_ident)


if __name__ == "__main__":
    cli()





