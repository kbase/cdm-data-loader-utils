"""

RefSeq annotation parser for transforming NCBI Datasets API JSON into CDM-formatted Delta Lake tables.

Usage:
uv run python src/cdm_data_loader_utils/parsers/refseq/api/annotation_report.py \
  --accession GCF_000869125.1 \
  --namespace refseq_api \
  --query

"""

import argparse
import json
from pathlib import Path

import requests
import logging

from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StructField, StringType

from cdm_data_loader_utils.model.kbase_cdm_schema import CDM_SCHEMA


logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------
# Accession-based annotation fetch
# ---------------------------------------------------------------------
def fetch_annotation_json(accession: str) -> dict:
    """Fetch annotation JSON from NCBI Datasets API."""
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/annotation_report"
    try:
        resp = requests.get(url, headers={"Accept": "application/json"}, timeout=60)
        resp.raise_for_status()
        return resp.json()
    except requests.exceptions.HTTPError as http_err:
        logger.error(f"HTTP error for accession {accession}: {http_err}")
        raise
    except requests.exceptions.RequestException as req_err:
        logger.error(f"Request failed for accession {accession}: {req_err}")
        raise


# ---------------------------------------------------------------------
# Spark initialization with Delta support
# ---------------------------------------------------------------------
def build_spark_session(app_name: str = "RefSeqAnnotationToCDM") -> SparkSession:
    """Build and return a SparkSession with Delta Lake and Hive support."""
    try:
        builder = (
            SparkSession.builder.appName(app_name)
            .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
            .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
            .enableHiveSupport()
        )
        return configure_spark_with_delta_pip(builder).getOrCreate()
    except Exception as e:
        logger.error(f"Failed to initialize SparkSession: {e}")
        raise


def init_spark_and_db(app_name: str, database: str) -> SparkSession:
    """Initialize Spark session and set active database."""
    spark = build_spark_session(app_name)
    try:
        spark.sql(f"CREATE DATABASE IF NOT EXISTS {database}")
        spark.sql(f"USE {database}")
    except Exception as e:
        logger.error(f"Failed to create/use database `{database}`: {e}")
        raise
    return spark


# ---------------------------------------------------------------------
# CDM PREFIX NORMALIZATION
# ---------------------------------------------------------------------
def apply_prefix(identifier: str | None) -> str | None:
    """Normalize identifier by applying CDM-compliant prefixes."""
    if not identifier:
        return None

    PREFIX_MAP = {
        "GeneID:": "ncbigene:",
        "YP_": "refseq:",
        "XP_": "refseq:",
        "WP_": "refseq:",
        "NP_": "refseq:",
        "NC_": "refseq:",
        "GCF_": "insdc.gcf:",
    }

    for prefix, replacement in PREFIX_MAP.items():
        if identifier.startswith(prefix):
            # Only replace "GeneID:" completely
            return (
                replacement + identifier[len(prefix) :]
                if ":" not in prefix
                else identifier.replace(prefix, replacement)
            )

    return identifier


# ---------------------------------------------------------------------
# Safe integer conversion
# ---------------------------------------------------------------------
def to_int(val: str) -> int | None:
    try:
        return int(val)
    except (ValueError, TypeError):
        return None


# ---------------------------------------------------------------------
# For repeat section markers
# ---------------------------------------------------------------------
def unique_annotations(data: dict):
    seen = set()
    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        gene_id = ann.get("gene_id")
        if gene_id and gene_id not in seen:
            seen.add(gene_id)
            yield gene_id, ann


# ---------------------------------------------------------------------
# IDENTIFIERS
# ---------------------------------------------------------------------
def load_identifiers(data: dict) -> list[tuple[str, str, str, str, str | None]]:
    """Extract Identifier table records."""
    out = []

    for gene_id, ann in unique_annotations(data):
        entity_id = f"ncbigene:{gene_id}"
        out.append((entity_id, gene_id, ann.get("name"), "RefSeq", ann.get("relationship")))
    return list({tuple(row) for row in out})  # deduplicate


# ---------------------------------------------------------------------
# NAME EXTRACTION
# ---------------------------------------------------------------------
def load_names(data: dict) -> list[tuple[str, str, str, str]]:
    """Extract Name table records."""
    out = []

    for gene_id, ann in unique_annotations(data):
        entity_id = f"ncbigene:{gene_id}"
        for label, desc in [
            ("symbol", "RefSeq gene symbol"),
            ("name", "RefSeq gene name"),
            ("locus_tag", "RefSeq locus tag"),
        ]:
            val = ann.get(label)
            if val:
                out.append((entity_id, val, desc, "RefSeq"))
    return list({tuple(row) for row in out})


# ---------------------------------------------------------------------
# FEATURE LOCATIONS
# ---------------------------------------------------------------------
def load_feature_records(data: dict) -> list[tuple]:
    """Extract Feature table records."""
    features = []

    for gene_id, ann in unique_annotations(data):
        feature_id = f"ncbigene:{gene_id}"
        for region in ann.get("genomic_regions", []):
            for r in region.get("gene_range", {}).get("range", []):
                strand = {
                    "plus": "positive",
                    "minus": "negative",
                    "unstranded": "unstranded",
                }.get(r.get("orientation"), "unknown")
                features.append((
                    feature_id,
                    None,
                    None,
                    None,
                    to_int(r.get("end")),
                    None,
                    to_int(r.get("begin")),
                    strand,
                    "RefSeq",
                    None,
                    "gene",
                ))
    return list({tuple(row) for row in features})


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION <-> FEATURE
# ---------------------------------------------------------------------
def load_contig_collection_x_feature(data: dict) -> list[tuple[str, str]]:
    """Parse ContigCollection Feature links."""
    links = []

    for gene_id, ann in unique_annotations(data):
        regions = ann.get("genomic_regions", [])

        if not regions:
            continue

        acc = regions[0].get("gene_range", {}).get("accession_version")
        if acc:
            links.append((apply_prefix(acc), f"ncbigene:{gene_id}"))

    return list(set(links))


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION <-> PROTEIN %%%
# ---------------------------------------------------------------------
def load_contig_collection_x_protein(data: dict) -> list[tuple[str, str]]:
    links = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        assembly = ann.get("annotations", [{}])[0].get("assembly_accession")
        if not assembly:
            continue

        contig_id = apply_prefix(assembly)
        for p in ann.get("proteins", []):
            pid = p.get("accession_version")
            if pid:
                links.append((contig_id, apply_prefix(pid)))

    return list(set(links))


# ---------------------------------------------------------------------
# PARSE FEATURE <-> PROTEIN
# ---------------------------------------------------------------------
def load_feature_x_protein(data: dict) -> list[tuple[str, str]]:
    links = []

    for gene_id, ann in unique_annotations(data):
        feature_id = f"ncbigene:{gene_id}"

        for p in ann.get("proteins", []):
            pid = p.get("accession_version")
            if pid:
                protein_id = apply_prefix(pid)
                links.append((feature_id, protein_id))

    return list(set(links))


# ---------------------------------------------------------------------
# PARSE CONTIGS
# ---------------------------------------------------------------------
def load_contigs(data: dict) -> list[tuple[str, str | None, float | None, int | None]]:
    contigs = {}

    for report in data.get("reports", []):
        for region in report.get("annotation", {}).get("genomic_regions", []):
            acc = region.get("gene_range", {}).get("accession_version")
            if acc:
                contig_id = apply_prefix(acc)
                # Only track first occurrence of each contig
                contigs.setdefault(contig_id, {"hash": None, "gc_content": None, "length": None})

    return [(cid, meta["hash"], meta["gc_content"], meta["length"]) for cid, meta in contigs.items()]


# ---------------------------------------------------------------------
# PARSE CONTIG <-> CONTIG_COLLECTION
# ---------------------------------------------------------------------
def load_contig_x_contig_collection(data: dict) -> list[tuple[str, str]]:
    links = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        regions = ann.get("genomic_regions", [])
        annotations = ann.get("annotations", [])

        if not regions or not annotations:
            continue

        contig = regions[0].get("gene_range", {}).get("accession_version")
        assembly = annotations[0].get("assembly_accession")

        if contig and assembly:
            links.append((
                f"refseq:{contig}",
                apply_prefix(assembly),
            ))

    return list(set(links))


def load_contig_x_feature(data: dict) -> list[tuple[str, str]]:
    """Extract (contig_id, feature_id) pairs."""
    links = []

    for gene_id, ann in unique_annotations(data):
        feature_id = f"ncbigene:{gene_id}"

        for region in ann.get("genomic_regions", []):
            acc = region.get("gene_range", {}).get("accession_version")
            if acc:
                contig_id = apply_prefix(acc)
                links.append((contig_id, feature_id))

    return list(set(links))


def load_contig_x_protein(data: dict) -> list[tuple[str, str]]:
    links = []

    for _, ann in unique_annotations(data):
        contig_id = None

        for region in ann.get("genomic_regions", []):
            acc = region.get("gene_range", {}).get("accession_version")
            if acc:
                contig_id = apply_prefix(acc)
                break  # only take first

        if contig_id:
            for p in ann.get("proteins", []):
                pid = p.get("accession_version")
                if pid:
                    links.append((contig_id, apply_prefix(pid)))

    return list({tuple(row) for row in links})


### contig collection has 34 rows
CONTIG_COLLECTION_MIN_SCHEMA = StructType([
    StructField("contig_collection_id", StringType(), nullable=False),
    StructField("hash", StringType(), nullable=True),
])


def load_contig_collections(data: dict) -> list[tuple]:
    records = []
    seen_ids = set()

    for report in data.get("reports", []):
        annotations = report.get("annotation", {}).get("annotations", [])
        for ann in annotations:
            raw_accession = ann.get("assembly_accession")
            prefixed_id = apply_prefix(raw_accession)
            if prefixed_id and prefixed_id not in seen_ids:
                seen_ids.add(prefixed_id)

                records.append((prefixed_id, None))

    return records


def load_protein(data: dict) -> list[tuple[str, str | None, str | None, str | None, int | None, str | None]]:
    out = []

    for _, ann in unique_annotations(data):
        for p in ann.get("proteins", []):
            pid = apply_prefix(p.get("accession_version"))
            name = p.get("name")
            length = p.get("length")
            out.append((pid, None, name, None, length, None))

    return list({tuple(row) for row in out})


# ---------------------------------------------------------------------
# DELTA TABLE
# ---------------------------------------------------------------------
def write_to_table(
    spark: SparkSession,
    records: list[tuple],
    table_name: str,
    database: str = "default",
) -> None:
    if not records:
        print(f"[DEBUG] {table_name}: no records to write.")
        return

    # Determine schema
    schema = CONTIG_COLLECTION_MIN_SCHEMA if table_name == "ContigCollection" else CDM_SCHEMA.get(table_name)
    if schema is None:
        raise ValueError(f"[ERROR] Unknown schema: {table_name}")

    df = spark.createDataFrame(records, schema)
    df.printSchema()
    df.show(truncate=False)

    # Write to Delta table
    df.write.format("delta").mode("overwrite").option("overwriteSchema", "true").saveAsTable(f"{database}.{table_name}")

    # Register as temp view for downstream use
    df.createOrReplaceTempView(table_name)


# ---------------------------------------------------------------------
# SQL PREVIEW
# ---------------------------------------------------------------------

CDM_TABLES = [
    "Identifier",
    "Name",
    "Feature",
    "ContigCollection_x_Feature",
    "ContigCollection_x_Protein",
    "Feature_x_Protein",
    "Contig",
    "Contig_x_ContigCollection",
    "ContigCollection",
    "Protein",
    "Contig_x_Feature",
    "Contig_x_Protein",
]


def run_sql_query(spark: SparkSession, database: str = "default", limit: int = 20) -> None:
    """Preview the contents of each CDM table in the given database."""
    spark.sql(f"USE {database}")
    print(f"\n[INFO] Running SQL preview in database: {database}\n")

    for table in CDM_TABLES:
        full_table_name = f"{database}.{table}"
        print(f"\n[SQL Preview] {full_table_name}")
        try:
            df = spark.sql(f"SELECT * FROM {full_table_name} LIMIT {limit}")
            if df.isEmpty():
                print(f"[INFO] Table {full_table_name} is empty.")
            else:
                df.show(truncate=False)
        except Exception as e:
            print(f"[WARNING] Could not query table {full_table_name}: {e}")


def parse_annotation_data(spark: SparkSession, datasets: list[dict], namespace: str) -> None:
    """Parse annotation data into CDM tables and write to Delta Lake."""
    # Mapping of table names to corresponding loader functions
    loader_map = {
        "Identifier": load_identifiers,
        "Name": load_names,
        "Feature": load_feature_records,
        "ContigCollection_x_Feature": load_contig_collection_x_feature,
        "ContigCollection_x_Protein": load_contig_collection_x_protein,
        "Feature_x_Protein": load_feature_x_protein,
        "Contig": load_contigs,
        "Contig_x_ContigCollection": load_contig_x_contig_collection,
        "Contig_x_Feature": load_contig_x_feature,
        "Contig_x_Protein": load_contig_x_protein,
        "ContigCollection": load_contig_collections,
        "Protein": load_protein,
    }

    for data in datasets:
        for table_name, loader_fn in loader_map.items():
            records = loader_fn(data)
            write_to_table(spark, records, table_name, namespace)


# ---------------------------------------------------------------------
# CLI ENTRY
# ---------------------------------------------------------------------
def load_input_data(args) -> list[dict]:
    datasets = []
    if args.accession:
        datasets.append(fetch_annotation_json(args.accession))

    if args.input_file:
        with open(args.input_file) as f:
            datasets.append(json.load(f))

    if args.input_dir:
        for path in Path(args.input_dir).rglob("*.json"):
            with open(path) as f:
                datasets.append(json.load(f))

    return datasets


def main():
    parser = argparse.ArgumentParser(description="RefSeq Annotation Parser to CDM")

    # ------------------------- Input -------------------------
    parser.add_argument("--accession", type=str, help="RefSeq genome accession (e.g. GCF_000869125.1)")
    parser.add_argument("--input_file", type=str, help="Path to a RefSeq annotation JSON file.")
    parser.add_argument("--input_dir", type=str, help="Directory containing RefSeq annotation JSON files.")

    # ------------------------- Output -------------------------
    parser.add_argument("--namespace", default="refseq_api", help="Database to write Delta tables.")
    parser.add_argument("--tenant", default=None, help="Tenant SQL warehouse to use.")
    parser.add_argument("--query", action="store_true", help="Preview SQL output after writing.")

    args = parser.parse_args()

    # ------------------------- Validation -------------------------
    if not args.accession and not args.input_file and not args.input_dir:
        raise ValueError("Provide --accession, --input_file, or --input_dir.")

    # ------------------------- Load data -------------------------
    datasets = load_input_data(args)
    if not datasets:
        raise RuntimeError("No valid annotation datasets were loaded.")

    # ------------------------- Spark  -------------------------
    spark = init_spark_and_db("RefSeq Annotation Parser", args.namespace)
    if args.tenant:
        spark.sql(f"USE CATALOG {args.tenant}")

    # ------------------------- Parse  -------------------------
    parse_annotation_data(spark, datasets, args.namespace)

    # ------------------------- SQL Preview -------------------------
    if args.query:
        run_sql_query(spark, args.namespace)

    spark.stop()


if __name__ == "__main__":
    main()
