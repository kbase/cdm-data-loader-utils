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
import os

from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StructField, StringType

from cdm_data_loader_utils.model.kbase_cdm_schema import CDM_SCHEMA
from cdm_data_loader_utils.utils.spark_delta import set_up_workspace
from cdm_data_loader_utils.utils.spark_delta import write_delta


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

DELTA_BASE_DIR = Path("/tmp/cdm_delta")


# ============================================================
# Spark init
# ============================================================
def init_spark(app_name: str, namespace: str, tenant: str | None):
    # ---------- BERDL  ----------
    if os.getenv("BERDL_POD_IP"):
        return set_up_workspace(
            app_name=app_name,
            namespace=namespace,
            tenant_name=tenant,
        )

    # ---------- local ----------
    builder = (
        SparkSession.builder.appName(app_name)
        .master("local[*]")
        .config(
            "spark.sql.extensions",
            "io.delta.sql.DeltaSparkSessionExtension",
        )
        .config(
            "spark.sql.catalog.spark_catalog",
            "org.apache.spark.sql.delta.catalog.DeltaCatalog",
        )
    )

    spark = configure_spark_with_delta_pip(builder).getOrCreate()
    return spark, namespace


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
# CDM PREFIX NORMALIZATION
# ---------------------------------------------------------------------
def apply_prefix(identifier: str) -> str:
    """Standardize ID with known prefixes. Preserve RefSeq accessions like YP_123456 and NC_012345."""
    if identifier.startswith(("YP_", "XP_", "WP_", "NP_", "NC_")):
        return f"refseq:{identifier}"
    if identifier.startswith(("GCF_", "GCA_")):
        return f"insdc.gcf:{identifier}"
    if identifier.isdigit():
        return f"ncbigene:{identifier}"
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
    records: set[tuple] = set()

    for _, ann in unique_annotations(data):
        # 1. Genome-level (GCF_)
        assembly_accession = ann.get("assembly_accession")
        if assembly_accession and assembly_accession.startswith("GCF_"):
            genome_id = apply_prefix(assembly_accession)
            records.add((
                genome_id,
                genome_id,
                "RefSeq genome ID",
                "RefSeq",
                None,
            ))

        # 2. Contig/Assembly (NC_) from gene_range
        for region in ann.get("genomic_regions", []):
            gene_range = region.get("gene_range", {})
            acc = gene_range.get("accession_version")
            if acc and acc.startswith("NC_"):
                contig_id = apply_prefix(acc)
                records.add((
                    contig_id,
                    contig_id,
                    "RefSeq assembly ID",
                    "RefSeq",
                    None,
                ))

        # 3. Gene ID
        gene_id = ann.get("gene_id")
        if gene_id:
            gene_entity = apply_prefix(gene_id)
            records.add((
                gene_entity,
                gene_entity,
                "NCBI gene ID",
                "RefSeq",
                None,
            ))

        # 4. Protein ID (YP)
        for p in ann.get("proteins", []):
            pid = p.get("accession_version")
            if pid:
                protein_id = apply_prefix(pid)
                records.add((
                    protein_id,
                    protein_id,
                    "RefSeq protein ID",
                    "RefSeq",
                    None,
                ))

    return list(records)


# ---------------------------------------------------------------------
# NAME EXTRACTION
# ---------------------------------------------------------------------
def load_names(data: dict) -> list[tuple[str, str, str, str]]:
    """Extract Name table records."""
    records = set()
    for gene_id, ann in unique_annotations(data):
        entity_id = f"ncbigene:{gene_id}"
        for label, desc in (
            ("symbol", "RefSeq gene symbol"),
            ("name", "RefSeq gene name"),
            ("locus_tag", "RefSeq locus tag"),
        ):
            val = ann.get(label)
            if val:
                records.add((entity_id, val, desc, "RefSeq"))
    return list(records)


# ---------------------------------------------------------------------
# FEATURE LOCATIONS
# ---------------------------------------------------------------------
def load_feature_records(data: dict) -> list[tuple]:
    """Extract Feature table records from RefSeq annotation JSON."""
    records = set()

    for gene_id, ann in unique_annotations(data):
        feature_id = f"ncbigene:{gene_id}"
        gene_type = ann.get("gene_type", "unknown")

        for region in ann.get("genomic_regions", []):
            for r in region.get("gene_range", {}).get("range", []):
                strand = {
                    "plus": "positive",
                    "minus": "negative",
                    "unstranded": "unstranded",
                }.get(r.get("orientation"), "unknown")

                records.add((
                    feature_id,  # feature_id
                    None,  # hash
                    None,  # cds_phase
                    None,  # e_value
                    to_int(r.get("end")),  # end
                    None,  # p_value
                    to_int(r.get("begin")),  # start
                    strand,  # strand
                    "ncbigene",  # source_database
                    None,  # protocol_id
                    gene_type,  # type: from JSON "protein-coding"
                ))

    return list(records)


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION <-> FEATURE
# ---------------------------------------------------------------------
def load_contig_collection_x_feature(data: dict) -> list[tuple[str, str]]:
    """Parse ContigCollection <-> Feature links from RefSeq annotations."""
    links = set()

    for gene_id, ann in unique_annotations(data):
        # extract assembly accession from annotations
        annotations = ann.get("annotations", [])
        if not annotations:
            continue

        accession = annotations[0].get("assembly_accession")
        if accession:
            contig_collection_id = f"insdc.gcf:{accession}"
            feature_id = f"ncbigene:{gene_id}"
            links.add((contig_collection_id, feature_id))

    return list(links)


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION <-> PROTEIN
# ---------------------------------------------------------------------
def load_contig_collection_x_protein(data: dict) -> list[tuple[str, str]]:
    """Extract links between ContigCollection and Protein."""
    links = set()

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        assembly = ann.get("annotations", [{}])[0].get("assembly_accession")
        if not assembly:
            continue

        contig_id = apply_prefix(assembly)
        for p in ann.get("proteins", []):
            pid = p.get("accession_version")
            if pid:
                links.add((contig_id, apply_prefix(pid)))

    return list(links)


# ---------------------------------------------------------------------
# PARSE FEATURE <-> PROTEIN
# ---------------------------------------------------------------------
def load_feature_x_protein(data: dict) -> list[tuple[str, str]]:
    """Extract Feature ↔ Protein links."""
    links = set()

    for gene_id, ann in unique_annotations(data):
        feature_id = f"ncbigene:{gene_id}"

        for protein in ann.get("proteins", []):
            pid = protein.get("accession_version")
            if pid:
                links.add((feature_id, apply_prefix(pid)))

    return list(links)


# ---------------------------------------------------------------------
# PARSE CONTIGS
# ---------------------------------------------------------------------
def load_contigs(data: dict) -> list[tuple[str, str | None, float | None, int | None]]:
    contigs = {}

    for report in data.get("reports", []):
        for region in report.get("annotation", {}).get("genomic_regions", []):
            accession = region.get("gene_range", {}).get("accession_version")
            if accession:
                contig_id = apply_prefix(accession)
                # Only track first occurrence of each contig
                contigs.setdefault(contig_id, {"hash": None, "gc_content": None, "length": None})

    return [(contig_id, meta["hash"], meta["gc_content"], meta["length"]) for contig_id, meta in contigs.items()]


# ---------------------------------------------------------------------
# PARSE CONTIG <-> CONTIG_COLLECTION
# ---------------------------------------------------------------------
def load_contig_x_contig_collection(data: dict) -> list[tuple[str, str]]:
    links = set()

    for report in data.get("reports", []):
        ann = report.get("annotation", {})

        regions = ann.get("genomic_regions", [])
        annotations = ann.get("annotations", [])

        if not regions or not annotations:
            continue

        contig = regions[0].get("gene_range", {}).get("accession_version")
        assembly = annotations[0].get("assembly_accession")

        if not contig or not assembly:
            continue

        links.add((
            f"refseq:{contig}",
            apply_prefix(assembly),
        ))

    return list(links)


# ---------------------------------------------------------------------
# PARSE CONTIG <-> FEATURE
# ---------------------------------------------------------------------
def load_contig_x_feature(data: dict) -> list[tuple[str, str]]:
    """Extract (contig_id, feature_id) pairs."""
    links = set()

    for gene_id, ann in unique_annotations(data):
        feature_id = f"ncbigene:{gene_id}"

        for region in ann.get("genomic_regions", []):
            acc = region.get("gene_range", {}).get("accession_version")

            if not acc:
                continue

            links.add((
                apply_prefix(acc),
                feature_id,
            ))

    return list(links)


# ---------------------------------------------------------------------
# PARSE CONTIG <-> PROTEIN
# ---------------------------------------------------------------------
def load_contig_x_protein(data: dict) -> list[tuple[str, str]]:
    links = set()

    for _, ann in unique_annotations(data):
        contig_id = None

        for region in ann.get("genomic_regions", []):
            acc = region.get("gene_range", {}).get("accession_version")
            if acc:
                contig_id = apply_prefix(acc)
                break

        if contig_id:
            for protein in ann.get("proteins", []):
                pid = protein.get("accession_version")
                if pid:
                    protein_id = apply_prefix(pid)
                    links.add((contig_id, protein_id))

    return list(links)


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION
# ---------------------------------------------------------------------
def load_contig_collections(data: dict) -> list[dict]:
    records = {}

    for report in data.get("reports", []):
        for ann in report.get("annotation", {}).get("annotations", []):
            accession = ann.get("assembly_accession")
            if not accession:
                continue

            collection_id = apply_prefix(accession)

            records[collection_id] = {
                "contig_collection_id": collection_id,
                "hash": None,
            }

    return list(records.values())


# ---------------------------------------------------------------------
# PARSE PROTEIN
# ---------------------------------------------------------------------
def load_protein(data: dict) -> list[tuple[str, None, str | None, None, int | None, None]]:
    """
    Extract Protein table records.
    """
    records = set()

    for _, ann in unique_annotations(data):
        for protein in ann.get("proteins", []):
            accession = protein.get("accession_version")
            if not accession:
                continue

            protein_id = apply_prefix(accession)
            name = protein.get("name")
            length = protein.get("length")

            records.add((protein_id, None, name, None, length, None))

    return list(records)


# ---------------------------------------------------------------------
# DELTA TABLE
# ---------------------------------------------------------------------
def write_to_table(
    spark,
    records,
    table_name,
    namespace,
    use_metastore: bool = False,
    mode: str = "overwrite",
):
    if not records:
        logger.info(f"[SKIP] {table_name}: empty")
        return

    schema = CDM_SCHEMA.get(table_name)
    if schema is None:
        raise ValueError(f"Unknown CDM table: {table_name}")

    df = spark.createDataFrame(records, schema)

    if use_metastore:
        full_table = f"{namespace}.{table_name}"

        # Explicitly drop to avoid Hive schema conflicts
        spark.sql(f"DROP TABLE IF EXISTS {full_table}")

        write_delta(
            spark=spark,
            sdf=df,
            delta_ns=namespace,
            table=table_name,
            mode=mode,
        )

        logger.info(
            "Recreated managed Delta table %s.%s via Hive metastore",
            namespace,
            table_name,
        )

    else:
        # Path-based Delta (local)
        table_path = DELTA_BASE_DIR / namespace / table_name
        table_path.parent.mkdir(parents=True, exist_ok=True)

        (df.write.format("delta").mode(mode).option("overwriteSchema", "true").save(str(table_path)))

        logger.info(f"{table_name} -> {table_path}")


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


def run_sql_query(spark, namespace, limit=20):
    print(f"\n[INFO] Preview Delta tables under {namespace}\n")

    for table in CDM_TABLES:
        path = DELTA_BASE_DIR / namespace / table
        print(f"\n[Preview] {table}")

        if not path.exists():
            print("(missing)")
            continue

        df = spark.read.format("delta").load(str(path))
        if df.isEmpty():
            print("(empty)")
        else:
            df.show(limit, truncate=False)


def parse_annotation_data(
    spark: SparkSession,
    datasets: list[dict],
    namespace: str,
    use_metastore: bool = False,
) -> None:
    """
    Parse annotation data into CDM tables and write to Delta Lake.
    """

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
            write_to_table(
                spark,
                records,
                table_name,
                namespace,
                use_metastore=use_metastore,
            )


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
    parser.add_argument("--accession", type=str)
    parser.add_argument("--input_file", type=str)
    parser.add_argument("--input_dir", type=str)

    # ------------------------- Output -------------------------
    parser.add_argument("--namespace", default="refseq_api")
    parser.add_argument("--tenant", default=None)
    parser.add_argument("--query", action="store_true")

    parser.add_argument(
        "--use-metastore",
        action="store_true",
        help="Write Delta tables using Hive metastore (BERDL / prod mode)",
    )

    args = parser.parse_args()

    if not args.accession and not args.input_file and not args.input_dir:
        raise ValueError("Provide --accession, --input_file, or --input_dir.")

    datasets = load_input_data(args)
    if not datasets:
        raise RuntimeError("No valid annotation datasets were loaded.")

    spark, delta_ns = init_spark(
        "RefSeq Annotation Parser",
        args.namespace,
        args.tenant,
    )

    parse_annotation_data(
        spark,
        datasets,
        delta_ns,
        use_metastore=args.use_metastore,
    )

    if args.query:
        run_sql_query(spark, delta_ns)

    spark.stop()


if __name__ == "__main__":
    main()
