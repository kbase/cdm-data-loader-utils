"""
Usage: 
python scripts/annotation_parse.py \
  --accession GCF_000869125.1 \
  --output-path output/refseq/GCF_000869125.1 \
  --query

"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import List, Tuple, Optional

import requests
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType
from delta import configure_spark_with_delta_pip

from cdm_data_loader_utils.parsers.kbase_cdm_pyspark import schema as cdm_schemas


# ---------------------------------------------------------------------
# Accession-based annotation fetch
# ---------------------------------------------------------------------
def fetch_annotation_json(accession: str) -> dict:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/annotation_report"
    resp = requests.get(url, headers={"Accept": "application/json"}, timeout=60)
    resp.raise_for_status()
    return resp.json()


# ---------------------------------------------------------------------
# SPARK SESSION
# ---------------------------------------------------------------------
def build_spark_session(app_name: str = "RefSeqAnnotationToCDM") -> SparkSession:
    builder = (
        SparkSession.builder.appName(app_name)
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config(
            "spark.sql.catalog.spark_catalog",
            "org.apache.spark.sql.delta.catalog.DeltaCatalog",
        )
    )
    return configure_spark_with_delta_pip(builder).getOrCreate()


# ---------------------------------------------------------------------
# CDM TABLE SCHEMAS
# ---------------------------------------------------------------------
# Using centralized schemas
IDENTIFIER_SCHEMA = cdm_schemas["Identifier"]
NAME_SCHEMA = cdm_schemas["Name"]
FEATURE_SCHEMA = cdm_schemas["Feature"]
CONTIG_COLLECTION_X_FEATURE_SCHEMA = cdm_schemas["ContigCollection_x_Feature"]
CONTIG_COLLECTION_X_PROTEIN_SCHEMA = cdm_schemas["ContigCollection_x_Protein"]
FEATURE_X_PROTEIN_SCHEMA = cdm_schemas["Feature_x_Protein"]
CONTIG_SCHEMA = cdm_schemas["Contig"]
CONTIG_X_CONTIG_COLLECTION_SCHEMA = cdm_schemas["Contig_x_ContigCollection"]


# ---------------------------------------------------------------------
# CDM PREFIX NORMALIZATION
# ---------------------------------------------------------------------
def apply_prefix(identifier: str) -> str:
    """
    Normalize identifiers to CDM-style prefixed IDs.
    """

    # Normalize GeneID
    if identifier.startswith("GeneID:"):
        return identifier.replace("GeneID:", "ncbigene:")

    # Normalize RefSeq protein accessions
    if identifier.startswith(("YP_", "XP_", "WP_", "NP_", "NC_")):
        return f"refseq:{identifier}"

    # Normalize RefSeq assembly or contig collection IDs
    if identifier.startswith("GCF_"):
        return f"insdc.gcf:{identifier}"

    return identifier


# ---------------------------------------------------------------------
# IDENTIFIERS
# ---------------------------------------------------------------------
def load_identifiers(
    input_json: Path,
) -> List[Tuple[str, str, str, str, Optional[str]]]:
    data = json.loads(input_json.read_text())
    out = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        gene_id = ann.get("gene_id")
        name = ann.get("name")  # used as description

        if not gene_id:
            continue

        entity_id = apply_prefix(f"GeneID:{gene_id}")
        relationship = ann.get("relationship", None)  # may be None

        out.append((entity_id, gene_id, name, "RefSeq", relationship))

    return out


# ---------------------------------------------------------------------
# NAME EXTRACTION
# ---------------------------------------------------------------------
def load_names(input_json: Path) -> list[tuple]:
    data = json.loads(input_json.read_text())
    out = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        gene_id = ann.get("gene_id")

        if not gene_id:
            continue

        entity_id = apply_prefix(f"GeneID:{gene_id}")
        if ann.get("symbol"):
            out.append((entity_id, ann["symbol"], "RefSeq gene symbol", "RefSeq"))
        if ann.get("name"):
            out.append((entity_id, ann["name"], "RefSeq gene name", "RefSeq"))
        if ann.get("locus_tag"):
            out.append((entity_id, ann["locus_tag"], "RefSeq locus tag", "RefSeq"))

    return out


# ---------------------------------------------------------------------
# FEATURE LOCATIONS
# ---------------------------------------------------------------------
def load_feature_records(input_json: Path) -> list[tuple]:
    def to_int(val):  # for safe int conversion
        try:
            return int(val)
        except Exception:
            return None

    data = json.loads(input_json.read_text())
    features = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        gene_id = ann.get("gene_id")

        if not gene_id:
            continue

        feature_id = apply_prefix(f"GeneID:{gene_id}")

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

    return features


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION <-> FEATURE
# ---------------------------------------------------------------------
def load_contig_collection_x_feature(input_json: Path) -> List[Tuple[str, str]]:
    data = json.loads(input_json.read_text())
    links = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        gene_id = ann.get("gene_id")
        regions = ann.get("genomic_regions", [])

        if not gene_id or not regions:
            continue

        acc = regions[0].get("gene_range", {}).get("accession_version")
        if acc:
            links.append((apply_prefix(acc), apply_prefix(f"GeneID:{gene_id}")))

    return links


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION <-> PROTEIN
# ---------------------------------------------------------------------
def load_contig_collection_x_protein(input_json: Path) -> list[tuple]:
    data = json.loads(input_json.read_text())
    links = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        proteins = ann.get("proteins", [])
        annotations = ann.get("annotations", [])

        if not proteins or not annotations:
            continue

        assembly = annotations[0].get("assembly_accession")

        if not assembly:
            continue

        contig_id = apply_prefix(assembly)

        for p in proteins:
            pid = p.get("accession_version")
            if pid:
                protein_id = apply_prefix(pid)
                links.append((contig_id, protein_id))

    return links


# ---------------------------------------------------------------------
# PARSE FEATURE <-> PROTEIN
# ---------------------------------------------------------------------
def load_feature_x_protein(input_json: Path) -> List[Tuple[str, str]]:
    data = json.loads(input_json.read_text())
    links = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        gene_id = ann.get("gene_id")
        proteins = ann.get("proteins", [])

        if not gene_id or not proteins:
            continue

        feature_id = apply_prefix(f"GeneID:{gene_id}")

        for p in proteins:
            pid = p.get("accession_version")
            if pid:
                protein_id = apply_prefix(pid)
                links.append((feature_id, protein_id))

    return links


# ---------------------------------------------------------------------
# PARSE CONTIGS
# ---------------------------------------------------------------------
def load_contigs(
    input_json: Path,
) -> List[Tuple[str, Optional[str], Optional[float], Optional[int]]]:
    data = json.loads(input_json.read_text())

    contigs: dict[str, dict[str, Optional[object]]] = {}

    for report in data.get("reports", []):
        annotation = report.get("annotation", {})

        for region in annotation.get("genomic_regions", []):
            accession = region.get("gene_range", {}).get("accession_version")
            if not accession:
                continue

            contig_id = apply_prefix(accession)

            if contig_id not in contigs:
                contigs[contig_id] = {
                    "hash": None,
                    "gc_content": None,
                    "length": None,
                }

    rows: List[Tuple[str, Optional[str], Optional[float], Optional[int]]] = []

    for contig_id, meta in contigs.items():
        rows.append((
            contig_id,
            meta["hash"],
            meta["gc_content"],
            meta["length"],
        ))

    return rows


# ---------------------------------------------------------------------
# PARSE CONTIG <-> CONTIG_COLLECTION
# ---------------------------------------------------------------------
def load_contig_x_contig_collection(input_json: Path) -> list[tuple]:
    data = json.loads(input_json.read_text())
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
            contig_id = f"refseq:{contig}"
            collection_id = apply_prefix(assembly)
            links.append((contig_id, collection_id))

    return links


# ---------------------------------------------------------------------
# DELTA TABLE
# ---------------------------------------------------------------------
def write_to_delta(
    spark: SparkSession,
    records: list[tuple],
    output_path: str,
    schema: StructType,
):
    if not records:
        return None

    df = spark.createDataFrame(records, schema=schema)
    df.write.format("delta").mode("overwrite").option("overwriteSchema", "true").save(output_path)

    return df


# ---------------------------------------------------------------------
# SQL PREVIEW
# ---------------------------------------------------------------------
def run_sql_query(spark: SparkSession, delta_path: str):
    for name in [
        "cdm_identifiers",
        "cdm_names",
        "cdm_features",
        "cdm_contig_collection_x_feature",
        "cdm_contig_collection_x_protein",
        "cdm_feature_x_protein",
        "cdm_contigs",
        "cdm_contig_x_contig_collection",
    ]:
        print(f"\n[SQL] {name}:")
        path = str(Path(delta_path) / name)
        spark.read.format("delta").load(path).createOrReplaceTempView(name)
        spark.sql(f"SELECT * FROM {name} LIMIT 20").show(truncate=False)


# ---------------------------------------------------------------------
# CLI ENTRY
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="RefSeq Annotation Parser to CDM")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output-path", required=True)
    parser.add_argument("--query", action="store_true")

    args = parser.parse_args()

    base_output = Path(args.output_path)
    base_output.mkdir(parents=True, exist_ok=True)
    data = fetch_annotation_json(args.accession)

    input_path = Path(f"/tmp/{args.accession}.json")
    input_path.write_text(json.dumps(data, indent=2))

    spark = build_spark_session()

    write_to_delta(
        spark,
        load_identifiers(input_path),
        str(base_output / "cdm_identifiers"),
        IDENTIFIER_SCHEMA,
    )
    write_to_delta(spark, load_names(input_path), str(base_output / "cdm_names"), NAME_SCHEMA)
    write_to_delta(
        spark,
        load_feature_records(input_path),
        str(base_output / "cdm_features"),
        FEATURE_SCHEMA,
    )
    write_to_delta(
        spark,
        load_contig_collection_x_feature(input_path),
        str(base_output / "cdm_contig_collection_x_feature"),
        CONTIG_COLLECTION_X_FEATURE_SCHEMA,
    )
    write_to_delta(
        spark,
        load_contig_collection_x_protein(input_path),
        str(base_output / "cdm_contig_collection_x_protein"),
        CONTIG_COLLECTION_X_PROTEIN_SCHEMA,
    )
    write_to_delta(
        spark,
        load_feature_x_protein(input_path),
        str(base_output / "cdm_feature_x_protein"),
        FEATURE_X_PROTEIN_SCHEMA,
    )

    write_to_delta(
        spark,
        load_contigs(input_path),
        str(base_output / "cdm_contigs"),
        CONTIG_SCHEMA,
    )

    write_to_delta(
        spark,
        load_contig_x_contig_collection(input_path),
        str(base_output / "cdm_contig_x_contig_collection"),
        CONTIG_X_CONTIG_COLLECTION_SCHEMA,
    )

    if args.query:
        run_sql_query(spark, str(base_output))

    spark.stop()


if __name__ == "__main__":
    main()
