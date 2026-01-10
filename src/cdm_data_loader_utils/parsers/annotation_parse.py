"""
Usage:
python3 scripts/annotation_parse.py \
  --input-json tests/data/refseq/GCF_000869125.1.annotation_report.json \
  --output-path output/refseq/GCF_000869125.1/cdm_identifiers \
  --query --join
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StructField, StringType
from delta import configure_spark_with_delta_pip


# ---------------------------------------------------------------------
# CDM TABLE SCHEMAS
# ---------------------------------------------------------------------

# ---------- Identifiers / Feature core ----------
CDM_SCHEMA = StructType([
    StructField("entity_id", StringType(), False),
    StructField("identifier", StringType(), False),
    StructField("name", StringType(), True),
    StructField("name_description", StringType(), True),
    StructField("symbol_description", StringType(), True),
    StructField("source", StringType(), True),
    StructField("tax_id", StringType(), True),
])

# ---------- ContigCollectionXFeature ----------
CONTIG_COLLECTION_X_FEATURE_SCHEMA = StructType([
    StructField("contig_collection_id", StringType(), False),
    StructField("feature_id", StringType(), False),
])

# ---------- genomic_regions ----------
FEATURE_LOCATION_SCHEMA = StructType([
    StructField("feature_id", StringType(), False),
    StructField("start", StringType(), True),
    StructField("end", StringType(), True),
    StructField("strand", StringType(), True),
    StructField("reference_accession", StringType(), True),
    StructField("reference_link", StringType(), True),
])

PROTEIN_IDENTIFIER_SCHEMA = StructType([
    StructField("entity_id", StringType(), False),  # protein_id
    StructField("identifier", StringType(), False),  # protein accession
    StructField("name", StringType(), True),
    StructField("name_description", StringType(), True),
    StructField("length", StringType(), True),
    StructField("isoform_name", StringType(), True),
    StructField("ensembl_protein", StringType(), True),
    StructField("mature_peptides", StringType(), True),
    StructField("source", StringType(), True),
])

# ContigCollectionXProtein
CONTIG_COLLECTION_X_PROTEIN_SCHEMA = StructType([
    StructField("contig_collection_id", StringType(), False),
    StructField("protein_id", StringType(), False),
])

# FeatureXProtein
FEATURE_X_PROTEIN_SCHEMA = StructType([
    StructField("feature_id", StringType(), False),
    StructField("protein_id", StringType(), False),
])


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
# FEATURE / IDENTIFIERS
# ---------------------------------------------------------------------
def parse_gene_identifiers(input_json: Path, prefix: str) -> list[tuple]:
    data = json.loads(input_json.read_text())
    out = []

    for report in data.get("reports", []):
        ann = report.get("annotation", {})
        gene_id = ann.get("gene_id")
        name = ann.get("name")
        tax_id = ann.get("tax_id")

        if not gene_id:
            continue

        out.append((
            f"{gene_id}",  # entity_id
            f"{prefix}:{gene_id}",  # identifier
            name,
            "RefSeq gene name" if name else None,
            "RefSeq gene symbol",
            "RefSeq",
            str(tax_id) if tax_id else None,
        ))

    return out


# ---------------------------------------------------------------------
# CONTIG_COLLECTION <-> FEATURE
# ---------------------------------------------------------------------
def load_contig_collection_links(
    input_json: Path, identifier_prefix: str
) -> list[tuple]:
    """

    Parse RefSeq JSON and generate links between:
    - CDM Contig entities (using accession_version from gene_range)
    - CDM Feature entities (based on gene_id)

    Returns:
        A list of (contig_id, feature_id) tuples

    """

    data = json.loads(input_json.read_text())
    links = []

    for report in data.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        if not gene_id:
            continue  # Skip records with no gene ID

        genomic_regions = annotation.get("genomic_regions", [])
        if not genomic_regions:
            continue  # Skip if no genomic region info

        gene_range = genomic_regions[0].get("gene_range", {})
        accession = gene_range.get("accession_version")
        if not accession:
            continue  # Skip if no contig reference

        # Create link tuple between Contig and Feature
        links.append((
            f"Contig:{accession}",  # CDM Contig ID
            f"{identifier_prefix}:{gene_id}",  # CDM Feature ID
        ))

    return links


# ---------------------------------------------------------------------
# genomic_regions
# ---------------------------------------------------------------------
def load_feature_locations(
    input_json: Path,
    identifier_prefix: str,
) -> list[tuple]:
    """

    Extracts genomic location information for features from a RefSeq annotation JSON.

    Each record corresponds to a genomic region (start, end, strand) for a gene feature,
    along with a reference accession and an optional NCBI link.

    Returns:
        List of tuples in the format: feature_id, start, end, strand, contig_accession, external_reference_url

    """

    data = json.loads(input_json.read_text())
    locations = []

    for report in data.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        if not gene_id:
            continue  # Skip if missing gene ID

        feature_id = f"{identifier_prefix}:{gene_id}"

        for region in annotation.get("genomic_regions", []):
            gene_range = region.get("gene_range", {})
            accession = gene_range.get("accession_version")
            ranges = gene_range.get("range", [])

            for rg in ranges:
                start = rg.get("begin")
                end = rg.get("end")
                orientation = rg.get("orientation")

                # Map orientation to CDM-compliant strand types
                strand = {"plus": "positive", "minus": "negative"}.get(
                    orientation, "unknown"
                )

                locations.append((
                    feature_id,
                    str(start) if start else None,
                    str(end) if end else None,
                    strand,
                    accession,
                    f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}"
                    if accession
                    else None,
                ))

    return locations


# ---------------------------------------------------------------------
# PARSE PROTEIN IDENTIFIERS
# ---------------------------------------------------------------------
def load_protein_identifiers(input_json: Path) -> list[tuple]:
    """
    Extract protein identifier metadata from RefSeq JSON.

    Each record includes the protein's accession, name, length, isoform name,
    Ensembl cross-reference, mature peptide info, and source metadata.

    """

    data = json.loads(input_json.read_text())
    records = []

    for report in data.get("reports", []):
        annotation = report.get("annotation", {})
        proteins = annotation.get("proteins", [])

        for protein in proteins:
            protein_id = protein.get("accession_version")
            if not protein_id:
                continue  # Skip entries without valid accession

            name = protein.get("name")
            length = str(protein.get("length")) if protein.get("length") else None
            isoform_name = protein.get("isoform_name")
            ensembl_prot = protein.get("ensembl_protein", {}).get("protein_id")

            # Aggregate all mature peptide names if available
            peptides = protein.get("mature_peptides", [])
            peptide_str = (
                ", ".join([
                    p.get("product", "")
                    for p in peptides
                    if isinstance(p, dict) and p.get("product")
                ])
                or None
            )

            records.append((
                protein_id,  # entity_id
                protein_id,  # identifier
                name,
                "RefSeq protein name" if name else None,
                length,
                isoform_name,
                ensembl_prot,
                peptide_str,
                "RefSeq",
            ))

    return records


# ---------------------------------------------------------------------
# PARSE CONTIG_COLLECTION <-> PROTEIN
# ---------------------------------------------------------------------
def load_contig_collection_x_protein(input_json: Path) -> list[tuple]:
    """
    Extract Contig ↔ Protein associations from RefSeq JSON.

    Each protein is expected to link to a contig collection (assembly accession),
    which is inferred from the first item in the `annotations` list.

    """
    data = json.loads(input_json.read_text())
    links = []

    for report in data.get("reports", []):
        annotation = report.get("annotation", {})

        proteins = annotation.get("proteins", [])
        annotations = annotation.get("annotations", [])

        if not proteins or not annotations:
            continue

        assembly_accession = annotations[0].get("assembly_accession")
        if not assembly_accession:
            continue

        contig_id = f"Contig:{assembly_accession}"

        for protein in proteins:
            protein_id = protein.get("accession_version")
            if protein_id:
                links.append((contig_id, protein_id))

    return links


# ---------------------------------------------------------------------
# Link between Feature and Protein
# ---------------------------------------------------------------------
def load_feature_x_protein(
    input_json: Path,
    identifier_prefix: str,
) -> list[tuple]:
    """
    Link gene features to their associated proteins using accession IDs.

    Each gene feature (from gene_id) may have one/more protein products.
    The function generates tuples to establishthese associations.
    """

    data = json.loads(input_json.read_text())
    links = []

    for report in data.get("reports", []):
        annotation = report.get("annotation", {})
        gene_id = annotation.get("gene_id")
        proteins = annotation.get("proteins", [])

        if not gene_id or not proteins:
            continue

        feature_id = f"{identifier_prefix}:{gene_id}"

        for protein in proteins:
            protein_id = protein.get("accession_version")
            if protein_id:
                links.append((feature_id, protein_id))

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
        print(f"No valid records to write for {output_path}")
        return None

    df = spark.createDataFrame(records, schema=schema)
    df.write.format("delta").mode("overwrite").option("overwriteSchema", "true").save(
        output_path
    )

    print(f"{df.count()} records written to {output_path}")
    return df


# ---------------------------------------------------------------------
# SQL PREVIEW
# ---------------------------------------------------------------------
def run_sql_query(
    spark: SparkSession,
    delta_path: str,
    run_join: bool = False,
):
    spark.read.format("delta").load(delta_path).createOrReplaceTempView(
        "cdm_identifiers"
    )
    print("\n[SQL] cdm_identifiers:")
    spark.sql("SELECT * FROM cdm_identifiers LIMIT 20").show(truncate=False)

    for name in [
        "contig_collection_x_feature",
        "feature_locations",
        "protein_identifiers",
        "contig_collection_x_protein",
        "feature_x_protein",
    ]:
        path = str(Path(delta_path).parent / name)
        spark.read.format("delta").load(path).createOrReplaceTempView(name)
        print(f"\n[SQL] {name}:")
        spark.sql(f"SELECT * FROM {name} LIMIT 20").show(truncate=False)


# ---------------------------------------------------------------------
# CLI ENTRY
# ---------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Convert RefSeq annotation_report.json to CDM Delta tables"
    )
    parser.add_argument("--input-json", required=True)
    parser.add_argument("--output-path", required=True)
    parser.add_argument("--identifier-prefix", default="GeneID")
    parser.add_argument("--query", action="store_true")
    parser.add_argument("--join", action="store_true")

    args = parser.parse_args()
    spark = build_spark_session()

    # ---------- feature identifiers ----------
    records = parse_gene_identifiers(Path(args.input_json), args.identifier_prefix)
    write_to_delta(spark, records, args.output_path, CDM_SCHEMA)

    # ---------- contig feature ----------
    links = load_contig_collection_links(Path(args.input_json), args.identifier_prefix)
    write_to_delta(
        spark,
        links,
        str(Path(args.output_path).parent / "contig_collection_x_feature"),
        CONTIG_COLLECTION_X_FEATURE_SCHEMA,
    )

    # ---------- feature locations ----------
    locations = load_feature_locations(Path(args.input_json), args.identifier_prefix)
    write_to_delta(
        spark,
        locations,
        str(Path(args.output_path).parent / "feature_locations"),
        FEATURE_LOCATION_SCHEMA,
    )

    # ---------- protein identifiers ----------
    protein_ids = load_protein_identifiers(Path(args.input_json))
    write_to_delta(
        spark,
        protein_ids,
        str(Path(args.output_path).parent / "protein_identifiers"),
        PROTEIN_IDENTIFIER_SCHEMA,
    )

    # ---------- contig to protein ----------
    ccxp = load_contig_collection_x_protein(Path(args.input_json))
    write_to_delta(
        spark,
        ccxp,
        str(Path(args.output_path).parent / "contig_collection_x_protein"),
        CONTIG_COLLECTION_X_PROTEIN_SCHEMA,
    )

    # ---------- feature to protein ----------
    fxp = load_feature_x_protein(Path(args.input_json), args.identifier_prefix)
    write_to_delta(
        spark,
        fxp,
        str(Path(args.output_path).parent / "feature_x_protein"),
        FEATURE_X_PROTEIN_SCHEMA,
    )

    if args.query:
        run_sql_query(spark, args.output_path)

    spark.stop()


if __name__ == "__main__":
    main()
