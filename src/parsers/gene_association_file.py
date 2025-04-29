"""
GO Gene Association File (GAF) parser.

This script processes Gene Ontology (GO) annotation data from GAF-style CSV files and normalizes it into a structured format compatible with ontology-driven data models.

Main Functionality:
--------------------
- Loads raw GO annotations from an input CSV file.
- Drops unused or redundant columns to simplify the dataset.
- Constructs a unified 'subject' identifier by combining 'DB' and 'DB_Object_ID'.
- Renames selected columns to match a target schema (e.g., predicate, object, supporting_objects).
- Converts annotation dates into standard YYYY-MM-DD format.
- Adds metadata fields such as 'aggregator' and 'protocol_id'.
- Determines predicate negation based on a whitelist of allowed GO relations.
- Maps evidence codes to ECO terms by merging with an external mapping file.
- Applies fallback evidence types when no direct evidence mapping is found.
- Outputs a clean, standardized CSV file for downstream processing.

Resources:
----------
- Gene Association File (GAF) format specification: https://geneontology.org/docs/go-annotation-file-gaf-format-2.2/
- Downloadable GAF files: http://current.geneontology.org/annotations/index.html
- ECO (Evidence and Conclusion Ontology) documentation: https://www.evidenceontology.org/annotation_resources/#eco_go_mappings
- GAF-to-ECO mapping file: http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt

Usage:
------
Run via CLI:
    python -m parsers.gene_association_file --input annotations_data100.csv --output normalized_annotation.csv

"""

import pandas as pd
import click
import os 
import sys

# --- Constants ---

SUBJECT = "subject"
PREDICATE = "predicate"
OBJECT = "object"
PUBLICATIONS = "publications"
EVIDENCE_CODE = "Evidence_Code"
SUPPORTING_OBJECTS = "supporting_objects"
ANNOTATION_DATE = "annotation_date"
PRIMARY_KNOWLEDGE_SOURCE = "primary_knowledge_source"
AGGREGATOR = "aggregator"
PROTOCOL_ID = "protocol_id"
NEGATED = "negated"
EVIDENCE_TYPE = "evidence_type"

ASSOCIATION_COL_HEADERS = [
    SUBJECT,
    PREDICATE,
    OBJECT,
    PUBLICATIONS,
    EVIDENCE_CODE,
    SUPPORTING_OBJECTS,
    ANNOTATION_DATE,
    PRIMARY_KNOWLEDGE_SOURCE,
    AGGREGATOR,
    PROTOCOL_ID,
    NEGATED,
    EVIDENCE_TYPE,
]

ASSOCIATION_COL_TYPES = {
    SUBJECT: str,
    PREDICATE: str,
    OBJECT: str,
    PUBLICATIONS: str,
    SUPPORTING_OBJECTS: str,
    ANNOTATION_DATE: str,
    PRIMARY_KNOWLEDGE_SOURCE: str,
    AGGREGATOR: str,
    PROTOCOL_ID: str,
    EVIDENCE_TYPE: str,
    NEGATED: bool,
}

ALLOWED_PREDICATES = [
    "enables", "contributes_to", "acts_upstream_of_or_within", "involved_in",
    "acts_upstream_of", "acts_upstream_of_positive_effect", "acts_upstream_of_negative_effect",
    "acts_upstream_of_or_within_negative_effect", "acts_upstream_of_or_within_positive_effect",
    "located_in", "part_of", "is_active_in", "colocalizes_with"
]

REQUIRED_INPUT_COLUMNS = {
    "DB", "DB_Object_ID", "Qualifier", "GO_ID",
    "DB_Reference", "Evidence_Code", "With_From",
    "Date", "Assigned_By"
}

DROP_COLUMNS = [
    "DB_Object_Symbol", "Aspect", "DB_Object_Name", "Synonym",
    "DB_Object_Type", "Taxon", "Annotation_Extension", "Gene_Product_Form_ID"
]


# --- Load ECO Mapping ---
ECO_MAPPING_URL = "http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt"
ECO_MAPPING_COLUMNS = ["Evidence_Code", "DB_Reference", "evidence_type"]

try:
    eco_mapping_df = pd.read_csv(ECO_MAPPING_URL, sep="\t", comment="#", header=None, names=ECO_MAPPING_COLUMNS)
except Exception as e:
    raise ConnectionError(f"Failed to load ECO mapping from {ECO_MAPPING_URL}: {e}")


# --- Helper Functions ---
    
"""
Validate that the DataFrame matches the required annotation schema.

This function checks:
- All required columns are present.
- Each column has the correct data type.
If any mismatch is found, it raises a ValueError.

"""

def validate_annotation_schema(df: pd.DataFrame):
    """Validates output DataFrame schema"""
    missing_cols = set(ASSOCIATION_COL_TYPES.keys()) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    for column, expected_type in ASSOCIATION_COL_TYPES.items():
        if not df[column].dropna().map(lambda x: isinstance(x, expected_type)).all():
            raise ValueError(f"Invalid data type detected in column: {column}. Expected {expected_type}.")


def merge_evidence_mapping(df: pd.DataFrame, evidence_df: pd.DataFrame) -> pd.DataFrame:

    """
    Merge GO annotation DataFrame with ECO evidence mapping.
    This function matches GO annotations with their corresponding evidence types based on both Evidence_Code and publications fields.

    """

    # Ensure DataFrames are copies to avoid modifying original data
    df = df.copy()
    evidence_df = evidence_df.copy()

    # Normalize evidence codes
    if df["publications"].apply(lambda x: isinstance(x, list)).any():
        df = df.explode("publications")

    # If publications are NaN, drop those rows
    df = df[df["publications"].notna()]

    # Normalize case and strip whitespace
    df["publications"] = df["publications"].str.strip().str.upper()
    evidence_df["DB_Reference"] = evidence_df["DB_Reference"].str.strip().str.upper()
    evidence_df["Evidence_Code"] = evidence_df["Evidence_Code"].str.strip().str.upper()
    
    # Merge with ECO mapping
    merged = df.merge(
        evidence_df,
        how="left",
        left_on=["Evidence_Code", "publications"],
        right_on=["Evidence_Code", "DB_Reference"]
    )

    # Handle unmatched evidence types
    unmatched_mask = merged["evidence_type"].isna()
    fallback_df = evidence_df[evidence_df["DB_Reference"] == "DEFAULT"]

    # Handle unmatched evidence types with fallback
    if not fallback_df.empty:
        fallback_match = merged.loc[unmatched_mask, ["Evidence_Code"]].merge(
            fallback_df, how="left", on="Evidence_Code")
        # Update unmatched rows with fallback evidence type
        merged.loc[unmatched_mask, "evidence_type"] = fallback_match["evidence_type"].values
    
    merged = merged.drop(columns=["DB_Reference"])
    return merged


# --- Core Processing ---

def process_go_annotations(input_path: str, output_path: str):
    """Processes GO annotation file and outputs normalized CSV."""
    go_annotations_df = pd.read_csv(input_path)

    # Validate input columns
    if not REQUIRED_INPUT_COLUMNS.issubset(go_annotations_df.columns):
        # Check for missing columns
        missing = REQUIRED_INPUT_COLUMNS - set(go_annotations_df.columns)
        raise ValueError(f"Missing required columns in input: {missing}")

    # Drop unnecessary columns
    go_annotations_df = go_annotations_df.drop(columns=[col for col in DROP_COLUMNS if col in go_annotations_df.columns])

    # format the subject column
    go_annotations_df["subject"] = (
        go_annotations_df["DB"].fillna("").astype(str) + ":" + go_annotations_df["DB_Object_ID"].fillna("").astype(str)
        )
    
    # Ensure subject is a string
    go_annotations_df["subject"] = go_annotations_df["subject"].astype(str)

    # Drop DB and DB_Object_ID columns
    go_annotations_df = go_annotations_df.drop(columns=["DB", "DB_Object_ID"])
    cols = ["subject"] + [col for col in go_annotations_df.columns if col != "subject"]

    # Reorder columns as final output
    go_annotations_df = go_annotations_df[cols]

    # Rename columns
    go_annotations_df.rename(columns={
        "Qualifier": "predicate",
        "GO_ID": "object",
        "DB_Reference": "publications",
        "With_From": "supporting_objects",
        "Date": "annotation_date",
        "Assigned_By": "primary_knowledge_source"
    }, inplace=True)
    
    try:
        go_annotations_df["annotation_date"] = pd.to_datetime(go_annotations_df["annotation_date"], format='%Y%m%d').dt.strftime('%Y-%m-%d')
    except Exception as e:
        raise ValueError(f"Error parsing annotation_date: {e}")
    
    # Ensure the date is in string format
    go_annotations_df["annotation_date"] = go_annotations_df["annotation_date"].astype(str)

    # Add metadata columns
    go_annotations_df["aggregator"] = "UniProt"
    go_annotations_df["protocol_id"] = None

    # Process multi-value fields
    go_annotations_df["publications"] = go_annotations_df["publications"].apply(lambda x: str(x).split('|') if pd.notna(x) else None)
    go_annotations_df["supporting_objects"] = go_annotations_df["supporting_objects"].apply(lambda x: str(x).split('|') if pd.notna(x) else None)

    # Detect negation
    go_annotations_df["negated"] = go_annotations_df["predicate"].str.contains(r"^NOT\|", regex=True)
    go_annotations_df["predicate"] = go_annotations_df["predicate"].str.replace(r"^NOT\|", "", regex=True)

    # Validate allowed predicates
    invalid_predicates = ~go_annotations_df["predicate"].isin(ALLOWED_PREDICATES)
    if invalid_predicates.any():
        raise ValueError(f"Invalid predicates found: {go_annotations_df.loc[invalid_predicates, 'predicate'].unique()}")

    # Merge with evidence types
    merged = merge_evidence_mapping(go_annotations_df, eco_mapping_df)

    # --- Fix multi-value columns ---
    for col in ["publications", "supporting_objects"]:
        if col in merged.columns:
            merged[col] = merged[col].apply(lambda x: '|'.join(x) if isinstance(x, list) else x)

    # --- Fix subject column ---
    if "subject" in merged.columns:
        merged["subject"] = merged["subject"].astype(str)

    # Validate schema
    validate_annotation_schema(merged)

    # Save
    merged = merged[ASSOCIATION_COL_HEADERS]  
    merged.to_csv(output_path, index=False)


# --- CLI Interface ---

@click.command()
@click.option('--input', '-i', required=True, help='Path to input GO annotation CSV file.')
@click.option('--output', '-o', required=True, help='Path to output normalized CSV file.')

def main(input, output):
    """CLI entry point to process GO annotations."""

    # Check input file exists
    if not os.path.isfile(input):
        print(f"Error: Input file '{input}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Check output directory exists
    output_dir = os.path.dirname(output) or '.'
    if not os.path.isdir(output_dir):
        print(f"Error: Output directory '{output_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Check ECO mapping file loads
    eco_mapping_url = "http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt"
    try:
        case = pd.read_csv(eco_mapping_url, sep="\t", comment="#", header=None)
    except Exception as e:
        print(f"Error: Failed to load ECO mapping from {eco_mapping_url}.", file=sys.stderr)
        print(str(e), file=sys.stderr)
        sys.exit(1)

    # If all checks passed, then process 
    process_go_annotations(input, output)


if __name__ == '__main__':
    main()

