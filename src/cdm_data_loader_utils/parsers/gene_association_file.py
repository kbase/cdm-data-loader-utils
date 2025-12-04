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
    PYTHONPATH=src python -m parsers.gene_association_file \
  --input tests/data/annotations_data100.csv \
  --output normalized_annotation.csv
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

# --- Data Types ---
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

# GAF Field Names
DB = "DB"
DB_OBJ_ID = "DB_Object_ID"
DB_OBJ_SYMBOL = "DB_Object_Symbol"
QUALIFIER = "Qualifier"
GO_ID = "GO_ID"
DB_REF = "DB_Reference"
WITH_FROM = "With_From"
ASPECT = "Aspect"
DB_OBJ_NAME = "DB_Object_Name"
SYNONYM = "Synonym"
DB_OBJ_TYPE = "DB_Object_Type"
TAXON = "Taxon"
DATE = "Date"
ASSIGNED_BY = "Assigned_By"
ANNOTATION_EXTENSION = "Annotation_Extension"
GENE_PRODUCT_FORM = "Gene_Product_Form_ID"

GAF_COLUMNS = {
    DB: "DB",
    DB_OBJ_ID: "DB_Object_ID",
    DB_OBJ_SYMBOL: "DB_Object_Symbol",
    QUALIFIER: "Qualifier",
    GO_ID: "GO_ID",
    DB_REF: "DB_Reference",
    EVIDENCE_CODE: "Evidence_Code",
    WITH_FROM: "With_From",
    ASPECT: "Aspect",
    DB_OBJ_NAME: "DB_Object_Name",
    SYNONYM: "Synonym",
    DB_OBJ_TYPE: "DB_Object_Type",
    TAXON: "Taxon",
    DATE: "Date",
    ASSIGNED_BY: "Assigned_By",
    ANNOTATION_EXTENSION: "Annotation_Extension",
    GENE_PRODUCT_FORM: "Gene_Product_Form_ID",
}

REQUIRED_INPUT_COLUMNS = {
    DB,
    DB_OBJ_ID,
    QUALIFIER,
    GO_ID,
    DB_REF,
    EVIDENCE_CODE,
    WITH_FROM,
    DATE,
    ASSIGNED_BY,
}

DROP_COLUMNS = {
    DB_OBJ_SYMBOL,
    ASPECT,
    DB_OBJ_NAME,
    SYNONYM,
    DB_OBJ_TYPE,
    TAXON,
    ANNOTATION_EXTENSION,
    GENE_PRODUCT_FORM,
}

# --- Load ECO Mapping ---
# The file is a headerless plain text TSV file, so we have to specify the column names artificially.
ECO_MAPPING_URL = "http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt"
ECO_MAPPING_COLUMNS = [EVIDENCE_CODE, DB_REF, EVIDENCE_TYPE]


def load_eco_mapping():
    try:
        eco_mapping_df = pd.read_csv(ECO_MAPPING_URL, sep="\t", comment="#", header=None, names=ECO_MAPPING_COLUMNS)
    except Exception as e:
        raise ConnectionError(f"Failed to load ECO mapping from {ECO_MAPPING_URL}: {e}")


# --- Helper Functions ---
def validate_annotation_schema(df: pd.DataFrame):
    """
    Validate that the DataFrame matches the required annotation schema.

    This function checks:
    - All required columns are present.
    - Each column has the correct data type.
    If any mismatch is found, it raises a ValueError.
    """
    # Find missing columns in the DataFrame
    missing_cols = set(ASSOCIATION_COL_TYPES.keys()) - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    for column, expected_type in ASSOCIATION_COL_TYPES.items():
        """
        Check all match the expected types
        Throw an exception if any of the elements in a column are wrong type 
        """
        if not df[column].map(type).eq(expected_type).all():
            raise ValueError(f"Invalid type in column '{column}': expected {expected_type}")


def check_header_data(df, header_names=GAF_COLUMNS):
    """
    Prevents duplicate header lines in a file if CSV spliced multiple times.
    Check if DataFrame contains any rows that exactly match the header names and remove them.
    """
    # Convert header names to a list for comparison
    header_names = list(header_names.values())
    # check the entire line is identical to the GAF_COLUMNS header
    mask = df.apply(lambda row: (row == header_names).all(), axis=1)

    if mask.any():
        return df[~mask]
    return df


def load_go_file(input_path):
    """
    Load and preprocess a Gene Association File (GAF) from the given input path.
    - Reads the GAF file into a pandas DataFrame using predefined GAF column names.
    - Removes rows that may be duplicated headers.
    - Extracts and returns only the relevant columns required for downstream processing,
      including DB and DB_Object_ID used to construct the 'subject' identifier.
    """
    # Read the GAF file with no header and assign column names
    df = pd.read_csv(input_path, header=None, names=GAF_COLUMNS)

    # Remove any rows that are mistakenly duplicated headers
    df = check_header_data(df)

    # Define relevant columns, required input and columns needed to construct 'subject'
    RELEVANT_COLUMNS = REQUIRED_INPUT_COLUMNS | {DB, DB_OBJ_ID}

    # Determine which of the relevant columns are actually present in the file
    actual_cols = list(RELEVANT_COLUMNS & set(df.columns))

    # Select only the relevant columns
    return df[actual_cols].copy()


def normalize_dates(df):
    """
    Converts the 'annotation_date' column from the string format YYYYMMDD to YYYY-MM-DD.
    Set to None if the date is illegal or in the wrong format.
    """
    # make sure the column is a string
    date_strs = df[ANNOTATION_DATE].astype(str)

    # detect invalid date formats (8 digits)
    invalid_mask = ~date_strs.str.match(r"^\d{8}$")
    if invalid_mask.any():
        print(df.loc[invalid_mask, ANNOTATION_DATE].unique())

    # Convert valid parts to standard date format
    df.loc[~invalid_mask, ANNOTATION_DATE] = pd.to_datetime(
        date_strs[~invalid_mask], format="%Y%m%d", errors="coerce"
    ).dt.strftime("%Y-%m-%d")
    df.loc[invalid_mask, ANNOTATION_DATE] = None
    return df


def process_predicates(df):
    """
    Detect negative relationships in 'predicate' (starting with 'NOT|')
    Clean the 'predicate' field into canonical GO relationship names
    """
    predicates = df[PREDICATE].astype(str)
    df[NEGATED] = predicates.str.startswith("NOT|")
    df[PREDICATE] = predicates.str.replace(r"^NOT\|", "", regex=True)
    return df


def transform_go_data(df):
    """
    Transform a raw GO annotation DataFrame into a standardized format for downstream use.
    - Renames specific GAF columns to semantic schema-compliant names.
    - Normalizes annotation dates.
    - Adds fixed metadata columns ('aggregator', 'protocol_id').
    - Splits pipe-delimited strings in 'publications' and 'supporting_objects' into lists.
    - Processes and validates predicates using domain rules.
    - Constructs the 'subject' field by combining 'DB' and 'DB_Object_ID'.
    - Ensures all required columns are present (inserting None where missing).
    - Reorders the columns to match the expected output format.
    """

    df.rename(
        columns={
            QUALIFIER: PREDICATE,
            GO_ID: OBJECT,
            DB_REF: PUBLICATIONS,
            WITH_FROM: SUPPORTING_OBJECTS,
            DATE: ANNOTATION_DATE,
            ASSIGNED_BY: PRIMARY_KNOWLEDGE_SOURCE,
        },
        inplace=True,
    )

    df = normalize_dates(df)
    df[AGGREGATOR] = "UniProt"
    df[PROTOCOL_ID] = None

    # Split '|' separated publication/supporting_objects into lists
    df[PUBLICATIONS] = df[PUBLICATIONS].apply(lambda x: str(x).split("|") if pd.notna(x) and isinstance(x, str) else x)
    df[SUPPORTING_OBJECTS] = df[SUPPORTING_OBJECTS].apply(
        lambda x: str(x).split("|") if pd.notna(x) and isinstance(x, str) else x
    )
    df = process_predicates(df)

    ALLOWED_PREDICATES = [
        "enables",
        "contributes_to",
        "acts_upstream_of_or_within",
        "involved_in",
        "acts_upstream_of",
        "acts_upstream_of_positive_effect",
        "acts_upstream_of_negative_effect",
        "acts_upstream_of_or_within_negative_effect",
        "acts_upstream_of_or_within_positive_effect",
        "located_in",
        "part_of",
        "is_active_in",
        "colocalizes_with",
    ]

    # Validate predicates
    invalid_predicates = ~df[PREDICATE].isin(ALLOWED_PREDICATES)
    if invalid_predicates.any():
        raise ValueError(f"Invalid predicates found: {df.loc[invalid_predicates, PREDICATE].unique()}")

    df[SUBJECT] = df[DB].astype(str) + ":" + df[DB_OBJ_ID].astype(str)

    DESIRED_COLUMN_ORDER = [
        "subject",
        "DB",
        "DB_Object_ID",
        "predicate",
        "object",
        "publications",
        "Evidence_Code",
        "supporting_objects",
        "annotation_date",
        "primary_knowledge_source",
        "aggregator",
        "protocol_id",
        "negated",
    ]

    for col in DESIRED_COLUMN_ORDER:
        if col not in df.columns:
            df[col] = None  # fill missing column with None
    return df


def merge_evidence_mapping(df: pd.DataFrame, evidence_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge GO annotation DataFrame with ECO evidence mapping.
    This function matches GO annotations with their corresponding evidence types based on both Evidence_Code and publications fields.

    """
    # Ensure DataFrames are copies to avoid modifying original data
    df = df.copy()
    evidence_df = evidence_df.copy()

    # Normalize evidence codes
    if df[PUBLICATIONS].apply(lambda x: isinstance(x, list)).any():
        df = df.explode(PUBLICATIONS)

    # If publications are NaN, drop those rows
    df = df[df[PUBLICATIONS].notna()]

    # Normalize case and strip whitespace
    df[PUBLICATIONS] = df[PUBLICATIONS].str.strip().str.upper()
    evidence_df[DB_REF] = evidence_df[DB_REF].str.strip().str.upper()
    evidence_df[EVIDENCE_CODE] = evidence_df[EVIDENCE_CODE].str.strip().str.upper()

    # Merge with ECO mapping
    merged = df.merge(
        evidence_df,
        how="left",
        left_on=[EVIDENCE_CODE, PUBLICATIONS],
        right_on=[EVIDENCE_CODE, DB_REF],
    )

    # Handle unmatched evidence types
    unmatched = merged[EVIDENCE_TYPE].isna()
    fallback = evidence_df[evidence_df[DB_REF] == "DEFAULT"].set_index(EVIDENCE_CODE)[EVIDENCE_TYPE]

    if unmatched.any():
        merged.loc[unmatched, EVIDENCE_TYPE] = merged.loc[unmatched, EVIDENCE_CODE].map(fallback)
    return merged.drop(columns=[DB_REF])


def process_go_annotations(input_path, output_path, eco_mapping_df, debug=False):
    """
    Full pipeline to process a GO annotation file and output a standardized, mapped result.
    - Loads a raw GAF file from the specified input path.
    - Transforms the data into a normalized format with semantic fields.
    - Merges the data with an ECO evidence mapping DataFrame.
    - Optionally prints debug information about intermediate steps.
    - Converts the annotation date to 'YYYY-MM-DD' format.
    - Writes the final result to the specified output CSV file.
    """
    try:
        raw_df = load_go_file(input_path)

        if debug:
            print("Raw data loaded:", raw_df.shape)
        transformed_df = transform_go_data(raw_df)

        if debug:
            print("Transformed data:", transformed_df.columns)
        merged = merge_evidence_mapping(transformed_df, eco_mapping_df)

        if merged.empty:
            print("Warning: Merging resulted in empty DataFrame.")
        else:
            merged[ANNOTATION_DATE] = pd.to_datetime(merged[ANNOTATION_DATE]).dt.strftime("%Y-%m-%d")
            merged.to_csv(output_path, index=False)
    except Exception as e:
        print(f"[ERROR] {e}")


# --- CLI Interface ---
@click.command()
@click.option("--input", "-i", required=True, help="Path to input GO annotation CSV file.")
@click.option("--output", "-o", required=True, help="Path to output normalized CSV file.")
def main(input, output):
    """CLI entry point to process GO annotations."""
    # Valid input file
    if not os.path.isfile(input):
        print(f"Error: Input file '{input}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Validate output directory
    output_dir = os.path.dirname(output) or "."
    if not os.path.isdir(output_dir):
        print(f"Error: Output directory '{output_dir}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Load the ECO mapping file once
    try:
        eco_mapping_df = pd.read_csv(ECO_MAPPING_URL, sep="\t", comment="#", header=None, names=ECO_MAPPING_COLUMNS)
    except Exception as e:
        click.echo(f"Error: Failed to load ECO mapping from {ECO_MAPPING_URL}: {str(e)}", err=True)
        sys.exit(1)

    # If all checks passed, then process annotations
    process_go_annotations(input, output, eco_mapping_df)


if __name__ == "__main__":
    main()
