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
    python gene_association_file.py --input annotations_data100.csv --output normalized_annotation.csv

"""

import pandas as pd
import click 

# --- Load ECO mapping ---
ECO_MAPPING_URL = "http://purl.obolibrary.org/obo/eco/gaf-eco-mapping.txt"
mapping_columns = ["Evidence_Code", "DB_Reference", "evidence_type"]

try:
    eco_mapping_df = pd.read_csv(ECO_MAPPING_URL, sep="\t", comment="#", header=None, names=mapping_columns)
except Exception as e:
    raise ConnectionError(f"Failed to load ECO mapping from {ECO_MAPPING_URL}: {e}")

def process_go_annotations(input_path: str, output_path: str):

    """
    Processes GO annotations by cleaning, transforming, and mapping evidence types.

    Args:
        input_path (str): Path to the input GO annotations CSV file.
        output_path (str): Path to save the normalized output CSV file.
    """
    
    # load the input data 
    go_annotations_df = pd.read_csv(input_path)
    
    # Check if essential columns exist in the metadata 
    requeired_columns = {"DB","DB_Object_ID","Qualifier","GO_ID",
                         "DB_Reference","Evidence_Code","With_From",
                         "Date","Assigned_By"}
    
    if not requeired_columns.issubset(go_annotations_df.columns):
        missing_columns = requeired_columns - set(go_annotations_df.columns)
        raise ValueError(f"Missing required columns in the input file: {', '.join(missing_columns)}")

    
    # Drop unnecessary columns to simplify the DataFrame
    drop_columns = [
        "DB_Object_Symbol", "Aspect", "DB_Object_Name", "Synonym", 
        "DB_Object_Type", "Taxon", "Annotation_Extension", "Gene_Product_Form_ID"
        ]
    
    go_annotations_df = go_annotations_df.drop(columns=[col for col in drop_columns if col in go_annotations_df.columns])

    # Create a unified 'subject' identifier
    go_annotations_df["subject"] = go_annotations_df["DB"].astype(str) + ":" + go_annotations_df["DB_Object_ID"].astype(str)
    go_annotations_df = go_annotations_df.drop(columns=["DB", "DB_Object_ID"])

    # Reorder columns to place 'subject' first
    cols = ["subject"] + [col for col in go_annotations_df.columns if col != "subject"]
    go_annotations_df = go_annotations_df[cols]
    
    # Rename columns to match the desired schema
    go_annotations_df.rename(columns={"Qualifier": "predicate", 
                                      "GO_ID": "object",
                                      "DB_Reference": "publications", 
                                      "With_From": "supporting_objects",
                                      "Date": "annotation_date",
                                      "Assigned_By": "primary_knowledge_source"}, inplace=True)
    
    # Normalize date format and convert 'annotation_date' to a standard format
    try:
        go_annotations_df["annotation_date"] = pd.to_datetime(go_annotations_df["annotation_date"], format='%Y%m%d').dt.strftime('%Y-%m-%d')
    except Exception as e:
        raise ValueError(f"Error parsing annotation_date: {e}")
    
    # Ensure 'annotation_date' is a string
    go_annotations_df["annotation_date"] = go_annotations_df["annotation_date"].astype(str)
    
    # Add constant metadata columns 
    go_annotations_df["aggregator"] = "UniProt"
    go_annotations_df["protocol_id"] = None

    # Split multi-value fields into lists
    go_annotations_df["publications"] = go_annotations_df["publications"].apply(lambda x: str(x).split('|') if pd.notna(x) else None)
    go_annotations_df["supporting_objects"] = go_annotations_df["supporting_objects"].apply(lambda x: str(x).split('|') if pd.notna(x) else None)

    # Detect negation
    go_annotations_df["negated"] = go_annotations_df["predicate"].str.contains(r"^NOT\|", regex=True)

    # Clean predicate string
    go_annotations_df["predicate"] = go_annotations_df["predicate"].str.replace(r"^NOT\|", "", regex=True)


    # Define allowed predicates, determine negation, and clean predicate names
    allowed_predicates = ["enables", "contributes_to", "acts_upstream_of_or_within", "involved_in",
                          "acts_upstream_of", "acts_upstream_of_positive_effect",
                          "acts_upstream_of_negative_effect", "acts_upstream_of_or_within_negative_effect",
                          "acts_upstream_of_or_within_positive_effect", "located_in", "part_of", 
                          "is_active_in", "colocalizes_with"]

    invalid_predicates = ~go_annotations_df["predicate"].isin(allowed_predicates)
    if invalid_predicates.any():
        raise ValueError(f"Invalid predicates found: {go_annotations_df.loc[invalid_predicates, 'predicate'].unique()}")

    # Normalize strings for consistent matching 
    df_exploded = go_annotations_df.explode("publications")
    df_exploded = df_exploded[df_exploded["publications"].notna()]
    df_exploded["publications"] = df_exploded["publications"].str.strip().str.upper()

    eco_mapping_df["DB_Reference"] = eco_mapping_df["DB_Reference"].str.strip().str.upper()
    eco_mapping_df["Evidence_Code"] = eco_mapping_df["Evidence_Code"].str.strip().str.upper()


    """
    merge two dataframes
    go_annotations_df: go annotations data
    eco_mapping_df: ECO mapping data
    eco_mapping_df DB_reference match with publications
    eco_mapping_df evidence code match with go_annotations_df evidence code 
    """

    # Merge input data with ECO mapping data
    merged = df_exploded.merge(
        eco_mapping_df, how="left",
        left_on=["Evidence_Code", "publications"],
        right_on=["Evidence_Code", "DB_Reference"])
    
    # Handle unmatched revidence_type by assigning default mappings
    unmatched_mask = merged["evidence_type"].isna()
    fallback_df = eco_mapping_df[eco_mapping_df["DB_Reference"] == "DEFAULT"]

    fallback_match = merged.loc[unmatched_mask, ["Evidence_Code"]].merge(fallback_df, how="left", on="Evidence_Code")
    merged.loc[unmatched_mask, "evidence_type"] = fallback_match["evidence_type"].values

    # Drop the unnecessary merged columns 
    merged = merged.drop(columns=["DB_Reference"])
    
    # Save the processed annotations to a new CSV file
    merged.to_csv(output_path, index=False)


@click.command()
@click.option('--input', '-i', required=True, help='Path to input GO annotation CSV file.')
@click.option('--output', '-o', required=True, help='Path to output normalized CSV file.')


def main(input, output):
    """
    CLI entry point to process GO annotations.

    Args:
        input (str): Input CSV file path.
        output (str): Output CSV file path.
    """

    process_go_annotations(input, output)


if __name__ == '__main__':
    main()

