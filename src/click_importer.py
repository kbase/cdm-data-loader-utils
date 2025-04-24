import pandas as pd
import click 


def process_go_annotations(input_path: str, output_path: str):
    """
    Processes GO annotations by cleaning, transforming, and mapping evidence types.

    Args:
        input_path (str): Path to the input GO annotations CSV file.
        output_path (str): Path to save the normalized output CSV file.
    """
    
    # load the input data 
    df = pd.read_csv(input_path)
    # df = pd.read_csv("annotations_data100.csv")
    
    # Drop unnecessary columns to simplify the DataFrame
    drop_columns = [
        "DB_Object_Symbol", "Aspect", "DB_Object_Name", "Synonym", 
        "DB_Object_Type", "Taxon", "Annotation_Extension", "Gene_Product_Form_ID"
        ]
    
    df = df.drop(columns = drop_columns)

    # Create a unified 'subject' identifier
    df["subject"] = df["DB"] + ":" + df["DB_Object_ID"]
    df = df.drop(columns=["DB", "DB_Object_ID"])

    # Reorder columns to place 'subject' first
    cols = ["subject"] + [col for col in df.columns if col != "subject"]
    df = df[cols]
    
    # Rename columns to match the desired schema
    df = df.rename(columns={
        "Qualifier": "predicate", 
        "GO_ID": "object",
        "DB_Reference": "publications", 
        "With_From": "supporting_objects",
        "Date": "annotation_date",
        "Assigned_By": "primary_knowledge_source"})
    
    # Convert 'annotation_date' to a standard format
    df["annotation_date"] = pd.to_datetime(df["annotation_date"], format='%Y%m%d').dt.strftime('%Y-%m-%d')
    df["annotation_date"] = df["annotation_date"].astype(str)
    
    # Add constant metadata columns 
    df["aggregator"] = "UniProt"
    df["protocol_id"] = "NULL"
    df["supporting_objects"] = df["supporting_objects"].fillna("NULL")
    
    # Define allowed predicates, determine negation, and clean predicate names
    allowed_predicates = ["enables", "contributes_to", "acts_upstream_of_or_within", "involved_in",
                          "acts_upstream_of", "acts_upstream_of_positive_effect",
                          "acts_upstream_of_negative_effect", "acts_upstream_of_or_within_negative_effect",
                          "acts_upstream_of_or_within_positive_effect", "located_in", "part_of", 
                          "is_active_in", "colocalizes_with"]

    df["negated"] = ~df["predicate"].isin(allowed_predicates)
    df["predicate"] = df["predicate"].str.replace(r"^NOT\|", "", regex=True)
    
    # Load the external ECO mapping data
    url = "https://raw.githubusercontent.com/evidenceontology/evidenceontology/master/gaf-eco-mapping.txt"
    mapping_columns = ["Evidence_Code", "DB_Reference", "evidence_type"]
    df2 = pd.read_csv(url, sep="\t", comment="#", header=None, names=mapping_columns)
    
    """
    merge two daraframe
    df1: go annotations data
    df2: ECO mapping data
    df2 DB_reference match with df1 publications
    df2 evidence code match with df1 evidence code 
    """

    # Normalize strings for consistent matching 
    df["publications"] = df["publications"].str.strip().str.upper()
    df2["DB_Reference"] = df2["DB_Reference"].str.strip().str.upper()
    df2["Evidence_Code"] = df2["Evidence_Code"].str.strip().str.upper()

    # Merge input data with ECO mapping data
    merged = df.merge(
        df2,how="left",
        left_on=["Evidence_Code", "publications"],
        right_on=["Evidence_Code", "DB_Reference"])
    
    # Handle unmatched revidence_type by assigning default mappings
    unmatched_mask = merged["evidence_type"].isna()
    default_mapping = df2[df2["DB_Reference"].str.upper() == "DEFAULT"]

    fallback_match = merged.loc[unmatched_mask, ["Evidence_Code"]].merge(default_mapping, how="left", on="Evidence_Code")
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

### python click_importer.py --input annotations_data100.csv --output normalized_annotation.csv ### 
