import os
import pandas as pd
from importer import process_go_annotations
import pytest 


def test_go_annotation_processed():
    assert os.path.exists("new_annotations1.csv"), "Output file not found" ## csv exists 
    df = pd.read_csv("new_annotations1.csv") 
    assert "subject" in df.columns
    assert "predicate" in df.columns
    assert "object" in df.columns 
    assert "publications" in df.columns 
    assert "Evidence_Code" in df.columns 
    assert "supporting_objects" in df.columns
    assert "annotation_date" in df.columns
    assert "primary_knowledge_source" in df.columns 
    assert "aggregator" in df.columns
    assert "protocol_id" in df.columns 
    assert "negated" in df.columns 
    assert "evidence_type" in df.columns 

## I/O Behavior 
example_csv = """
subject	predicate	object	publications	Evidence_Code	supporting_objects	annotation_date	primary_knowledge_source	aggregator	protocol_id	negated	evidence_type
UniProtKB:Q04604	enables	GO:0005507	GO_REF:0000024	ISS	UniProtKB:Q9ZP19	2016-11-18	UniProt	UniProt	NULL	FALSE	ECO:0000250
"""

@pytest.fixture
def temp_csv_files(tmp_path):
    input_path = tmp_path / "go_annotations.csv"
    output_path = tmp_path / "new_annotations1.csv"

    # Write the input CSV to a file
    input_path.write_text(example_csv)
    return input_path, output_path


def test_io_behavior(temp_csv_files):
    input_path, output_path = temp_csv_files

    # Run the function
    process_go_annotations(str(input_path), str(output_path))

    # Ensure output file was created
    assert output_path.exists()

    # Read the output and confirm it's not empty
    df_out = pd.read_csv(output_path)
    assert not df_out.empty

    # Check all the expected columns exist or not 
    expected_columns = ["subject", "predicate", "object", "publications", "Evidence_Code", 
                        "supporting_objects", "annotation_date", "primary_knowledge_source",
                        "aggregator", "protocol_id", "negated", "evidence_type"]
    
    for col in expected_columns:
        assert col in df_out.columns


## check data types and column names 
def test_column_types(tmp_path):
    output_file = tmp_path / "test_output.csv"

    # Call function
    process_go_annotations(
        input_path="go_annotation.csv",  
        output_path=str(output_file)
    )

    df = pd.read_csv(output_file)

    # Define expected types
    expected_types = {
        "subject": str,
        "predicate": str,
        "object": str,
        "publications": str,
        "supporting_objects": str,
        "annotation_date": str,  
        "primary_knowledge_source": str,
        "aggregator": str,
        "protocol_id": str,
        "evidence_type": str,
        "negated": bool
    }

  
    for column, expected_type in expected_types.items():
        invalid = df[column].dropna().apply(lambda x: not isinstance(x, expected_type))
        assert not invalid.any(), f"Column {column} has values not of type {expected_type.__name__}"


## Negation logic test 
def test_negated_logic(temp_csv_files):
    input_path, output_path = temp_csv_files
    process_go_annotations(str(input_path), str(output_path))

    df = pd.read_csv(output_path)
    assert df["negated"].iloc[0] == False, "Expected predicate to be allowed"


#######################################################################

## Copy the code from importer.py 
## Does it properly merge evidence_type from the mapping file?
## Are fallback values applied when no direct match is found?

#######################################################################

def merge_evidence_mapping(df, df2):
    df["publications"] = df["publications"].str.strip().str.upper()
    df2["DB_Reference"] = df2["DB_Reference"].str.strip().str.upper()
    df2["Evidence_Code"] = df2["Evidence_Code"].str.strip().str.upper()

    merged = df.merge(
        df2,
        how="left",
        left_on=["Evidence_Code", "publications"],
        right_on=["Evidence_Code", "DB_Reference"]
    )

    unmatched_mask = merged["evidence_type"].isna()
    fallback_df = df2[df2["DB_Reference"] == "DEFAULT"]

    fallback_match = merged.loc[unmatched_mask, ["Evidence_Code"]].merge(
        fallback_df,
        how="left",
        on="Evidence_Code"
    )

    merged.loc[unmatched_mask, "evidence_type"] = fallback_match["evidence_type"].values
    return merged.drop(columns=["DB_Reference"])


## Test the partial code 
def test_evidence_mapping_merge_and_fallback():
    df = pd.DataFrame([
        {"Evidence_Code": "ISS", "publications": "GO_REF:0000024"},
        {"Evidence_Code": "ND", "publications": "PMID:UNKNOWN"},
    ])

    mapping_df = pd.DataFrame([
        {"Evidence_Code": "ISS", "DB_Reference": "GO_REF:0000024", "evidence_type": "ECO:0000250"},
        {"Evidence_Code": "ND", "DB_Reference": "DEFAULT", "evidence_type": "ECO:0000307"},
    ])
   
    merged = merge_evidence_mapping(df, mapping_df)

    assert merged.loc[0, "evidence_type"] == "ECO:0000250"
    assert merged.loc[1, "evidence_type"] == "ECO:0000307"
    assert "DB_Reference" not in merged.columns

