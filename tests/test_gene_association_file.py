import pandas as pd
import pytest
from gene_association_file import process_go_annotations


@pytest.fixture
def temp_csv_files(tmp_path):
    """Create temporary CSV files for input and output testing."""

    input_path = tmp_path / "annotations_data100.csv"
    output_path = tmp_path / "normalized_annotation.csv"

    example_csv = """DB,DB_Object_ID,DB_Object_Symbol,Qualifier,GO_ID,DB_Reference,Evidence_Code,With_From,Aspect,DB_Object_Name,Synonym,DB_Object_Type,Taxon,Date,Assigned_By,Annotation_Extension,Gene_Product_Form_ID
UniProtKB,Q04604,TYR,enables,GO:0005507,GO_REF:0000024,ISS,UniProtKB:Q9ZP19,F,Tyrosinase,TYR|TYRS,protein,taxon:8409,20161118,UniProt"""

    input_path.write_text(example_csv)
    return input_path, output_path


@pytest.fixture
def evidence_mapping_df():
    """Fixture providing an example mapping DataFrame for testing."""
    return pd.DataFrame([
        {"Evidence_Code": "ISS", "DB_Reference": "GO_REF:0000024", "evidence_type": "ECO:0000250"},
        {"Evidence_Code": "ND", "DB_Reference": "DEFAULT", "evidence_type": "ECO:0000307"},
    ])


@pytest.fixture
def annotation_df():
    """Fixture providing an example annotations DataFrame for testing evidence mapping."""
    return pd.DataFrame([
        {"Evidence_Code": "ISS", "publications": "GO_REF:0000024"},
        {"Evidence_Code": "ND", "publications": "PMID:UNKNOWN"},
    ])


@pytest.mark.unit
def test_go_annotation_processed(tmp_path):
    """Test that the processed GO annotation CSV file has correct columns."""
    output_file = tmp_path / "new_annotations1.csv"
    process_go_annotations("annotations_data100.csv", str(output_file))

    assert output_file.exists(), "Output file not created."
    df = pd.read_csv(output_file)

    expected_columns = [
        "subject", "predicate", "object", "publications", "Evidence_Code", "supporting_objects",
        "annotation_date", "primary_knowledge_source", "aggregator", "protocol_id", "negated", "evidence_type"
    ]

    for col in expected_columns:
        assert col in df.columns, f"Column {col} missing in output."


@pytest.mark.integration
def test_io_behavior(temp_csv_files):
    """Integration test for I/O behavior of process_go_annotations."""
    input_path, output_path = temp_csv_files

    process_go_annotations(str(input_path), str(output_path))

    assert output_path.exists(), "Output file not created."
    df_out = pd.read_csv(output_path)

    assert not df_out.empty, "Output DataFrame is empty."

    expected_columns = [
        "subject", "predicate", "object", "publications", "Evidence_Code",
        "supporting_objects", "annotation_date", "primary_knowledge_source",
        "aggregator", "protocol_id", "negated", "evidence_type"
    ]

    for col in expected_columns:
        assert col in df_out.columns, f"Column {col} missing in output."


@pytest.mark.unit
def test_column_types(temp_csv_files):
    """Test the column data types of the processed annotations."""
    input_path, output_path = temp_csv_files

    process_go_annotations(str(input_path), str(output_path))
    df = pd.read_csv(output_path)

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
        assert not invalid.any(), f"Column {column} has invalid data types."


@pytest.mark.unit
def test_negated_logic(temp_csv_files):
    """Test negation logic of predicate column."""
    input_path, output_path = temp_csv_files

    process_go_annotations(str(input_path), str(output_path))
    df = pd.read_csv(output_path)

    assert df.loc[0, "negated"] == False, "Negation logic failed for allowed predicate."


@pytest.mark.unit
def merge_evidence_mapping(df, df2):
    """Merge annotations DataFrame with evidence mapping DataFrame."""
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
        fallback_df, how="left", on="Evidence_Code"
    )

    merged.loc[unmatched_mask, "evidence_type"] = fallback_match["evidence_type"].values
    return merged.drop(columns=["DB_Reference"])


@pytest.mark.unit
def test_evidence_mapping_merge_and_fallback(annotation_df, evidence_mapping_df):
    """Test the evidence mapping merge and fallback logic."""
    merged = merge_evidence_mapping(annotation_df, evidence_mapping_df)

    assert merged.loc[0, "evidence_type"] == "ECO:0000250", "Direct match failed."
    assert merged.loc[1, "evidence_type"] == "ECO:0000307", "Fallback match failed."
    assert "DB_Reference" not in merged.columns, "Unnecessary column 'DB_Reference' present."


### PYTHONPATH=. pytest test_gene_association_file.py ###
