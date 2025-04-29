import pandas as pd
import pytest
from parsers.gene_association_file import (
    process_go_annotations, 
    ASSOCIATION_COL_TYPES, 
    ASSOCIATION_COL_HEADERS,
    merge_evidence_mapping, )


# --- Fixtures ---

@pytest.fixture
def temp_csv_files(tmp_path):
    """Create temporary CSV files for input and output testing."""

    input_path = tmp_path / "annotations_data100.csv"
    output_path = tmp_path / "normalized_annotation.csv"

    example_csv = (
        "DB,DB_Object_ID,DB_Object_Symbol,Qualifier,GO_ID,DB_Reference,Evidence_Code,With_From,Aspect,DB_Object_Name,Synonym,DB_Object_Type,Taxon,Date,Assigned_By,Annotation_Extension,Gene_Product_Form_ID\n"
        "UniProtKB,Q04604,TYR,enables,GO:0005507,GO_REF:0000024,ISS,UniProtKB:Q9ZP19,F,Tyrosinase,TYR|TYRS,protein,taxon:8409,20161118,UniProt,,\n"
        "UniProtKB,Q97823,ABC,NOT|enables,GO:0001234,GO_REF:0000001,ISS,UniProtKB:Q9ZP19,F,ABC_Protein,,protein,taxon:9606,20170101,UniProt,,\n"
        "UniProtKB,Q88219,XYZ,located_in,GO:0044464,,EXP,,F,XYZ Protein,,protein,taxon:9606,20210101,UniProt,,\n"
        "UniProtKB,Q79372,FOO,enables,GO:0000001,GO_REF:1234567,IMP,UniProtKB:P12345,F,FOO Protein,,protein,taxon:9606,20190101,UniProt,,\n"
        )

    input_path.write_text(example_csv)
    return input_path, output_path


@pytest.fixture
def evidence_mapping_df():
    """Example evidence mapping DataFrame for testing."""
    return pd.DataFrame([
        {"Evidence_Code": "ISS", "DB_Reference": "GO_REF:0000024", "evidence_type": "ECO:0000250"},
        {"Evidence_Code": "ND", "DB_Reference": "DEFAULT", "evidence_type": "ECO:0000307"},
    ])


@pytest.fixture
def annotation_df():
    """Fixture providing an example annotations DataFrame for testing evidence mapping."""
    return pd.DataFrame([
        {"Evidence_Code": "ISS", "publications": "GO_REF:0000024"}, # Direct match
        {"Evidence_Code": "ND", "publications": "PMID:UNKNOWN"}, # Fallback match
        {"Evidence_Code": "EXP", "publications": "GO_REF:0009999"}, # No match, no fallback
        {"Evidence_Code": "IDA", "publications": "GO_REF:0000001|PMID:123456"}, # Multiple publications (explode test)
        {"Evidence_Code": "IMP", "publications": None}, # Null publications
    ])


# --- Tests ---

@pytest.mark.unit
def test_go_annotation_processed(temp_csv_files):
    input_path, output_path = temp_csv_files
    process_go_annotations(str(input_path), str(output_path))

    assert output_path.exists(), "Output file not created."
    df = pd.read_csv(output_path)

    assert set(df.columns) == set(ASSOCIATION_COL_HEADERS), (
        f"Output DataFrame columns do not match expected columns.\nExpected: {ASSOCIATION_COL_HEADERS}\nGot: {list(df.columns)}"
    )


@pytest.mark.integration
def test_io_behavior(temp_csv_files):
    """Integration test for I/O behavior of process_go_annotations."""
    input_path, output_path = temp_csv_files

    # Before process_go_annotations runs, make sure that output_path does not exist # 
    assert not output_path.exists(), "Output file already exist before processing."
    process_go_annotations(str(input_path), str(output_path))

    assert output_path.exists(), "Output file not created."
    df_out = pd.read_csv(output_path)

    assert not df_out.empty, "Output DataFrame is empty."

    for col in ASSOCIATION_COL_HEADERS:
        assert col in df_out.columns, f"Column {col} missing in output."


@pytest.mark.unit
def test_column_types(temp_csv_files):
    """Test the column data types of the processed annotations."""
    input_path, output_path = temp_csv_files

    process_go_annotations(str(input_path), str(output_path))
    df = pd.read_csv(output_path)

    for column, expected_type in ASSOCIATION_COL_TYPES.items():
        invalid = df[column].dropna().apply(lambda x: not isinstance(x, expected_type))
        assert not invalid.any(), f"Column {column} has invalid data types."


@pytest.mark.unit
def test_negated_logic(temp_csv_files):
    """Test negation logic of predicate column."""
    input_path, output_path = temp_csv_files

    process_go_annotations(str(input_path), str(output_path))
    df = pd.read_csv(output_path)

     # At least one negated == False (normal enables)
    assert (df["negated"] == False).any(), "Expected at least one False negation not found."

    # At least one negated == True (NOT| enables)
    assert (df["negated"] == True).any(), "Expected at least one True negation not found."


# --- Helper function for evidence mapping test --- 

@pytest.mark.unit
def test_evidence_mapping_merge_and_fallback(annotation_df, evidence_mapping_df):
    """Test the evidence mapping merge and fallback logic."""
    merged = merge_evidence_mapping(annotation_df, evidence_mapping_df)

    assert merged.loc[0, "evidence_type"] == "ECO:0000250", "Direct match failed."
    assert merged.loc[1, "evidence_type"] == "ECO:0000307", "Fallback match failed."
    assert "DB_Reference" not in merged.columns, "Unnecessary column 'DB_Reference' present."


### PYTHONPATH=. pytest tests/test_gene_association_file.py ###
