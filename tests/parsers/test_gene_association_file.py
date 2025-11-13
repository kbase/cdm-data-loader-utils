import pandas as pd
import pytest
from cdm_data_loader_utils.parsers.gene_association_file import (
    process_go_annotations,
    ASSOCIATION_COL_TYPES,
    merge_evidence_mapping,
    normalize_dates,
    process_predicates,
)


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
    return pd.DataFrame(
        [
            {
                "Evidence_Code": "ISS",
                "DB_Reference": "GO_REF:0000024",
                "evidence_type": "ECO:0000250",
            },
            {"Evidence_Code": "ND", "DB_Reference": "DEFAULT", "evidence_type": "ECO:0000307"},
            {
                "Evidence_Code": "IDA",
                "DB_Reference": "PMID:12345678",
                "evidence_type": "ECO:0000314",
            },
            {
                "Evidence_Code": "EXP",
                "DB_Reference": "PMID:87654321",
                "evidence_type": "ECO:0000269",
            },
            {
                "Evidence_Code": "IMP",
                "DB_Reference": "GO_REF:0000033",
                "evidence_type": "ECO:0000315",
            },
        ]
    )


@pytest.fixture
def annotation_df():
    """Fixture providing an example annotations DataFrame for testing evidence mapping."""
    return pd.DataFrame(
        [
            {"Evidence_Code": "ISS", "publications": "GO_REF:0000024"},  # Direct match
            {"Evidence_Code": "ND", "publications": "PMID:UNKNOWN"},  # Fallback match
            {"Evidence_Code": "EXP", "publications": "GO_REF:0009999"},  # No match, no fallback
            {
                "Evidence_Code": "IDA",
                "publications": "GO_REF:0000001|PMID:123456",
            },  # Multiple publications (explode test)
            {"Evidence_Code": "IMP", "publications": None},  # Null publications
        ]
    )


# --- Tests ---
def test_go_annotation_processed(temp_csv_files, evidence_mapping_df):
    input_path, output_path = temp_csv_files
    process_go_annotations(str(input_path), str(output_path), evidence_mapping_df)
    assert output_path.exists(), "Output file not created."
    df = pd.read_csv(output_path)
    expected_cols = set(ASSOCIATION_COL_TYPES.keys())
    assert expected_cols.issubset(set(df.columns)), (
        f"Output columns do not match expected.\nExpected: {expected_cols}\nGot: {set(df.columns)}"
    )


def test_io_behavior(temp_csv_files, evidence_mapping_df):
    """Integration test for I/O behavior of process_go_annotations."""
    input_path, output_path = temp_csv_files
    assert not output_path.exists(), "Output file already exist before processing."
    process_go_annotations(str(input_path), str(output_path), evidence_mapping_df)
    assert output_path.exists(), "Output file not created."
    df_out = pd.read_csv(output_path)
    assert not df_out.empty, "Output DataFrame is empty."
    for col in ASSOCIATION_COL_TYPES.keys():
        assert col in df_out.columns, f"Column {col} missing in output."


def test_column_types(temp_csv_files, evidence_mapping_df):
    input_path, output_path = temp_csv_files
    process_go_annotations(str(input_path), str(output_path), evidence_mapping_df)
    df = pd.read_csv(output_path)
    for column, expected_type in ASSOCIATION_COL_TYPES.items():
        not_null = df[column].dropna()
        if not not_null.empty:
            mismatches = ~not_null.map(lambda x: isinstance(x, expected_type))
            assert not mismatches.any(), f"Column {column} has incorrect types."


def test_negated_logic(temp_csv_files, evidence_mapping_df):
    input_path, output_path = temp_csv_files
    process_go_annotations(str(input_path), str(output_path), evidence_mapping_df)
    df = pd.read_csv(output_path)
    assert df["negated"].any(), "Expected at least one negated is True row, it present NOT|enables"
    assert (~df["negated"]).any(), "Expected at least one negated is False row, it present enables"


# --- Helper function for evidence mapping test ---


def test_merge_and_fallback(annotation_df, evidence_mapping_df):
    """Test the evidence mapping merge and fallback logic."""
    merged = merge_evidence_mapping(annotation_df, evidence_mapping_df)
    assert merged.loc[0, "evidence_type"] == "ECO:0000250", "Direct match failed."
    assert merged.loc[1, "evidence_type"] == "ECO:0000307", "Fallback match failed."
    assert "DB_Reference" not in merged.columns, "'DB_Reference' should have been dropped after merge."


def test_no_fallback():
    annotation_df = pd.DataFrame([{"Evidence_Code": "IEA", "publications": "PMID:999999"}])
    evidence_df = pd.DataFrame(
        [{"Evidence_Code": "EXP", "DB_Reference": "GO_REF:0000024", "evidence_type": "ECO:0000269"}]
    )
    merged = merge_evidence_mapping(annotation_df, evidence_df)
    assert pd.isna(merged.loc[0, "evidence_type"]), (
        "Expected evidence_type to be NaN when no match and no fallback are available."
    )


def test_fallback_after_explode():
    annotation_df = pd.DataFrame([{"Evidence_Code": "ND", "publications": ["PMID:111111", "PMID:222222"]}])
    evidence_df = pd.DataFrame([{"Evidence_Code": "ND", "DB_Reference": "DEFAULT", "evidence_type": "ECO:0000307"}])
    merged = merge_evidence_mapping(annotation_df, evidence_df)
    assert len(merged) == 2, "Expected two rows after exploding publication list"
    assert all(merged["evidence_type"] == "ECO:0000307"), "Fallback evidence_type not correctly applied after explode"


def test_normalize_dates():
    df = pd.DataFrame({"annotation_date": ["20200101", "20221301", "abcd1234", "20181231", None]})
    result = normalize_dates(df.copy())
    expected = ["2020-01-01", None, None, "2018-12-31", None]
    actual = result["annotation_date"].tolist()
    for a, e in zip(actual, expected):
        if e is None:
            assert pd.isna(a)
        else:
            assert a == e


def test_process_predicates():
    df = pd.DataFrame({"predicate": ["NOT|enables", "involved_in", "NOT|located_in", "part_of"]})
    result = process_predicates(df.copy())
    assert result["negated"].tolist() == [True, False, True, False], "Negation detection incorrect."
    assert result["predicate"].tolist() == ["enables", "involved_in", "located_in", "part_of"], (
        "Predicate cleaning incorrect."
    )
