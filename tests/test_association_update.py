"""

Unit tests for the association_update module.

Run with:
    python3 -m pytest test_association_update.py

"""

import pytest
from pyspark.sql import SparkSession, Row
from pyspark.sql.functions import col
from pyspark.sql import functions as F
from delta import configure_spark_with_delta_pip

from association_update import (
    load_annotation, normalize_dates, process_predicates, add_metadata,
    load_eco_mapping, merge_evidence, write_output)

from association_update import (
    ANNOTATION_DATE, PREDICATE, NEGATED,
    AGGREGATOR, PROTOCOL_ID, SUBJECT, DB, DB_OBJ_ID,
    EVIDENCE_CODE, DB_REF, PUBLICATIONS, EVIDENCE_TYPE
)



# ---------------------------
# Spark Session Fixture
# ---------------------------
@pytest.fixture(scope="session")
def spark():
    builder = SparkSession.builder \
        .master("local[1]") \
        .appName("TestAssociationUpdate") \
        .config("spark.sql.warehouse.dir", "/tmp/spark-warehouse") \
        .config("spark.ui.enabled", "false") \
        .config("spark.driver.bindAddress", "127.0.0.1") \
        .config("spark.driver.host", "127.0.0.1") \
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension") \
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    return configure_spark_with_delta_pip(builder).getOrCreate()


# ---------------------------
# Test: Load Annotation
# ---------------------------
def test_load_annotation(spark, tmp_path):
    test_csv = tmp_path / "test_input.csv"
    test_csv.write_text("""DB,DB_Object_ID,Qualifier,GO_ID,DB_Reference,Evidence_Code,With_From,Date,Assigned_By
UniProtKB,P12345,enables,GO:0008150,PMID:123456,ECO:0000313,GO_REF:0000033,20240101,GO_Curator
""")

    df = load_annotation(spark, str(test_csv))
    row = df.collect()[0]

    assert row["predicate"] == "enables"
    assert row["object"] == "GO:0008150"
    assert row["publications"] == ["PMID:123456"]              
    assert row["supporting_objects"] == ["GO_REF:0000033"]     
    assert str(row["annotation_date"]) == "20240101"               
    assert row["primary_knowledge_source"] == "GO_Curator"


# ---------------------------
# Test: Normalize Dates
# ---------------------------
@pytest.mark.parametrize("date_input,expected", [
    ("20240101", "2024-01-01"),
    ("notadate", None)
])
def test_normalize_dates(spark, date_input, expected):
    df = spark.createDataFrame([(date_input,)], [ANNOTATION_DATE])
    result = normalize_dates(df).collect()[0][ANNOTATION_DATE]

    if expected is None:
        assert result is None
    else:
        assert str(result) == expected


# ---------------------------
# Test: Process Predicates
# ---------------------------
@pytest.mark.parametrize("predicate_val,expected_negated,expected_cleaned_predicate", [
    ("NOT|enables", True, "enables"),
    ("involved_in", False, "involved_in")
])
def test_process_predicates(spark, predicate_val, expected_negated, expected_cleaned_predicate):
    df = spark.createDataFrame([(predicate_val,)], ["Qualifier"])
    df = df.withColumn(PREDICATE, col("Qualifier"))
    result_df = process_predicates(df)
    row = result_df.collect()[0]

    assert row[NEGATED] == expected_negated
    assert row[PREDICATE] == expected_cleaned_predicate


# ---------------------------
# Test: Add Metadata
# ---------------------------
@pytest.mark.parametrize("db, db_obj_id, expected_subject", [
    ("UniProtKB", "P12345", "UniProtKB:P12345"),
    ("TAIR", "AT1G01010", "TAIR:AT1G01010"),
    ("MGI", "MGI:87938", "MGI:MGI:87938")
])
def test_add_metadata(spark, db, db_obj_id, expected_subject):
    df = spark.createDataFrame([(db, db_obj_id)], [DB, DB_OBJ_ID])
    result_df = add_metadata(df)
    row = result_df.collect()[0]

    assert row[AGGREGATOR] == "UniProt"
    assert row[PROTOCOL_ID] is None
    assert row[SUBJECT] == expected_subject


# ---------------------------
# Test: Load ECO Mapping
# ---------------------------
@pytest.mark.parametrize("eco_content,expected_rows", [
    ("ECO:0000313\tPMID:123456\tIEA\n", [("ECO:0000313", "PMID:123456", "IEA")]),
    ("ECO:0000256\tPMID:789012\tEXP\nECO:0000244\tDEFAULT\tTAS\n", [
        ("ECO:0000256", "PMID:789012", "EXP"),
        ("ECO:0000244", "DEFAULT", "TAS")
    ])
])
def test_load_eco_mapping_from_file(spark, tmp_path, eco_content, expected_rows):
    eco_file = tmp_path / "gaf-eco-mapping.txt"
    eco_file.write_text(eco_content)

    df = load_eco_mapping(spark, local_path=str(eco_file))
    result = [(row[EVIDENCE_CODE], row[DB_REF], row[EVIDENCE_TYPE]) for row in df.collect()]
    assert result == expected_rows


# ---------------------------
# Test: Merge Evidence Mapping
# ---------------------------
@pytest.mark.parametrize("annotation_rows,eco_rows,expected", [
    (
        # annotation df
        [Row(evidence_code="ECO:0000313", publications=["PMID:123456"])],
        # eco df
        [Row(evidence_code="ECO:0000313", db_ref="PMID:123456", evidence_type="IEA")],
        # expected result
        [("ECO:0000313", "PMID:123456", "IEA")]
    ),
    (
        # Fallback case
        [Row(evidence_code="ECO:0000256", publications=["PMID:999999"])],
        [Row(evidence_code="ECO:0000256", db_ref="DEFAULT", evidence_type="EXP")],
        [("ECO:0000256", "PMID:999999", "EXP")]
    )
])
def test_merge_evidence(spark, annotation_rows, eco_rows, expected):
    annotation_df = spark.createDataFrame(annotation_rows).select(
        F.col("evidence_code").alias(EVIDENCE_CODE),
        F.col("publications").alias(PUBLICATIONS)
    )
    eco_df = spark.createDataFrame(eco_rows).select(
        F.col("evidence_code").alias(EVIDENCE_CODE),
        F.col("db_ref").alias(DB_REF),
        F.col("evidence_type").alias(EVIDENCE_TYPE)
    )

    result_df = merge_evidence(annotation_df, eco_df)
    result = [(row[EVIDENCE_CODE], row[PUBLICATIONS], row[EVIDENCE_TYPE]) for row in result_df.collect()]
    assert result == expected


# ---------------------------
# Test: Delta Write + Read
# ---------------------------
def test_write_output_and_read_back(spark, tmp_path):

    # Sample test data
    data = [("GO:0008150", "UniProtKB", "2024-01-01")]
    columns = ["object", "db", "annotation_date"]
    df = spark.createDataFrame(data, columns)

    # Write to temporary Delta location
    output_path = str(tmp_path / "delta_table")
    write_output(df, output_path)

    # Read back and validate
    result_df = spark.read.format("delta").load(output_path)
    result = result_df.collect()

    assert len(result) == 1
    assert result[0]["object"] == "GO:0008150"
    assert result[0]["db"] == "UniProtKB"
    assert str(result[0]["annotation_date"]) == "2024-01-01"

