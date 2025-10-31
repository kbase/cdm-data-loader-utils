from refseq_pipeline.core.hashes_diff import diff_hash_and_get_changed_taxids
from refseq_pipeline.core.spark_delta import build_spark
from pyspark.sql import Row

### pytest -v refseq_tests/test_hashes_diff.py ###

def test_diff_hashes(tmp_path):
    spark = build_spark("test_db")

    # ---------------------------------------------------
    #  Create mock old + new snapshot entries
    # ---------------------------------------------------
    old_rows = [Row(accession="A1", kind="genome", content_sha256="sha_old", tag="20250101")]
    new_rows = [Row(accession="A1", kind="genome", content_sha256="sha_new", tag="20250201")]

    # ---------------------------------------------------
    #  Write to a real Delta directory
    # ---------------------------------------------------
    delta_dir = tmp_path / "refseq_api" / "assembly_hashes"
    spark.createDataFrame(old_rows).write.format("delta").mode("overwrite").save(str(delta_dir))
    spark.createDataFrame(new_rows).write.format("delta").mode("append").save(str(delta_dir))

    # ---------------------------------------------------
    #  Register the Delta table 
    # ---------------------------------------------------
    spark.sql("CREATE DATABASE IF NOT EXISTS refseq_api")
    spark.sql(f"""
        CREATE TABLE IF NOT EXISTS refseq_api.assembly_hashes
        USING DELTA
        LOCATION '{str(delta_dir)}'
    """)

    # ---------------------------------------------------
    #  Fake accession → taxid metadata
    # ---------------------------------------------------
    acc_index = {
        "A1": {"taxid": "99999"}   # matches function requirement: dict of dict
    }

    # ---------------------------------------------------
    #  Run the diff function
    # ---------------------------------------------------
    changed_taxids = diff_hash_and_get_changed_taxids(
        spark=spark,
        database="refseq_api",
        hash_table="assembly_hashes",
        acc_index=acc_index,
        tag_old="20250101",
        tag_new="20250201",
        debug=False
    )

    # ---------------------------------------------------
    #  Assert correctness
    # ---------------------------------------------------
    assert changed_taxids == ["99999"], f"Expected ['99999'], got {changed_taxids}"

    spark.stop()


    
