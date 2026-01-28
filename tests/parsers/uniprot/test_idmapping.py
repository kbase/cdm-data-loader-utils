"""Tests for UniProt ID Mapping importer."""

import datetime
from pathlib import Path

import pytest
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as sf

from cdm_data_loader_utils.audit.schema import (
    AUDIT_SCHEMA,
    METRICS,
    N_INVALID,
    N_READ,
    N_VALID,
    PIPELINE,
    REJECTS,
    RUN_ID,
    SOURCE,
    VALIDATION_ERRORS,
)
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.parsers.uniprot.idmapping import ingest
from cdm_data_loader_utils.utils.spark_delta import write_delta
from tests.audit.conftest import (
    INIT_TIMESTAMP_FIELDS,
)
from tests.helpers import create_empty_delta_table

SAVE_DIR = "spark.sql.warehouse.dir"


@pytest.mark.requires_spark
@pytest.mark.parametrize("file_name", ["ECOLI_83333_idmapping.txt", "ECOLI_83333_idmapping.tsv.gz"])
@pytest.mark.parametrize("parse_only", [True, False])
def test_ingest(
    test_data_dir: Path,
    pipeline_run: PipelineRun,
    spark: SparkSession,
    file_name: str,
    parse_only: bool,
) -> None:
    """Test basic ingestion of UniProt ID Mapping data."""
    # prep the appropriate tables
    for t in AUDIT_SCHEMA:
        create_empty_delta_table(spark, pipeline_run.namespace, t, AUDIT_SCHEMA[t])

    test_mapping = test_data_dir / "uniprot" / "idmapping" / file_name

    if parse_only:
        df = ingest(spark, pipeline_run, str(test_mapping))
    else:
        table_name = "uniprot_identifier"
        # the test will write the delta table; we will check the results
        write_delta(spark, ingest(spark, pipeline_run, str(test_mapping)), pipeline_run.namespace, table_name, "append")
        df = spark.table(f"{pipeline_run.namespace}.{table_name}")

    assert isinstance(df, DataFrame)
    assert df.columns == ["uniprot_id", "db", "xref", "description", "source", "relationship"]
    n_rows = 331510
    assert df.count() == n_rows
    assert df.select("uniprot_id").filter(~sf.col("uniprot_id").startswith("UniProt:")).count() == 0

    # a few sanity checks on the data
    n_entries = 4598
    # 4598 NCBI_TaxID, UniProtKB-ID, CRC64
    for db in ["NCBI_TaxID", "CRC64", "UniProtKB-ID"]:
        assert df.filter(df.db == db).count() == n_entries
    # there are 4625 unique IDs, including uniparc isoforms
    n_unique_ids = 4625
    assert df.select("uniprot_id").distinct().count() == n_unique_ids
    assert df.filter(sf.coalesce(sf.col("db"), sf.lit("")) == "").count() == 0
    assert df.filter(sf.coalesce(sf.col("xref"), sf.lit("")) == "").count() == 0

    # count the number of IDs without a hyphen
    assert df.select("uniprot_id").filter(~sf.col("uniprot_id").contains("-")).distinct().count() == n_entries
    n_isoforms = n_unique_ids - n_entries
    # and the number with a hyphen
    assert df.select("uniprot_id").filter(sf.col("uniprot_id").contains("-")).distinct().count() == n_isoforms

    n_entities_without_uniprot = (
        df.groupBy("uniprot_id")
        .agg(sf.sum(sf.when(sf.col("db") == "UniProtKB-ID", 1).otherwise(0)).alias("has_uniprot"))
        .filter(sf.col("has_uniprot") == 0)
        .select("uniprot_id")
    ).count()
    # the UniParc isoforms do not have UniProt IDs
    assert n_entities_without_uniprot == n_isoforms

    # check the audit tables
    rejects_df = spark.table(f"{pipeline_run.namespace}.{REJECTS}")
    assert rejects_df.count() == 0
    metrics_df = spark.table(f"{pipeline_run.namespace}.{METRICS}")
    assert metrics_df.count() == 1
    row = metrics_df.collect()[0].asDict()
    # filter out timestamp rows
    for f in INIT_TIMESTAMP_FIELDS.get(METRICS, []):
        assert isinstance(row[f], datetime.datetime)
        del row[f]
    assert row == {
        N_READ: n_rows,
        N_VALID: n_rows,
        N_INVALID: 0,
        VALIDATION_ERRORS: [],
        PIPELINE: pipeline_run.pipeline,
        RUN_ID: pipeline_run.run_id,
        SOURCE: pipeline_run.source_path,
    }
