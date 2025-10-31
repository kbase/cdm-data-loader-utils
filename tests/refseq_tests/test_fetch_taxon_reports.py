import json
import pytest
import subprocess
from pyspark.sql import Row, SparkSession

### pytest -v refseq_tests/test_fetch_taxon_reports.py ### 


@pytest.fixture
def mock_parse_and_fetch(monkeypatch):
    # Mock fetch_reports_by_taxon → pretend API returns one fake report
    monkeypatch.setattr(
        "refseq_pipeline.core.datasets_api.fetch_reports_by_taxon",
        lambda taxon: [{"fake": "data"}]
    )

    # Mock parse_reports → return a tiny DataFrame
    def fake_parse_reports(reports, return_spark, spark: SparkSession):
        return spark.createDataFrame([Row(cdm_id="CDM001", taxid=None, n_contigs=10)])

    monkeypatch.setattr("refseq_pipeline.core.cdm_parse.parse_reports", fake_parse_reports)


@pytest.fixture
def taxid_file(tmp_path):
    path = tmp_path / "taxids.json"
    path.write_text(json.dumps(["11111", "22222"]))
    return str(path)


@pytest.mark.parametrize("prefer_spark", [True, False])
@pytest.mark.parametrize("mode", ["append", "overwrite"])

def test_cli_fetch_taxon_reports(tmp_path, monkeypatch, mock_parse_and_fetch, taxid_file, prefer_spark, mode):

    # Avoid writing delta tables during test
    monkeypatch.setattr("refseq_pipeline.core.spark_delta.write_delta_table", lambda *args, **kwargs: None)
    monkeypatch.setattr("refseq_pipeline.core.spark_delta.cleanup_after_write", lambda *args, **kwargs: None)

    monkeypatch.setenv("SPARK_SQL_WAREHOUSE_DIR", str(tmp_path / "warehouse"))

    cmd = [
        "python", "-m", "refseq_pipeline.cli.fetch_taxon_reports",
        "--database", "refseq_api",
        "--table", "assembly_stats",
        "--taxids-json", taxid_file,
        "--mode", mode,
    ]

    if prefer_spark:
        cmd.append("--prefer-spark")

    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, f"CLI failed:\n{result.stdout}\n{result.stderr}"

    assert "Summary" in result.stdout or result.stderr


