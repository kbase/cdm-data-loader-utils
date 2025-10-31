import json
import pytest
from click.testing import CliRunner
from refseq_pipeline.cli.diff_changed_taxids import main

### pytest -v refseq_tests/test_diff_changed_taxids.py ###


@pytest.mark.parametrize(
    "fake_taxids, expect_write",
    [
        (["11111", "22222"], True),  
        ([], False),                  
    ],
)
@pytest.mark.parametrize(
    "subdir, filename",
    [
        ("", "changed_taxids.json"),                   
        ("nested/output", "changed_taxids_custom.json")  
    ],
)
def test_cli_with_out_path(monkeypatch, tmp_path, fake_taxids, expect_write, subdir, filename):

    monkeypatch.setattr(
        "refseq_pipeline.cli.diff_changed_taxids.run_diff_changed_taxids",
        lambda database, hash_table, new_tag, old_tag: fake_taxids
    )

    out_dir = tmp_path / subdir if subdir else tmp_path
    out_path = out_dir / filename

    runner = CliRunner()
    args = [
        "--database", "refseq_api",
        "--hash-table", "assembly_hashes",
        "--old-tag", "20250101",
        "--new-tag", "20250201",
        "--out-path", str(out_path),
    ]
    result = runner.invoke(main, args)

    assert result.exit_code == 0, f"CLI failed:\n{result.output}"

    if expect_write:
        assert out_path.exists(), f"Expected output not created: {out_path}"
        data = json.loads(out_path.read_text())
        assert data == fake_taxids, f"JSON content mismatch: {data} != {fake_taxids}"
        assert f"{len(fake_taxids)} changed TaxIDs written" in result.output
    else:
        assert not out_path.exists(), f"Unexpected file created: {out_path}"
        assert "No changed TaxIDs found between 20250101 and 20250201" in result.output


        