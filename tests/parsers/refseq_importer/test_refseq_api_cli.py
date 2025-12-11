from unittest.mock import MagicMock, patch

from click.testing import CliRunner

from cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli import (
    cli,
    main,
    parse_taxid_args,
)

# -------------------------------------------------
# test_parse_taxid_args
# -------------------------------------------------


def test_parse_taxid_args_basic() -> None:
    taxids = parse_taxid_args("123, 456x, abc789", None)
    assert taxids == ["123", "456", "789"]


def test_parse_taxid_args_file(tmp_path) -> None:
    file = tmp_path / "ids.txt"
    file.write_text("111\n222x\n333\n")

    taxids = parse_taxid_args("123", str(file))
    assert taxids == ["123", "111", "222", "333"]


# -------------------------------------------------
# test main()
# -------------------------------------------------


@patch("cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli.write_and_preview")
@patch("cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli.finalize_tables")
@patch("cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli.process_taxon")
@patch("cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli.write_delta")
@patch("cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli.build_spark")
def test_main_end_to_end(
    mock_build,
    mock_write_delta,
    mock_process,
    mock_finalize,
    mock_preview,
) -> None:
    mock_spark = MagicMock()
    mock_build.return_value = mock_spark
    mock_process.return_value = (["E"], ["C"], ["N"], ["I"])
    mock_finalize.return_value = ("EE", "CC", "NN", "II")

    main(
        taxid="123",
        api_key=None,
        database="refseq_api",
        mode="overwrite",
        debug=False,
        allow_genbank_date=False,
        unique_per_taxon=False,
        data_dir="/tmp",
    )

    mock_build.assert_called_once()
    mock_process.assert_called_once()
    mock_write_delta.assert_called()
    mock_preview.assert_called_once()


# -------------------------------------------------
# test cli() wrapper
# -------------------------------------------------


def test_cli_invocation() -> None:
    with patch("cdm_data_loader_utils.parsers.refseq_importer.cli.refseq_api_cli.main") as mock_main:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "--taxid",
                "123",
                "--data-dir",
                "/tmp",
            ],
        )
        assert result.exit_code == 0
        mock_main.assert_called_once()
