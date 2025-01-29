"""Tests for the MultiGenomeDataFileCreator."""

from pathlib import Path

import pytest

from src.parsers.genome_loader import MultiGenomeDataFileCreator


@pytest.mark.parametrize(
    "genome_paths_file",
    [
        [None, "Missing genome_paths_file"],
        ["", "Missing genome_paths_file"],
        ["    ", "Missing genome_paths_file"],
        ["does/not/exist", "genome_paths_file 'does/not/exist' does not exist"],
        ["tests/data/genome_paths_file/valid.json"],
    ],
)
@pytest.mark.parametrize(
    "output_dir",
    [
        [None, "Missing output_dir"],
        ["", "Missing output_dir"],
        ["    ", "Missing output_dir"],
        ["valid_dir"],
    ],
)
@pytest.mark.parametrize(
    "checkm2_path",
    [
        [None, None],
        ["", None],
        ["    ", None],
        ["valid_dir", "valid_dir"],
        [" some valid_dir   \t", "some valid_dir"],
    ],
)
@pytest.mark.parametrize(
    "stats_path",
    [
        [None, None],
        ["", None],
        ["    ", None],
        ["valid_dir", "valid_dir"],
        [" some valid_dir   \t", "some valid_dir"],
    ],
)
def test_init_multigenomedatafilecreator(
    genome_paths_file: list[str | None],
    output_dir: list[str | None],
    checkm2_path: list[str | None],
    stats_path: list[str | None],
) -> None:
    """Test that MGDFC init fails without required params."""
    errs = []
    if len(genome_paths_file) == 2:
        errs.append(genome_paths_file[1])
    if len(output_dir) == 2:
        errs.append(output_dir[1])

    params = [x[0] for x in [genome_paths_file, output_dir, checkm2_path, stats_path]]

    if errs:
        err_msg = "MultiGenomeDataFileCreator init error:\n" + "\n".join(errs)
        with pytest.raises(RuntimeError, match=err_msg):
            MultiGenomeDataFileCreator(*params)
    else:
        mgdfc = MultiGenomeDataFileCreator(*params)
        assert mgdfc.genome_paths_file == Path(genome_paths_file[0])
        assert mgdfc.output_dir == Path(output_dir[0])
        if checkm2_path[1] is None:
            assert mgdfc.checkm2_file is None
        else:
            assert mgdfc.checkm2_file == Path(checkm2_path[1])
        if stats_path[1] is None:
            assert mgdfc.stats_file is None
        else:
            assert mgdfc.stats_file == Path(stats_path[1])


@pytest.mark.parametrize("use_checkm2", [True, False])
def test_file_creation(use_checkm2: bool, test_data_dir: Path, tmp_path: Path) -> None:
    """Check files are created."""
    test_dir = tmp_path / "test_directory"
    expected_dir_name = "file_creation_checkm2" if use_checkm2 else "file_creation"
    expected_dir = test_data_dir / expected_dir_name
    genome_paths_file = test_data_dir / "genome_paths_file" / "valid.json"

    # Initialize the creators for each test case
    feature_protein_creator = MultiGenomeDataFileCreator(
        str(genome_paths_file),
        str(test_dir),
        "tests/data/example_files/checkm2/quality_report.tsv" if use_checkm2 else None,
        "tests/data/example_files/stats.json",
    )
    feature_protein_creator.create_all_tables()

    # Define file paths and expected line counts
    files_and_expected_lines = {
        test_dir / "contig.tsv": 89,
        test_dir / "contigset.tsv": 3,
        test_dir / "feature.tsv": 12028,
        test_dir / "feature_association.tsv": 12028,
        test_dir / "structural_annotation.tsv": 3,
    }

    err_list = []
    for file, n_lines in files_and_expected_lines.items():
        assert file.exists()
        parsed_lines = file.read_text().strip().split("\n")
        # parse file, check number of lines
        assert len(parsed_lines) == n_lines
        # read in the expected file and check they are identical
        expected_file = expected_dir / file.name
        expected_lines = expected_file.read_text().strip().split("\n")
        if parsed_lines == expected_lines:
            assert parsed_lines == expected_lines
        else:
            err_list.append(file)
            print(f"Error: {expected_file!s} differs from {file!s}")

    if err_list:
        assert err_list == []
