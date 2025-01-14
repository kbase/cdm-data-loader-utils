"""Tests for the MultiGenomeDataFileCreator."""

import csv
import unittest
from pathlib import Path

import pytest

from genome_loader_scripts.genome_loader import MultiGenomeDataFileCreator


def test_file_creation(genome_paths_file: str, tmp_path: Path) -> None:
    """Check files are created."""

    test_dir = Path("tests") / "test_file_creation"
    # Initialize the creators for each test case
    feature_protein_creator = MultiGenomeDataFileCreator(genome_paths_file, test_dir, None)
    feature_protein_creator.create_all_tables()

    # Define file paths and expected line counts
    files_and_expected_lines = {
        test_dir / "contig.tsv": 89,
        test_dir / "contigset.tsv": 3,
        test_dir / "feature.tsv": 12028,
        test_dir / "feature_association.tsv": 12028,
        test_dir / "structural_annotation.tsv": 3,
    }

    for file, n_lines in files_and_expected_lines.items():
        assert file.exists()
        parsed_file = file.read_text().split("\n")
        # parse file, check number of lines
        assert len(parsed_file) == n_lines


# @pytest.skip("Skipping checkm2 test")
def test_checkm2(genome_paths_file: str, tmp_path: Path) -> None:
    # check file creation
    # Initialize the creators for each test case

    test_dir = Path("tests") / "test_checkm2"

    print("\nTest: Includes Checkm2 run")
    feature_protein_creator = MultiGenomeDataFileCreator(genome_paths_file, test_dir, 1)
    feature_protein_creator.create_all_tables()

    expected_scores = {
        "423e40e5b4056069f9b0bfb71a3c682b41a8b68200b617a6e19902fa5dac7e94": {
            "contamination": 1.09,
            "completeness": 99.99,
        },
        "b32625f62ae333d2290a989cef9b3db75f462aeee8543aa11190af7f41c1d931": {
            "contamination": 2.27,
            "completeness": 99.97,
        },
    }
    contigset_file = test_dir / "contigset.tsv"

    with contigset_file.open() as tsvfile:
        reader = csv.DictReader(tsvfile, delimiter="\t")

        for row in reader:
            contigset_hash = row["contigset_hash"]
            if contigset_hash in expected_scores:
                for c in ["contamination", "completeness"]:
                    c_value = float(row[f"checkm2_{c}"])
                    # Assert contamination and completeness match the expected values
                    assert c_value == pytest.approx(expected_scores[contigset_hash][c], rel=1e-3)
            else:
                pytest.fail(f"Unexpected contigset_hash found in TSV file: {contigset_hash}")


if __name__ == "__main__":
    unittest.main()
