import os
import unittest
import shutil
import csv

from genome_loader_scripts.genome_loader import MultiGenomeDataFileCreator

# Define paths for the input and output files
genome_paths_file = "tests/data/genome_paths.json"


class TestGenomeDataFileCreation(unittest.TestCase):
    def test_file_creation(self):
        # check file creation
        # Initialize the creators for each test case
        feature_protein_creator = MultiGenomeDataFileCreator(
            genome_paths_file, "tests/test_file_creation", None
        )
        feature_protein_creator.create_all_tables()
        # Define file paths and expected line counts
        files_and_expected_lines = {
            "tests/test_file_creation/contig.tsv": 89,
            "tests/test_file_creation/contigset.tsv": 3,
            "tests/test_file_creation/feature.tsv": 12028,
            "tests/test_file_creation/feature_association.tsv": 12028,
            "tests/test_file_creation/structural_annotation.tsv": 3,
        }

        print("\nTest: Check file creation")
        for file in files_and_expected_lines:
            expected_lines = files_and_expected_lines[file]
            print(f"checking if number of file lines in {file} is equal to  {expected_lines}")
            self.assertTrue(os.path.exists(file), f"{file} was not created.")

    # @unittest.skip("Skipping checkm2 test")
    def test_checkm2(self):
        # check file creation
        # Initialize the creators for each test case

        print("\nTest: Includes Checkm2 run")
        feature_protein_creator = MultiGenomeDataFileCreator(
            genome_paths_file, "tests/test_checkm2", 1
        )
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
        contigset_file = os.path.join("tests/test_checkm2", "contigset.tsv")

        with open(contigset_file) as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")

            for row in reader:
                contigset_hash = row["contigset_hash"]

                contamination = float(row["checkm2_contamination"])
                completeness = float(row["checkm2_completeness"])

                print(f"contamination for {contigset_hash} is {contamination}")
                print(f"completenessfor {contigset_hash} is {completeness}")

                if contigset_hash in expected_scores:
                    expected_contamination = expected_scores[contigset_hash]["contamination"]
                    expected_completeness = expected_scores[contigset_hash]["completeness"]

                    # Assert contamination and completeness match the expected values
                    self.assertAlmostEqual(
                        contamination,
                        expected_contamination,
                        places=2,
                        msg=f"Contamination score mismatch for contigset hash {contigset_hash}",
                    )
                    self.assertAlmostEqual(
                        completeness,
                        expected_completeness,
                        places=2,
                        msg=f"Completeness score mismatch for contigset hash {contigset_hash}",
                    )
                else:
                    self.fail(f"Unexpected contig_hash found in TSV file: {contigset_hash}")


if __name__ == "__main__":
    unittest.main()
