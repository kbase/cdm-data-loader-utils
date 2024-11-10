import os
import unittest
import shutil

from  genome_loader_scripts.genome_loader import MultiGenomeDataFileCreator



class TestGenomeDataFileCreation(unittest.TestCase):
    def setUp(self):
        # Define paths for the input and output files
        self.genome_paths_file = "tests/data/genome_paths.json"
        self.output_dir = "tests/output"
        # Remove the output directory if it exists
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)


        # Define file paths and expected line counts
        self.files_and_expected_lines = {
            "tests/output/contig.tsv": 89,
            "tests/output/contigset.tsv": 3,
            "tests/output/feature.tsv": 12028,
            "tests/output/feature_association.tsv": 12028,
            "tests/output/structural_annotation.tsv": 3,
        }



    def test_file_creation(self):
        # check file creation
        # Initialize the creators for each test case
        self.feature_protein_creator = MultiGenomeDataFileCreator(
            self.genome_paths_file, self.output_dir, None
        )
        self.feature_protein_creator.create_all_tables()

        for file in self.files_and_expected_lines:
            expected_lines = self.files_and_expected_lines[file]
            print (f"checking if number of file lines in {file} is equal to  {expected_lines}")
            self.assertTrue(os.path.exists(file), f"{file} was not created.")


if __name__ == "__main__":
    unittest.main()

