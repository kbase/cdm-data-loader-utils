import os
import unittest
import shutil
import csv

from  genome_loader_scripts.genome_loader import MultiGenomeDataFileCreator



class TestGenomeDataFileCreation(unittest.TestCase):
    def setUp(self):
        # Define paths for the input and output files
        self.genome_paths_file = "tests/data/genome_paths.json"
        self.output_dir = "tests/output"

        # Remove the output directory if it exists
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)




    def test_file_creation(self):
        # check file creation
        # Initialize the creators for each test case
        self.feature_protein_creator = MultiGenomeDataFileCreator(
            self.genome_paths_file, self.output_dir, None
        )
        self.feature_protein_creator.create_all_tables()
        # Define file paths and expected line counts
        self.files_and_expected_lines = {
            "tests/output/contig.tsv": 89,
            "tests/output/contigset.tsv": 3,
            "tests/output/feature.tsv": 12028,
            "tests/output/feature_association.tsv": 12028,
            "tests/output/structural_annotation.tsv": 3,
        }

        print ("\nTest: Check file creation")
        for file in self.files_and_expected_lines:
            expected_lines = self.files_and_expected_lines[file]
            print (f"checking if number of file lines in {file} is equal to  {expected_lines}")
            self.assertTrue(os.path.exists(file), f"{file} was not created.")

    @unittest.skip("Skipping checkm2  test")
    def test_checkm2(self):
        # check file creation
        # Initialize the creators for each test case

        print ("\nTest: Includes Checkm2 run")
        self.feature_protein_creator = MultiGenomeDataFileCreator(
            self.genome_paths_file, self.output_dir, 1
        )
        self.feature_protein_creator.create_all_tables()

        self.expected_scores = {
            "423e40e5b4056069f9b0bfb71a3c682b41a8b68200b617a6e19902fa5dac7e94": {"contamination": 1.09, "completeness": 99.99},
            "b32625f62ae333d2290a989cef9b3db75f462aeee8543aa11190af7f41c1d931": {"contamination": 2.27, "completeness": 99.97},
        }
        self.contigset_file = os.path.join(self.output_dir, "contigset.tsv")

        with open(self.contigset_file, 'r') as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter='\t')
            
            for row in reader:
                contigset_hash = row['contigset_hash']
                
                contamination = float(row['checkm2_contamination'])
                completeness = float(row['checkm2_completeness'])

                print (f"contamination for {contigset_hash} is {contamination}")
                print (f"completenessfor {contigset_hash} is {completeness}")

                if contigset_hash in self.expected_scores:
                    expected_contamination = self.expected_scores[contigset_hash]['contamination']
                    expected_completeness = self.expected_scores[contigset_hash]['completeness']

                    # Assert contamination and completeness match the expected values
                    self.assertAlmostEqual(contamination, expected_contamination, places=2, 
                                           msg=f"Contamination score mismatch for contigset hash {contigset_hash}")
                    self.assertAlmostEqual(completeness, expected_completeness, places=2, 
                                           msg=f"Completeness score mismatch for contigset hash {contigset_hash}")
                else:
                    self.fail(f"Unexpected contig_hash found in TSV file: {contig_hash}")



if __name__ == "__main__":
    unittest.main()

