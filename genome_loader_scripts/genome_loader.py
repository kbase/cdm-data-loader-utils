"""Loader for genome files."""

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
import uuid

from Bio import SeqIO

from . import calculate_hash as ch

NUM_THREADS_CHECKM2_RUN = 10


# Define SO terms mapping
so_terms = {
    "gene": "SO:0000704",
    "pseudogene": "SO:0000336",
    "ncRNA_gene": "SO:0001263",
    "mRNA": "SO:0000234",
    "CDS": "SO:0000316",
    "exon": "SO:0000147",
    "five_prime_UTR": "SO:0000204",
    "three_prime_UTR": "SO:0000205",
    "ncRNA": "SO:0000655",
    "rRNA": "SO:0000252",
    "tRNA": "SO:0000253",
    "SRP_RNA": "SO:0000590",
    "RNase_P_RNA": "SO:0000386",
    "riboswitch": "SO:0000035",
    "direct_repeat": "SO:0000319",
    "origin_of_replication": "SO:0000296",
    "CRISPR": "SO:0001459",
    "mobile_genetic_element": "SO:0001037",
    "region": "SO:0000001",
    "sequence_feature": "SO:0000110",
}


class GenomeDataFileCreator:
    def __init__(
        self, contigset_file, gff_file, protein_file, output_dir, run_checkm2_option
    ) -> None:
        self.contigset_file = contigset_file
        self.gff_file = gff_file
        self.protein_file = protein_file
        self.run_checkm2_option = run_checkm2_option
        self.output_dir = output_dir

        self.features = []
        self.feature_associations = []
        self.feature_protein_associations = []
        self.gff_hash = None
        self.protein_count = 0
        self.feature_with_protein_mapping = 0

        # TODO
        # Put self.genome_id and use in contigset

        self.contigset = {}
        self.contigs = []
        self.structural_annotation = {}

        self.contigset_hash, self.contig_id_map, self.protein_id_map = self.compute_hash(
            self.contigset_file, self.protein_file
        )

    # TODO: check protein_id_map issue
    @staticmethod
    def compute_hash(contigset_file, protein_file):
        """
        Compute the hash of the entire contigset and contigs
        """
        contig_hash_dict = dict()
        contigs = ch.read_fasta2(contigset_file)
        contigset_hash = ch.contig_set_hash(contigs)
        for f in contigs:
            contig_hash_dict[f.id] = ch.HashSeq(f.seq).hash_value
        contig_id_map = contig_hash_dict

        protein_hash_dict = dict()
        proteins = ch.read_fasta2(protein_file)
        for f in proteins:
            protein_hash_dict[f.id] = ch.HashSeq(f.seq).hash_value
        protein_id_map = protein_hash_dict

        return contigset_hash, contig_id_map, protein_id_map

    @staticmethod
    def generate_file_sha256(filepath, blocksize=65536):
        """Generate the SHA-256 checksum of a file's decompressed content."""
        sha256 = hashlib.sha256()
        open_func = gzip.open if filepath.endswith(".gz") else open
        try:
            with open_func(filepath, "rt", encoding="utf-8", errors="ignore") as f:
                for block in iter(lambda: f.read(blocksize), ""):
                    sha256.update(block.encode("utf-8"))
            return sha256.hexdigest()
        except Exception as e:
            print(f"Error generating SHA-256 for {filepath}: {e}")
            return None

    @staticmethod
    def generate_hash_id(*args):
        """Generate a hash-based ID from the given arguments."""
        unique_string = "".join(map(str, args))
        return hashlib.sha256(unique_string.encode("utf-8")).hexdigest()

    @staticmethod
    def parse_attributes(attributes_str):
        """Parse the GFF3 attributes field into a dictionary."""
        attributes = {}
        for attribute in attributes_str.strip(";").split(";"):
            if "=" in attribute:
                key, value = attribute.split("=", 1)
                key = key.strip('"')
                value = value.strip('"')
                attributes[key] = value
        return attributes

    def run_bbmap(self, fna_path):
        """
        Runs BBMap on a given genome and returns the mapping output.
        """
        bbmap_command = [
            "stats.sh",
            f"in={fna_path}",  # Input genome file (FASTA format)
            "--format=8",
        ]

        try:
            result = subprocess.run(
                bbmap_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True
            )
            json_output = result.stdout.decode("utf-8")
            parsed_json = json.loads(json_output)
            return parsed_json
        except subprocess.CalledProcessError as e:
            print(f"BBMap failed with error: {e.stderr.decode('utf-8')}")
            return None

    def run_checkm2(self, contigset_path, output_dir):
        # Run CheckM2 'predict' command on the specified contigset file
        print("Running checkm2\n")
        try:
            result = subprocess.run(
                [
                    "checkm2",
                    "predict",
                    "--threads",
                    str(NUM_THREADS_CHECKM2_RUN),
                    "--input",
                    contigset_path,
                    "--output_directory",
                    output_dir,
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            # Open and read the TSV file
            tsv_file_path = os.path.join(output_dir, "quality_report.tsv")
            with open(tsv_file_path, "r") as file:
                reader = csv.DictReader(file, delimiter="\t")
                # Ensure required columns are present
                if (
                    "Completeness" not in reader.fieldnames
                    or "Contamination" not in reader.fieldnames
                ):
                    err_msg = "TSV file does not contain 'Completeness' or 'Contamination' columns."
                    raise ValueError(err_msg)

                # Loop through each row and extract the metrics
                for row in reader:
                    # Convert completeness and contamination to floats
                    completeness = float(row["Completeness"])
                    contamination = float(row["Contamination"])
                    return {
                        "checkm2_contamination": contamination,
                        "checkm2_completeness": completeness,
                    }

        except subprocess.CalledProcessError as e:
            print("CheckM2 failed with the following error:")
            print(e.stderr)
            return None

    # TODO: make it static method
    def calculate_sha256_checksums(self):
        """Calculate  checksums for the contigset and its contigs."""

        print(f"Calculating sha256 for GFF: {self.gff_file}")
        self.gff_hash = self.generate_file_sha256(self.gff_file)
        if not self.gff_hash:
            print(f"Error calculating sha256 for GFF file {self.gff_file}")
            return

    # TODO: Update to use prepare_feature_data
    def prepare_gff3_data(self):
        """Prepare data for insertion into the database."""
        print(f"Preparing GFF3 data from: {self.gff_file}")
        if not self.gff_hash:
            print(f"Error: GFF file hash not calculated.")
            return

        open_func = gzip.open if self.gff_file.endswith(".gz") else open

        try:
            with open_func(self.gff_file, "rt", encoding="utf-8", errors="ignore") as file:
                reader = csv.reader(file, delimiter="\t")
                for row in reader:
                    if row[0].startswith("#") or len(row) < 9:
                        continue

                    seq_id = row[0]
                    source = row[1]
                    feature_type = row[2]
                    start = int(row[3])
                    end = int(row[4])
                    score = row[5] if row[5] != "." else None
                    strand = row[6] if row[6] in ["+", "-"] else None
                    phase = row[7] if row[7] in ["0", "1", "2"] else None
                    attributes_str = row[8]

                    feature_ontology = so_terms.get(feature_type, "")

                    # Parse attributes
                    attributes = self.parse_attributes(attributes_str)
                    feature_name = attributes.get("ID", None)
                    protein_accession = attributes.get("protein_id", None)
                    protein_hash = None

                    # Mapping between features in gff and proteins in protein file
                    if feature_type == "CDS":
                        protein_hash = self.protein_id_map.get(feature_name, None)

                    # Note: This is to handle cases where protein_id in protein file are stashed in
                    # protein_id field of attributes
                    # TODO: check if we need more ways to handle this
                    if not protein_hash and "protein_id" in attributes:
                        protein_hash = self.protein_id_map.get("protein_id", None)

                    if protein_hash:
                        self.feature_with_protein_mapping += 1

                    contig_hash = self.contig_id_map.get(seq_id, "")
                    contigset_hash = self.contigset_hash

                    # Generate a unique hash ID for each feature
                    # TODO: May be switch this to use feature_ontology
                    # If you are working with contigsets of very close strains
                    # the same contig_hash may appear in each.
                    # Including contigset_hash ensures that features
                    #  are uniquely identified across different contigsets.

                    feature_hash = self.generate_hash_id(
                        contigset_hash, contig_hash, start, end, feature_type
                    )

                    # Prepare feature data including hash of the contigset and contig
                    feature_data = {
                        "contigset_hash": contigset_hash,
                        "feature_hash": feature_hash,
                        "feature_type": feature_type,
                        "feature_ontology": feature_ontology,
                        "start": start,
                        "end": end,
                        "strand": strand,
                        "phase": phase,
                        "contig_hash": self.contig_id_map.get(seq_id, ""),
                        "protein_hash": protein_hash,
                    }
                    self.features.append(feature_data)

                    # Add all other attributes to the associations

                    self.feature_associations.append(
                        {
                            "gff_hash": self.gff_hash,
                            "feature_hash": feature_hash,
                            "feature_attributes": json.dumps(attributes).replace('"', ""),
                        }
                    )

            print(f"Finished preparing GFF3 data. Total features: {len(self.features)}")
        except Exception as e:
            print(f"Error reading GFF file {self.gff_file}: {e}")

    def prepare_contig_data(self):
        """
        Calculate statistics for each contig in the assembly file.
        Handles both compressed (.gz) and uncompressed files.
        """
        if not self.contigset_hash:
            self.contigset_hash, self.contig_id_map, self.protein_id_map = self.compute_hash(
                self.contigset_file, self.protein_file
            )

        open_func = gzip.open if self.contigset_file.endswith(".gz") else open

        with open_func(self.contigset_file, "rt") as handle:
            sequences = SeqIO.parse(handle, "fasta")
            for seq_record in sequences:
                contig_name = seq_record.id
                sequence = str(seq_record.seq).upper()
                length = len(sequence)
                gc_content = (
                    (sequence.count("G") + sequence.count("C")) / length if length > 0 else 0
                )
                # TODO: check if contig_id_map generated or not
                self.contigs.append(
                    {
                        "contig_hash": self.contig_id_map[contig_name],
                        "contig_name": contig_name,
                        "length": length,
                        "gc_content": gc_content,
                        "contigset_hash": self.contigset_hash,
                        "contigset_file": self.contigset_file,
                    }
                )

    def prepare_contigset_data(self):
        """
        Calculate statistics for each contig in the assembly file.
        Handles both compressed (.gz) and uncompressed files.
        """
        if not self.contigset_hash:
            self.contigset_hash, self.contig_id_map, self.protein_id_map = self.compute_hash(
                self.contigset_file, self.protein_file
            )

        print(self.contigset_file)

        if self.contigset_hash:
            # Get bbmap data if available
            bbmap_info = self.run_bbmap(self.contigset_file)

            checkm2_info = {"checkm2_contamination": "", "checkm2_completeness": ""}

            if self.run_checkm2_option is not None:
                unique_id = str(uuid.uuid4())
                temp_dir = os.path.join(self.output_dir, unique_id)
                checkm2_info = self.run_checkm2(self.contigset_file, temp_dir)

                # Note: You can keep this directory for checkm debugging
                if os.path.exists(temp_dir):
                    shutil.rmtree(temp_dir)

            # Create a QC entry hash bbmap data
            self.contigset = {
                "contigset_hash": self.contigset_hash,
                **(bbmap_info or {}),
                **(checkm2_info or {}),
            }

    def prepare_structural_annotation(self):
        self.structural_annotation = {
            "contigset_hash": self.contigset_hash,
            "gff_hash": self.gff_hash,
            "contigset_file": self.contigset_file,
            "gff_file": self.gff_file,
        }


class MultiGenomeDataFileCreator:
    def __init__(self, genome_paths_file, output_dir, run_checkm2_option):
        self.genome_paths_file = genome_paths_file
        self.output_dir = output_dir
        self.run_checkm2_option = run_checkm2_option

    def process_files(self, contigset_file, gff_file, protein_file):
        print(f"\n==Processing contigset {contigset_file}==\n")
        parser = GenomeDataFileCreator(
            contigset_file, gff_file, protein_file, self.output_dir, self.run_checkm2_option
        )
        parser.calculate_sha256_checksums()
        parser.prepare_gff3_data()
        parser.prepare_contig_data()
        parser.prepare_contigset_data()
        parser.prepare_structural_annotation()

        return (
            parser.contigset,
            parser.contigs,
            parser.features,
            parser.feature_associations,
            parser.structural_annotation,
        )

    def write_to_tsv(
        self, contigset, contigs, features, associations, structural_annotation, headers_written
    ):
        """Write data to TSV files incrementally."""
        try:
            # Contigset
            contigset_out_file = os.path.join(self.output_dir, "contigset.tsv")
            write_header = not headers_written["contigset"]

            contigset_fields = [
                "contigset_hash",
                "checkm2_contamination",
                "checkm2_completeness",
                "scaffolds",
                "contigs",
                "scaf_bp",
                "contig_bp",
                "gap_pct",
                "scaf_N50",
                "scaf_L50",
                "ctg_N50",
                "ctg_L50",
                "scaf_N90",
                "scaf_L90",
                "ctg_N90",
                "ctg_L90",
                "scaf_logsum",
                "scaf_powsum",
                "ctg_logsum",
                "ctg_powsum",
                "asm_score",
                "scaf_max",
                "ctg_max",
                "scaf_n_gt50K",
                "scaf_l_gt50k",
                "scaf_pct_gt50K",
                "gc_avg",
                "gc_std",
                "filename",
            ]

            with open(contigset_out_file, "a", newline="") as f_out:
                writer = csv.DictWriter(f_out, fieldnames=contigset_fields, delimiter="\t")
                if write_header:
                    writer.writeheader()
                    headers_written["contigset"] = True
                writer.writerow(contigset)
            print(f"Contigset appended to {contigset_out_file}")

            # Contig
            contig_file = os.path.join(self.output_dir, "contig.tsv")
            write_header = not headers_written["contigs"]
            with open(contig_file, "a", newline="") as f_out:
                writer = csv.DictWriter(
                    f_out,
                    fieldnames=[
                        "contig_hash",
                        "contig_name",
                        "length",
                        "gc_content",
                        "contigset_hash",
                        "contigset_file",
                    ],
                    delimiter="\t",
                )
                if write_header:
                    writer.writeheader()
                    headers_written["contig"] = True
                writer.writerows(contigs)
            print(f"Contig appended to {contig_file}")

            # Structural annotation
            sa_file = os.path.join(self.output_dir, "structural_annotation.tsv")
            write_header = not headers_written["structural_annotation"]
            with open(sa_file, "a", newline="") as f_out:
                writer = csv.DictWriter(
                    f_out,
                    fieldnames=["contigset_hash", "gff_hash", "contigset_file", "gff_file"],
                    delimiter="\t",
                )
                if write_header:
                    writer.writeheader()
                    headers_written["structural_annotation"] = True
                writer.writerow(structural_annotation)
            print(f"Structural annotation appended to {sa_file}")

            # Features
            features_file = os.path.join(self.output_dir, "feature.tsv")
            write_header = not headers_written["features"]
            with open(features_file, "a", newline="") as f_out:
                contig_fieldnames = [
                    "contigset_hash",
                    "contig_hash",
                    "feature_hash",
                    "feature_type",
                    "feature_ontology",
                    "start",
                    "end",
                    "strand",
                    "phase",
                    "protein_hash",
                ]
                writer = csv.DictWriter(f_out, fieldnames=contig_fieldnames, delimiter="\t")
                if write_header:
                    writer.writeheader()
                    headers_written["features"] = True
                writer.writerows(features)
            print(f"Features appended to {features_file}")

            # Feature associations
            associations_file = os.path.join(
                self.output_dir, "feature_association.tsv"
            )  # Renamed file
            write_header = not headers_written["associations"]
            with open(associations_file, "a", newline="") as f_out:
                fieldnames = ["feature_hash", "gff_hash", "feature_attributes"]
                writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
                if write_header:
                    writer.writeheader()
                    headers_written["associations"] = True
                writer.writerows(associations)
            print(f"Feature associations appended to {associations_file}")

        except Exception as e:
            print(f"Error writing to TSV files: {e}")

    def create_all_tables(self):
        try:
            os.makedirs(self.output_dir, exist_ok=False)
        except FileExistsError:
            print(f"Error: The directory {self.output_dir} already exists.")
            print(f"Error: Remove directory {self.output_dir} to continue")
            sys.exit(1)  # Exit with a non-zero code to indicate an error
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            sys.exit(1)  # Exit with a non-zero code to indicate an error

        with open(self.genome_paths_file, "r") as json_file:
            genome_paths = json.load(json_file)

        genome_ids = list(genome_paths.keys())
        headers_written = {
            "contigs": False,
            "contigset": False,
            "structural_annotation": False,
            "features": False,
            "associations": False,
        }

        for gid in genome_ids:
            paths = genome_paths.get(gid)
            if paths:
                contigset_file = paths.get("fna")
                gff_file = paths.get("gff")
                protein_file = paths.get("protein")
                if contigset_file and gff_file and protein_file:
                    contigset, contigs, features, associations, structural_annotation = (
                        self.process_files(contigset_file, gff_file, protein_file)
                    )

                    # Write data to TSV files incrementally
                    self.write_to_tsv(
                        contigset,
                        contigs,
                        features,
                        associations,
                        structural_annotation,
                        headers_written,
                    )

                    # Clear data to free memory
                    contigset.clear()
                    contigs.clear()
                    features.clear()
                    associations.clear()
                else:
                    print(f"Missing file paths for genome ID: {gid}")
            else:
                print(f"No paths found for genome ID: {gid}")

    @classmethod
    def from_args(cls):
        """Parse command-line arguments and create an instance of MultiGenomeDataFileCreator."""
        parser = argparse.ArgumentParser(description="Create tables for a specified genome ID.")
        parser.add_argument(
            "--genome_paths_file", type=str, required=True, help="Path to genome paths JSON file"
        )
        parser.add_argument(
            "--output_dir", type=str, required=True, help="Output path to save the table files"
        )
        parser.add_argument(
            "--run_checkm2_option", type=str, default=None, help="Run checkm2, use 1 to run checkm2"
        )

        args = parser.parse_args()

        return cls(
            genome_paths_file=args.genome_paths_file,
            output_dir=args.output_dir,
            run_checkm2_option=args.run_checkm2_option,
        )

    def run(self):
        """Main method to execute the creation of all the  table."""
        self.create_all_tables()


if __name__ == "__main__":
    # Create an instance of MultiGenomeDataFileCreator using command-line arguments
    creator = MultiGenomeDataFileCreator.from_args()
    creator.run()