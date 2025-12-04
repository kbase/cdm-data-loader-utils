"""Loader for genome files."""

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path
from typing import Any

from Bio import SeqIO

from cdm_data_loader_utils.parsers.bbmap_stats import get_bbmap_stats
from cdm_data_loader_utils.parsers.checkm2 import get_checkm2_data
from cdm_data_loader_utils.parsers.genome_paths import get_genome_paths
from cdm_data_loader_utils.utils import calculate_hash as ch

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
    def __init__(self, contigset_file, gff_file, protein_file, output_dir) -> None:
        self.contigset_file = contigset_file
        self.gff_file = gff_file
        self.protein_file = protein_file
        self.contigset_dir = Path(contigset_file).parent
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

        self.contigset_hash, self.contig_id_map, self.protein_id_map = ch.compute_hash(
            self.contigset_file, self.protein_file
        )

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

    # TODO: Update to use prepare_feature_data
    def prepare_gff3_data(self):
        """Prepare data for insertion into the database."""
        print(f"Preparing GFF3 data from: {self.gff_file}")
        if not self.gff_hash:
            print("Error: GFF file hash not calculated.")
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

                    feature_hash = ch.generate_hash_id(contigset_hash, contig_hash, start, end, feature_type)

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
            self.contigset_hash, self.contig_id_map, self.protein_id_map = ch.compute_hash(
                self.contigset_file, self.protein_file
            )

        open_func = gzip.open if self.contigset_file.endswith(".gz") else open

        with open_func(self.contigset_file, "rt") as handle:
            sequences = SeqIO.parse(handle, "fasta")
            for seq_record in sequences:
                contig_name = seq_record.id
                sequence = str(seq_record.seq).upper()
                length = len(sequence)
                gc_content = (sequence.count("G") + sequence.count("C")) / length if length > 0 else 0
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
            self.contigset_hash, self.contig_id_map, self.protein_id_map = ch.compute_hash(
                self.contigset_file, self.protein_file
            )

        print(self.contigset_file)

        if self.contigset_hash:
            self.contigset = {
                "contigset_hash": self.contigset_hash,
            }

    def prepare_structural_annotation(self):
        self.structural_annotation = {
            "contigset_hash": self.contigset_hash,
            "gff_hash": self.gff_hash,
            "contigset_file": self.contigset_file,
            "gff_file": self.gff_file,
        }

    def run(self) -> tuple[Any]:
        """Import the genomic, feature, and protein data from a set of files."""
        self.gff_hash = ch.calculate_sha256_checksums(self.gff_file)
        self.prepare_gff3_data()
        self.prepare_contig_data()
        self.prepare_contigset_data()
        self.prepare_structural_annotation()

        return (
            self.contigset,
            self.contigs,
            self.features,
            self.feature_associations,
            self.structural_annotation,
        )


class MultiGenomeDataFileCreator:
    """Parser that takes in GFF, FAA, and FNA files and generates CDM table data."""

    def __init__(
        self: "MultiGenomeDataFileCreator",
        genome_paths_file: str,
        output_dir: str,
        checkm2_file: str | None = None,
        stats_file: str | None = None,
    ) -> None:
        """Initialise the MultiGenomeDataFileCreator."""
        errs = []
        if not genome_paths_file or not genome_paths_file.strip():
            errs.append("Missing genome_paths_file")
        else:
            self.genome_paths_file = Path(genome_paths_file.strip())
            if not self.genome_paths_file.exists():
                errs.append(f"genome_paths_file '{self.genome_paths_file}' does not exist")

        if not output_dir or not output_dir.strip():
            errs.append("Missing output_dir")
        else:
            self.output_dir = Path(output_dir.strip())

        if errs:
            err_msg = f"MultiGenomeDataFileCreator init error:\n{'\n'.join(errs)}"
            raise RuntimeError(err_msg)

        self.checkm2_file = None
        self.stats_file = None
        if checkm2_file and checkm2_file.strip():
            self.checkm2_file = Path(checkm2_file.strip())
        if stats_file and stats_file.strip():
            self.stats_file = Path(stats_file.strip())

    def write_to_tsv(
        self: "MultiGenomeDataFileCreator",
        contigset,
        contigs,
        features,
        associations,
        structural_annotation,
        headers_written,
    ):
        """Write data to TSV files incrementally."""
        try:
            # Contigset
            contigset_out_file = self.output_dir / "contigset.tsv"
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
            ]

            with contigset_out_file.open("a", newline="") as f_out:
                writer = csv.DictWriter(f_out, fieldnames=contigset_fields, delimiter="\t")
                if write_header:
                    writer.writeheader()
                    headers_written["contigset"] = True
                writer.writerow(contigset)
            print(f"Contigset appended to {contigset_out_file}")

            # Contig
            contig_file = self.output_dir / "contig.tsv"
            write_header = not headers_written["contigs"]
            with contig_file.open("a", newline="") as f_out:
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
            sa_file = self.output_dir / "structural_annotation.tsv"
            write_header = not headers_written["structural_annotation"]
            with sa_file.open("a", newline="") as f_out:
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
            features_file = self.output_dir / "feature.tsv"
            write_header = not headers_written["features"]
            with features_file.open("a", newline="") as f_out:
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
            associations_file = self.output_dir / "feature_association.tsv"  # Renamed file
            write_header = not headers_written["associations"]
            with associations_file.open("a", newline="") as f_out:
                fieldnames = ["feature_hash", "gff_hash", "feature_attributes"]
                writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
                if write_header:
                    writer.writeheader()
                    headers_written["associations"] = True
                writer.writerows(associations)
            print(f"Feature associations appended to {associations_file}")

        except Exception as e:
            print(f"Error writing to TSV files: {e}")

    def create_all_tables(self: "MultiGenomeDataFileCreator"):
        checkm2_data = None
        stats_data = None

        genome_paths = get_genome_paths(self.genome_paths_file)
        # check for checkm2 and stats data
        if self.checkm2_file:
            checkm2_data = get_checkm2_data(self.checkm2_file)
        if self.stats_file:
            stats_data = get_bbmap_stats(self.stats_file)

        if self.output_dir.exists():
            print(f"Error: The directory {self.output_dir} already exists.")
            print(f"Error: Remove directory {self.output_dir} to continue")
            sys.exit(1)  # Exit with a non-zero code to indicate an error

        # create the output directory
        self.output_dir.mkdir(parents=True)

        genome_ids = sorted(genome_paths.keys())
        headers_written = {
            "contigs": False,
            "contigset": False,
            "structural_annotation": False,
            "features": False,
            "associations": False,
        }

        for gid in genome_ids:
            paths = genome_paths.get(gid)
            # this should not happen
            if not paths:
                print(f"No paths found for genome ID: {gid}")
                continue

            # _scaffolds.fna
            contigset_file = paths.get("fna")
            # _genes.gff
            gff_file = paths.get("gff")
            # _genes.faa
            protein_file = paths.get("protein")
            if not contigset_file or not gff_file or not protein_file:
                print(f"Missing file paths for genome ID: {gid}")
                continue

            print(f"\n==Processing contigset {contigset_file}==\n")
            parser = GenomeDataFileCreator(contigset_file, gff_file, protein_file, self.output_dir)
            contigset, contigs, features, associations, structural_annotation = parser.run()

            # get the stem of the contigset file name and check the checkm2 data for it
            if checkm2_data and Path(contigset_file).stem in checkm2_data:
                contigset.update(checkm2_data[Path(contigset_file).stem])
            # check for the contigset file name
            if stats_data and Path(contigset_file).name in stats_data:
                contigset.update(stats_data[Path(contigset_file).name])
                del contigset["filename"]

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

    @classmethod
    def from_args(cls: type["MultiGenomeDataFileCreator"]) -> "MultiGenomeDataFileCreator":
        """Parse command-line arguments and create an instance of MultiGenomeDataFileCreator."""
        parser = argparse.ArgumentParser(description="Create tables for a specified genome ID.")
        parser.add_argument("--genome_paths_file", type=str, required=True, help="Path to genome paths JSON file")

        parser.add_argument("--output_dir", type=str, required=True, help="Output path to save the table files")

        parser.add_argument(
            "--stats",
            type=str,
            default="",
            help="Path to the directory containing the `stats.json` fna stats file",
        )

        parser.add_argument(
            "--checkm2",
            type=str,
            default="",
            help="Path to the directory containing the `quality_report.tsv` file produced by checkm2",
        )

        args = parser.parse_args()

        return cls(
            genome_paths_file=args.genome_paths_file.strip(),
            stats_file=args.stats_file.strip(),
            checkm2_file=args.checkm2_file.strip(),
            output_dir=args.output_dir.strip(),
        )

    def run(self: "MultiGenomeDataFileCreator") -> None:
        """Main method to parse input files and generate TSV output."""
        self.create_all_tables()


if __name__ == "__main__":
    # Create an instance of MultiGenomeDataFileCreator using command-line arguments
    creator = MultiGenomeDataFileCreator.from_args()
    creator.run()
