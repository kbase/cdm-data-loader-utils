"""Utils for calculating hashes."""

import gzip
import hashlib
from typing import Any

from src.parsers.fasta import read_fasta


def _hash_string(s: str) -> str:
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


class HashSeq(str):
    def __new__(cls, v) -> "HashSeq":
        return super().__new__(cls, v.upper())

    @property
    def hash_value(self):
        return _hash_string(self)


class HashSeqList(list):
    def append(self: "HashSeqList", o: str | HashSeq) -> None:
        if isinstance(o, str):
            super().append(HashSeq(o))
        elif isinstance(o, HashSeq):
            super().append(o)
        else:
            err_msg = f"Invalid type: {type(o)}"
            raise TypeError(err_msg)

    @property
    def hash_value(self: "HashSeqList") -> str:
        h_list = [x.hash_value for x in self]
        hash_seq = "_".join(sorted(h_list))
        return _hash_string(hash_seq)


def contig_set_hash(features) -> str:
    hl = HashSeqList()
    for contig in features:
        seq = HashSeq(contig.seq)
        hl.append(seq)
    return hl.hash_value


# TODO: check protein_id_map issue
def compute_hash(
    contigset_file: str, protein_file: str
) -> tuple[str, dict[str, Any], dict[str, Any]]:
    """
    Compute the hash of the entire contigset and contigs.
    """
    contig_hash_dict = {}
    contigs = read_fasta(contigset_file)
    contigset_hash = contig_set_hash(contigs)
    for f in contigs:
        contig_hash_dict[f.id] = HashSeq(f.seq).hash_value
    contig_id_map = contig_hash_dict

    protein_hash_dict = {}
    proteins = read_fasta(protein_file)
    for f in proteins:
        protein_hash_dict[f.id] = HashSeq(f.seq).hash_value
    protein_id_map = protein_hash_dict

    return contigset_hash, contig_id_map, protein_id_map


def generate_file_sha256(filepath: str, blocksize=65536) -> str | None:
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


def generate_hash_id(*args):
    """Generate a hash-based ID from the given arguments."""
    unique_string = "".join(map(str, args))
    return hashlib.sha256(unique_string.encode("utf-8")).hexdigest()


def calculate_sha256_checksums(gff_file: str) -> str | None:
    """Calculate checksums for the contigset and its contigs."""

    print(f"Calculating sha256 for GFF: {gff_file}")
    gff_hash = generate_file_sha256(gff_file)
    if not gff_hash:
        print(f"Error calculating sha256 for GFF file {gff_file}")
        return None
    return gff_hash
