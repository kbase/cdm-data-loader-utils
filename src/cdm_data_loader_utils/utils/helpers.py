"""Misc database-related helpers."""


def _ensembl_type(xref: str) -> str:
    """Given an Ensembl ID, return the likely type."""
    prefixes = {
        "ENST": "transcript",
        "ENSP": "protein sequence",
        "ENSG": "gene sequence",
        "ENSE": "exon",
        "ENSFM": "protein family",
        "ENSGT": "gene tree",
        "ENSR": "regulatory feature",
    }

    return prefixes.get(xref[0:4], prefixes.get(xref[0:5], "sequence"))


def _refseq_type(xref: str) -> str | None:
    prefixes = {
        "AC_": "Genomic: complete genomic molecule, usually alternate assembly",
        "NC_": "Genomic: complete genomic molecule, usually reference assembly",
        "NG_": "Genomic: incomplete genomic region",
        "NT_": "Genomic: contig or scaffold, clone-based or Whole Genome Shotgun sequence data",
        "NW_": "Genomic: contig or scaffold, primarily Whole Genome Shotgun sequence data",
        "NZ_": "Genomic: complete genomes and unfinished WGS data; an ordered collection of WGS sequence for a genome",
        "NM_": "mRNA: protein-coding transcripts (usually curated)",
        "NR_": "RNA: non-protein-coding transcripts",
        "XM_": "mRNA: predicted model protein-coding transcript (computed)",
        "XR_": "RNA: predicted model non-protein-coding transcript (computed)",
        "AP_": "Protein: annotated on AC_ alternate assembly",
        "NP_": "Protein: associated with an NM_ or NC_ accession",
        "YP_": "Protein: annotated on genomic molecules without an instantiated transcript record (computed)",
        "XP_": "Protein: predicted model, associated with an XM_ accession (computed)",
        "WP_": "Protein: non-redundant across multiple strains and species",
    }
    first_three_letters = xref[0:3]
    return prefixes.get(first_three_letters)
