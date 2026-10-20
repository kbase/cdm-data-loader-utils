"""Shared fixtures for ncbi_ftp tests."""

from pathlib import PurePosixPath
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

import cdm_data_loaders.ncbi_ftp.manifest as manifest_mod
from tests.conftest import TEST_BUCKET as TEST_BUCKET_STR

AWS_REGION = "us-east-1"
TEST_BUCKET: PurePosixPath = PurePosixPath(TEST_BUCKET_STR)


# Minimal assembly_summary_refseq.txt content (tab-separated, 20+ columns)
SAMPLE_SUMMARY = (
    "# assembly_accession\tbioproject\tbiosample\twgs_master\trefseq_category\t"
    "taxid\tspecies_taxid\torganism_name\tinfraspecific_name\tisolate\t"
    "version_status\tassembly_level\trelease_type\tgenome_rep\tseq_rel_date\t"
    "asm_name\t16\t17\t18\tftp_path\n"
    "GCF_000001215.4\tPRJNA13812\tSAMN02803731\t\treference genome\t7227\t7227\t"
    "Drosophila melanogaster\t\t\tlatest\tChromosome\tMajor\tFull\t2014/10/21\t"
    "Release_6_plus_ISO1_MT\t\t\t\t"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT\n"
    "GCF_000001405.40\tPRJNA168\tna\t\treference genome\t9606\t9606\t"
    "Homo sapiens\t\t\tlatest\tChromosome\tPatch\tFull\t2022/02/03\t"
    "GRCh38.p14\t\t\t\t"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14\n"
    "GCF_000005845.2\tPRJNA57779\tSAMN02604091\t\trepresentative genome\t511145\t562\t"
    "Escherichia coli\t\t\treplaced\tComplete Genome\tMajor\tFull\t2013/09/26\t"
    "ASM584v2\t\t\t\t"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2\n"
    "GCF_000009999.1\tPRJNA999\tSAMN999\t\tna\t0\t0\t"
    "Test organism\t\t\tsuppressed\tScaffold\tMajor\tFull\t2010/01/01\t"
    "ASM999v1\t\t\t\t"
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/999/GCF_000009999.1_ASM999v1\n"
    "GCF_000099999.1\tPRJNA888\tSAMN888\t\tna\t0\t0\t"
    "Test organism 2\t\t\tlatest\tContig\tMajor\tFull\t2023/06/15\t"
    "ASM9999v1\t\t\t\tna\n"
)

# Relative paths for test accessions
ACC_PATH_215 = PurePosixPath("GCF") / "000" / "001" / "215" / "GCF_000001215.4_Release_6_plus_ISO1_MT"
ACC_PATH_405 = PurePosixPath("GCF") / "000" / "001" / "405" / "GCF_000001405.40_GRCh38.p14"
ACC_PATH_845 = PurePosixPath("GCF") / "000" / "005" / "845" / "GCF_000005845.2_ASM584v2"
ACC_PATH_999 = PurePosixPath("GCF") / "000" / "009" / "999" / "GCF_000009999.1_ASM999v1"

_MD5_CHECKSUMS_TXT = (
    "d41d8cd98f00b204e9800998ecf8427e  ./GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz\n"
    "a1b2c3d4e5f6a1b2c3d4e5f6a1b2c3d4  ./GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa.gz\n"
    "ffffffffffffffffffffffffffffffff  ./GCF_000001215.4_Release_6_plus_ISO1_MT_assembly_report.txt\n"
    "0000000000000000000000000000dead  ./GCF_000001215.4_Release_6_plus_ISO1_MT_README.txt\n"
)


@pytest.fixture
def manifest_transfer_mocks(monkeypatch: pytest.MonkeyPatch) -> SimpleNamespace:
    """Mocked manifest transfer functions for retrieving MD5 checksum text and listing matching S3 objects."""
    mock_retrieve = MagicMock(return_value=_MD5_CHECKSUMS_TXT)
    mock_list = MagicMock(return_value=["one key"])

    monkeypatch.setattr(manifest_mod, "ftp_retrieve_text", mock_retrieve)
    monkeypatch.setattr(manifest_mod, "list_objects", mock_list)

    return SimpleNamespace(retrieve=mock_retrieve, list=mock_list)
