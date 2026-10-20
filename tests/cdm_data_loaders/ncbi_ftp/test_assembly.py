"""Tests for ncbi_ftp.assembly module — path helpers, file filtering, checksum parsing."""

from pathlib import PurePosixPath

import pytest

from cdm_data_loaders.ncbi_ftp.assembly import (
    build_accession_path,
    parse_assembly_path,
    parse_md5_checksums_file,
)

from .conftest import ACC_PATH_215

# Path helpers


@pytest.mark.parametrize(
    ("assembly_dir", "expected"),
    [
        pytest.param(
            "GCF_000001215.4_Release_6_plus_ISO1_MT",
            PurePosixPath("raw_data") / ACC_PATH_215,
            id="gcf",
        ),
        pytest.param(
            "GCA_012345678.1_ASM1234v1",
            PurePosixPath("raw_data") / "GCA" / "012" / "345" / "678" / "GCA_012345678.1_ASM1234v1",
            id="gca",
        ),
    ],
)
def test_build_accession_path(assembly_dir: PurePosixPath, expected: str) -> None:
    """Verify accession path construction for various inputs."""
    assert build_accession_path(assembly_dir) == expected


def test_build_accession_path_invalid() -> None:
    """Verify ValueError on invalid assembly name."""
    with pytest.raises(ValueError, match="Cannot parse"):
        build_accession_path(PurePosixPath("invalid_name"))


@pytest.mark.parametrize(
    ("path", "expected"),
    [
        pytest.param(
            PurePosixPath("/") / "genomes" / "all" / ACC_PATH_215,
            ("GCF", PurePosixPath("GCF_000001215.4_Release_6_plus_ISO1_MT"), "GCF_000001215.4"),
            id="with_trailing_slash",
        ),
        pytest.param(
            PurePosixPath("/") / "genomes" / "all" / ACC_PATH_215,
            ("GCF", PurePosixPath("GCF_000001215.4_Release_6_plus_ISO1_MT"), "GCF_000001215.4"),
            id="without_trailing_slash",
        ),
    ],
)
def test_parse_assembly_path(path: PurePosixPath, expected: tuple[str, PurePosixPath, str]) -> None:
    """Verify db, assembly_dir, and accession are parsed correctly."""
    assert parse_assembly_path(path) == expected


def test_parse_assembly_path_invalid() -> None:
    """Verify ValueError on invalid path."""
    with pytest.raises(ValueError, match="Cannot parse"):
        parse_assembly_path(PurePosixPath("/") / "random" / "path")


# parse_md5_checksums_file


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        pytest.param(
            "abc123  ./GCF_000001215.4_genomic.fna.gz\ndef456  ./GCF_000001215.4_genomic.gff.gz\n",
            {
                PurePosixPath("GCF_000001215.4_genomic.fna.gz"): "abc123",
                PurePosixPath("GCF_000001215.4_genomic.gff.gz"): "def456",
            },
            id="dot_slash_prefix",
        ),
        pytest.param(
            "abc123  GCF_000001215.4_genomic.fna.gz\n",
            {PurePosixPath("GCF_000001215.4_genomic.fna.gz"): "abc123"},
            id="no_dot_slash_prefix",
        ),
        pytest.param("", {}, id="empty_string"),
        pytest.param("   \n  \n", {}, id="whitespace_only"),
        pytest.param(
            "abc123  file1.txt\n\n\ndef456  file2.txt\n",
            {PurePosixPath("file1.txt"): "abc123", PurePosixPath("file2.txt"): "def456"},
            id="blank_lines_ignored",
        ),
    ],
)
def test_parse_md5_checksums_file(text: str, expected: dict[PurePosixPath, str]) -> None:
    """Verify parse_md5_checksums_file handles various input formats."""
    assert parse_md5_checksums_file(text) == expected
