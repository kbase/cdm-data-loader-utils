"""Tests for the UniProt relnotes (release notes) parser."""

import datetime
from pathlib import Path

import pytest

from cdm_data_loaders.parsers.uniprot.relnotes import parse, parse_relnotes


def test_parse_relnotes(test_data_dir: Path) -> None:
    """Test parsing some generic release notes."""
    relnotes_path = test_data_dir / "uniprot" / "relnotes_2025_03.txt"

    parsed = parse_relnotes(relnotes_path)
    assert parsed == {
        "version": "2025_03",
        "date_published": datetime.datetime(2025, 6, 18, 0, 0, 0),  # noqa: DTZ001
        "UniProtKB/Swiss-Prot": "573,661",
        "UniProtKB/TrEMBL": "253,061,697",
        "UniProtKB": "253,635,358",
        "UniRef100": "465,330,530",
        "UniRef90": "208,005,650",
        "UniRef50": "70,198,728",
        "UniParc": "982,121,738",
    }


def test_format_change() -> None:
    """Test parsing a file with the important intro stuff (version/date information) missing."""
    content = """# Some unrelated header line\nUniProtKB Release 2025_03 consists of 1,000 entries (UniProtKB/Swiss-Prot:\n100 entries and UniProtKB/TrEMBL: 900 entries)\nUniRef100 Release 2025_03 consists of 5,000 entries"""

    with pytest.raises(RuntimeError, match=r"Could not find double line break. Relnotes file format may have changed"):
        parse(content)


def test_missing_intro_regex() -> None:
    """Test parsing a file with the important intro stuff (version/date information) missing."""
    content = """# Some unrelated header line\n\nUniProtKB Release 2025_03 consists of 1,000 entries (UniProtKB/Swiss-Prot:\n100 entries and UniProtKB/TrEMBL: 900 entries)\n\nUniRef100 Release 2025_03 consists of 5,000 entries"""

    with pytest.raises(RuntimeError) as excinfo:
        parse(content)
    assert "Could not find text matching the release version date regex" in str(excinfo.value)


def test_missing_trembl_stats() -> None:
    """Test parsing a file where the UniProt/TrEMBL stats block is malformed."""
    content = """# The UniProt consortium European Bioinformatics Institute (EBI), SIB Swiss Institute of Bioinformatics and Protein Information Resource (PIR), is pleased to announce UniProt Knowledgebase (UniProtKB) Release 2025_03 (18-Jun-2025).\n\n# Some other unrelated content that does not match the expected UNIPROT_TREMBL_STATS pattern\n\nUniRef100 Release 2025_03 consists of 5,000 entries"""

    with pytest.raises(RuntimeError) as excinfo:
        parse(content)
    # Both the missing TrEMBL stats and the missing release stats (UniRef) should be reported
    msg = str(excinfo.value)
    assert "Could not find text matching the UniProt/TrEMBL stats regex" in msg
    # UniRef stats are also missing
    assert "No stats for UniRef50 found" in msg
    assert "No stats for UniRef90 found" in msg


def test_missing_release_stats_and_uniref() -> None:
    """Test parsing a file where there are no release stats for UniRefs or UniParc."""
    content = """# The UniProt consortium European Bioinformatics Institute (EBI), SIB Swiss Institute of Bioinformatics and Protein Information Resource (PIR), is pleased to announce UniProt Knowledgebase (UniProtKB) Release 2025_03.\n\nUniProtKB Release 2025_03 consists of 1,000 entries (UniProtKB/Swiss-Prot:\n100 entries and UniProtKB/TrEMBL: 900 entries)\n\n# No further release stats"""
    with pytest.raises(RuntimeError) as excinfo:
        parse(content)
    msg = str(excinfo.value)
    # The parser should flag missing release stats regex and then report each missing UniRef
    assert "Could not find text matching the release version date regex." in msg
    assert "No stats for UniRef50 found" in msg
    assert "No stats for UniRef90 found" in msg
    assert "No stats for UniRef100 found" in msg
