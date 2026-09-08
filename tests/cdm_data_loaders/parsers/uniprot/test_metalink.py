"""Tests for UniProt metalink file parser."""

import logging
from pathlib import Path

import pytest

from cdm_data_loaders.parsers.uniprot.metalink import generate_data_source_table, get_files


def test_generate_data_source_table(test_data_dir: Path) -> None:
    """Test the parsing of the data source table."""
    metalink_path = test_data_dir / "uniprot" / "metalink" / "uniprot_metalink.xml"

    # monkeypatch the datetime
    data_src = generate_data_source_table(str(metalink_path))

    assert data_src == {
        "license": "Creative Commons Attribution 4.0 International (CC BY 4.0)",
        "publisher": "UniProt Consortium",
        "resource_type": "dataset",
        "version": "2025_04",
    }


@pytest.mark.parametrize(
    ("file", "missing"),
    [
        ("wrong_root", ["license", "publisher", "version"]),
        ("root_only", ["license", "publisher", "version"]),
        ("missing", ["license", "publisher"]),
    ],
)
def test_generate_data_source_error_missing_fields(test_data_dir: Path, file: str, missing: list[str]) -> None:
    """Test the parsing of the data source table."""
    metalink_path = test_data_dir / "uniprot" / "metalink" / "invalid" / f"{file}_metalink.xml"

    missing_fields = ", ".join(missing)
    with pytest.raises(
        RuntimeError,
        match=f"Missing required elements from metalink file: {missing_fields}",
    ):
        generate_data_source_table(str(metalink_path))


UNIPROT_METALINK_DATA = {
    "uniprot_trembl.fasta.gz": {
        "name": "uniprot_trembl.fasta.gz",
        "size": "51514683718",
        "checksum_fn": "md5",
        "checksum": "e3cd39d0c48231aa5abb3eca81b3c62a",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz",
    },
    "uniprot_sprot_varsplic.fasta.gz": {
        "name": "uniprot_sprot_varsplic.fasta.gz",
        "size": "8558360",
        "checksum_fn": "md5",
        "checksum": "70a3be6000cfacc9bd1899a47da8255b",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz",
    },
    "uniprot.xsd": {
        "name": "uniprot.xsd",
        "size": "55045",
        "checksum_fn": "md5",
        "checksum": "2a315fabb1b626efd66b4b5189e1f404",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd",
    },
    "README": {
        "name": "README",
        "size": "3918",
        "checksum_fn": "md5",
        "checksum": "c9ba41e8a3bb47c3cd69ad7f3d69810c",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/README",
    },
    "reldate.txt": {
        "name": "reldate.txt",
        "size": "151",
        "checksum_fn": "md5",
        "checksum": "8e4497ca94ff1006fea4712fa48d63a8",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/reldate.txt",
    },
    "uniprot_sprot.fasta.gz": {
        "name": "uniprot_sprot.fasta.gz",
        "size": "93100075",
        "checksum_fn": "md5",
        "checksum": "54e5460f7950fe1a66df187ca9dda1df",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
    },
    "uniprot_sprot.xml.gz": {
        "name": "uniprot_sprot.xml.gz",
        "size": "925885124",
        "checksum_fn": "md5",
        "checksum": "d277ed56e099ad7ea46b4837440f6065",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.xml.gz",
    },
    "uniprot_sprot.dat.gz": {
        "name": "uniprot_sprot.dat.gz",
        "size": "687123790",
        "checksum_fn": "md5",
        "checksum": "ecfb866a5de8f27497af396735f09b30",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz",
    },
    "LICENSE": {
        "name": "LICENSE",
        "size": "384",
        "checksum_fn": "md5",
        "checksum": "0fd67cdc8823f4bbb0b9fa0d66b8ad11",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/LICENSE",
    },
    "uniprot_trembl.xml.gz": {
        "name": "uniprot_trembl.xml.gz",
        "size": "183500710160",
        "checksum_fn": "md5",
        "checksum": "75c6580650966db18edf7f84b37d4513",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.xml.gz",
    },
    "uniprot_trembl.dat.gz": {
        "name": "uniprot_trembl.dat.gz",
        "size": "149230209630",
        "checksum_fn": "md5",
        "checksum": "0acf7438a87e86c34d85acb8fe3d1113",
        "url": "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz",
    },
}


@pytest.mark.parametrize(
    ("files_to_find", "msg"),
    [
        (["README", "uniprot.xsd"], None),
        (["uniprot.xsd"], None),
        (["file one", "file two"], "The following files were not found: file one, file two"),
        ([], "Empty file list supplied to get_files: aborting."),
        (None, None),
    ],
)
def test_get_files(
    test_data_dir: Path, files_to_find: list[str] | None, msg: str | None, caplog: pytest.LogCaptureFixture
) -> None:
    """Test retrieval of file information from a metalink file."""
    metalink_path = test_data_dir / "uniprot" / "metalink" / "uniprot_metalink.xml"
    if files_to_find is None:
        expected = UNIPROT_METALINK_DATA
    else:
        expected = {k: v for k, v in UNIPROT_METALINK_DATA.items() if k in files_to_find}

    assert get_files(metalink_path, files_to_find) == expected

    if msg:
        assert len(caplog.records) == 1
        assert caplog.records[0].levelno == logging.WARNING
        assert caplog.records[0].message == msg
    else:
        assert not caplog.records
