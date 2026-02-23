"""Tests for the gz module."""

import gzip
import logging
from pathlib import Path

import pytest
from click.testing import CliRunner

from cdm_data_loader_utils.utils.gz import compress_file, compress_files, main


@pytest.fixture
def temporary_file(tmp_path: Path) -> Path:
    """Create a temporary text file and return its path."""
    file_path = tmp_path / "sample.txt"
    content = b"Hello, world!\nThis is a test."
    file_path.write_bytes(content)
    return file_path


def read_gz(path: Path) -> bytes:
    """Utility function to read a ``.gz`` file and return its decompressed bytes."""
    with gzip.open(path, "rb") as f:
        return f.read()


def test_compress_file_creates_gzip(temporary_file: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Test that running compress file actually compresses a file."""
    # Ensure no .gz exists initially
    gz_path = temporary_file.with_suffix(temporary_file.suffix + ".gz")
    if gz_path.exists():
        gz_path.unlink()

    compress_file(temporary_file)

    # Verify the .gz file was created
    assert gz_path.exists()

    # Decompress and compare content
    decompressed = read_gz(gz_path)
    original = temporary_file.read_bytes()
    assert decompressed == original

    assert len(caplog.records) == 1
    assert caplog.records[0].levelno == logging.INFO
    assert f"Created output file {gz_path!s}" in caplog.records[0].message


def test_compress_file_skips_existing_gzip(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Ensure that an existing gz file prevents compress_file from running compression."""
    file_path = tmp_path / "data.bin"
    file_path.write_bytes(b"binary\x00data")
    # fake gzip file
    gz_path = file_path.with_suffix(file_path.suffix + ".gz")
    gz_path.write_bytes(b"some old crap")
    original_mtime = gz_path.stat().st_mtime

    # attempt to compress data.bin
    compress_file(file_path)

    # Ensure the .gz file was not overwritten (mtime unchanged)
    assert gz_path.stat().st_mtime == original_mtime
    assert gz_path.read_bytes() == b"some old crap"

    assert len(caplog.records) == 1
    assert caplog.records[0].levelno == logging.INFO
    assert f"Found existing file {file_path!s}.gz: skipping gz operation" in caplog.records[0].message


def test_compress_files_multiple(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Test that compressing a directory with a glob pattern to match works."""
    # Create three files matching the glob pattern
    n_files = 3
    for i in range(n_files):
        (tmp_path / f"file{i}.txt").write_text(f"content {i}")
    # Create a non-matching file
    (tmp_path / "ignore.md").write_text("should be ignored")

    compress_files(tmp_path, "*.txt")

    # Verify .gz files were created for the three txt files only
    for i in range(n_files):
        gz = tmp_path / f"file{i}.txt.gz"
        assert gz.exists()
        with gzip.open(gz, "rt") as f:
            assert f.read() == f"content {i}"
    # Ensure the non-matching file was not gzipped
    assert not (tmp_path / "ignore.md.gz").exists()

    assert len(caplog.records) == n_files + 2
    for r in caplog.records:
        assert r.levelno == logging.INFO
    assert "Found 3 file(s) to compress" in caplog.records[0].message
    for r in caplog.records[1:4]:
        assert f"Created output file {tmp_path!s}/file" in r.message
    assert caplog.records[-1].message == "Work complete!"


def test_compress_files_no_matches(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Test that nothing happens (no error, etc.) if a directory doesn't contain matching files."""
    compress_files(tmp_path, "*.doesnotexist")
    # No .gz files should be present
    assert list(tmp_path.rglob("*.gz")) == []

    assert len(caplog.records) == 1
    assert caplog.records[0].levelno == logging.INFO
    assert "Found 0 file(s) to compress" in caplog.records[0].message


def test_compress_files_accepts_str_path(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """The `directory` argument can be a plain string."""
    f = tmp_path / "x.txt"
    f.write_text("payload")

    # Pass directory as a string
    compress_files(str(tmp_path), "*.txt")

    # Verify outcome is the same as when a Path is supplied
    gz = f.with_name(f.name + ".gz")
    assert gz.exists()
    assert read_gz(gz) == b"payload"
    # Logger should still contain the expected messages
    assert len(caplog.records) == 3
    msgs = [rec.message for rec in caplog.records]
    assert msgs[0] == "Found 1 file(s) to compress"
    assert f"Created output file {f!s}.gz" in msgs[1]
    assert msgs[-1] == "Work complete!"


def test_compress_files_invalid_directory(tmp_path: Path) -> None:
    """Check the appropriate error is thrown if the directory is not valid."""
    non_existent = tmp_path / "no_such_dir"
    with pytest.raises(FileNotFoundError, match=f"Directory {non_existent!s} not found"):
        compress_files(non_existent, "*")


@pytest.mark.skip("CliRunner conflicts with logging, causing a ValueError")
def test_cli_single_file(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Test the CLI with a valid file."""
    # ensure that the logger does not produce any output
    caplog.set_level(logging.CRITICAL)
    runner = CliRunner()
    # Create a single file and test CLI compression of a file
    file_path = tmp_path / "single.txt"
    file_path.write_text("single file content")
    result = runner.invoke(main, ["--source", str(file_path)])
    assert result.exit_code == 0
    assert (file_path.with_suffix(file_path.suffix + ".gz")).exists()


@pytest.mark.skip("CliRunner conflicts with logging, causing a ValueError")
def test_cli_directory_with_glob(tmp_path: Path, caplog: pytest.LogCaptureFixture) -> None:
    """Test the CLI with a directory and some globbery action."""
    caplog.set_level(logging.CRITICAL)
    sub_dir = tmp_path / "sub"
    sub_dir.mkdir()
    for name in ["a.txt", "b.log"]:
        (sub_dir / name).write_text(name)
    runner = CliRunner()
    result = runner.invoke(main, ["--source", str(sub_dir), "--file-glob", "*.txt"])
    assert result.exit_code == 0
    # Only the .txt file should be gzipped
    assert (sub_dir / "a.txt.gz").exists()
    assert not (sub_dir / "b.log.gz").exists()


def test_cli_main_nonexistent_source(tmp_path: Path) -> None:
    """Test the CLI interface with a missing source."""
    runner = CliRunner()
    result = runner.invoke(main, ["--source", str(tmp_path / "doesnotexist")])
    assert result.exit_code != 0
    assert isinstance(result.exception, RuntimeError)
