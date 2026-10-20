"""File compression utils."""

import gzip
import shutil
from collections.abc import Generator
from contextlib import contextmanager
from io import BytesIO
from logging import Logger, getLogger
from pathlib import Path
from typing import BinaryIO

import click

from cdm_data_loaders.core.fields import GZIP_SUFFIX

logger: Logger = getLogger(__name__)


@contextmanager
def open_maybe_gzip(file: str | Path, mode: str = "rb") -> Generator[BinaryIO | BytesIO]:
    """Open a file with gzip transparency for binary reads.

    Opens the file with gzip for names ending in .gz, and normally otherwise.
    Only binary read modes are supported.

    :param file: file to open; the file can be gzipped or not.
    :type file: str | Path
    :param mode: mode to open the file with; must be a binary read mode.
    :type mode: str
    :raises ValueError: if mode is not a binary read mode.
    :yield: an open binary file handle.
    :rtype: Generator[BinaryIO | BytesIO, None, None]
    """
    if mode != "rb":
        err_msg = f"open_maybe_gzip only supports binary reads, got mode: {mode}"
        raise ValueError(err_msg)

    path = file if isinstance(file, Path) else Path(file)
    if path.suffix == GZIP_SUFFIX:
        with gzip.open(path, "rb") as f_in:
            yield f_in
        return
    with path.open("rb") as f_in:
        yield f_in


def decompress_file(file: str | Path) -> None:
    """Decompress a gzip file.

    :param file: file to decompress
    :type file: str | Path
    """
    if not isinstance(file, Path):
        file = Path(file)

    if file.suffix != GZIP_SUFFIX:
        logger.info("File %s does not end with .gz: skipping decompression", str(file))
        return

    if not file.exists():
        logger.warning("File %s does not exist: skipping decompression", str(file))
        return

    if not file.is_file():
        logger.warning("%s is not a file: skipping decompression", str(file))
        return

    # remove .gz suffix
    output_file = file.with_suffix("")
    if output_file.exists():
        logger.info("Found existing file %s: skipping decompression", output_file)
        return

    with gzip.open(file, "rb") as f_in, output_file.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    logger.info("Created output file %s", output_file)


def compress_files(directory: Path | str, file_glob: str) -> None:
    """Compress all files matching a certain pattern in a directory.

    :param directory: directory to look in
    :type directory: Path | str
    :param file_glob: pattern to match
    :type file_glob: str
    :raises ValueError: if the directory is not found
    """
    if not isinstance(directory, Path):
        directory = Path(directory)

    if not directory.exists() or not directory.is_dir():
        msg = f"Directory {directory!s} not found: check the path is correct"
        raise FileNotFoundError(msg)

    to_compress = list(directory.glob(file_glob))
    logger.info("Found %d file(s) to compress", len(to_compress))
    if not to_compress:
        return

    for f in to_compress:
        compress_file(f)
    logger.info("Work complete!")


def compress_file(file: Path) -> None:
    """Compress a file using gzip.

    :param file: file to compress
    :type file: Path
    """
    if Path(f"{file!s}.gz").exists():
        logger.info("Found existing file %s: skipping gz operation", str(file) + GZIP_SUFFIX)
        return

    with file.open("rb") as f_in, gzip.open(str(file) + GZIP_SUFFIX, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    logger.info("Created output file %s.gz", str(file))


@click.command()
@click.option("--source", "-i", required=True, help="input file(s) to process")
@click.option("--file-glob", "-f", default="*", help="glob for files in a directory")
def main(source: str, file_glob: str) -> None:
    """Compress a file or directory from the command line."""
    source_path = Path(source)
    if not source_path.exists():
        msg = f"Path {source} does not exist"
        raise RuntimeError(msg)

    if source_path.is_dir():
        compress_files(source_path, file_glob or "*")
    else:
        compress_file(source_path)


if __name__ == "__main__":
    main()
