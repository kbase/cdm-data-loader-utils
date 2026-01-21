"""File compression utils."""

import gzip
import shutil
from pathlib import Path

import click

from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger

logger = get_cdm_logger()


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
        raise ValueError(msg)

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
        logger.info("Found existing file %s: skipping gz operation", str(file) + ".gz")
        return

    with file.open("rb") as f_in, gzip.open(str(file) + ".gz", "wb") as f_out:
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
