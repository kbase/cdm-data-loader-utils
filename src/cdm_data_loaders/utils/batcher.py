"""Utilities for splitting a selection into batches."""

import re
from collections.abc import Generator
from logging import Logger, getLogger
from pathlib import Path
from typing import Annotated, Any, Final, Self

from pydantic import BaseModel, DirectoryPath, Field, model_validator

from cdm_data_loaders.core.fields import MIN_START_AT
from cdm_data_loaders.core.settings import CtsSettings

# Matches files like: name_00001.ext or name_00001.ext.gz
FILE_NAME_REGEX: re.Pattern[str] = re.compile(r"^\w+_(\d+)(\.\w+)+$")

MIN_BATCH_SIZE: Final[int] = 1
MIN_END_AT: Final[int] = 0
DEFAULT_PIPELINE_BATCH_SIZE: Final[int] = 50

logger: Logger = getLogger(__name__)


class NumericFileSequenceBatcher(BaseModel):
    """A batcher that can be used to retrieve batches of files from a directory.

    Searches the directory for files with names ending in `_[0-9]+`
    (underscore followed by one or more digits) with one or more extensions (e.g. `.gz`, `.txt`, `.tar.gz`).

    Assumes that file names only contain alphanumeric characters, underscore, and extension(s).

    Re-scans the directory each time the next batch is requested.

    :param directory:   directory to retrieve files from; strings will be coerced to Paths
    :type  directory:   Path, required
    :param batch_size:  number of files to return per invocation, defaults to 1
    :type  batch_size:  int, optional
    :param start_at:    file number to start at, defaults to 1.
                        Must be greater than or equal to 1.
    :type  start_at:    int, optional
    :param end_at:      file number to end at, inclusive (i.e. if set to 15, file_0015.txt will be the last file)
                        If set to 0, the end_at parameter is ignored -- i.e. there is no maximum value.
                        Must be greater than or equal to 0
                        Defaults to 0
    :type  end_at:      int, optional
    :param file_regex:  pattern that the file names should match. Capture group 1 is expected to be numeric.
    :type  file_regex:  re.Pattern, optional
    """

    directory: DirectoryPath
    batch_size: Annotated[int, Field(ge=MIN_BATCH_SIZE, strict=True)] = MIN_BATCH_SIZE
    start_at: Annotated[int, Field(ge=MIN_START_AT, strict=True)] = MIN_START_AT
    end_at: Annotated[int, Field(ge=MIN_END_AT, strict=True)] = MIN_END_AT
    file_regex: re.Pattern[str] = Field(default=FILE_NAME_REGEX)

    @model_validator(mode="after")
    def check_start_at_end_at(self) -> Self:
        """Ensure that the end_at value, if set and non-zero, is greater than the start_at value.

        :raises ValueError: if end_at is set and start_at exceeds end_at
        :return: self
        :rtype: Self
        """
        if self.end_at != 0 and self.end_at < self.start_at:
            err_msg = "end_at must be greater than start_at"
            raise ValueError(err_msg)
        return self

    def _get_sequence_number(self, path: Path) -> int:
        match = self.file_regex.match(path.name)
        return int(match.group(1))  # pyright: ignore[reportOptionalMemberAccess]

    def get_batch(self) -> list[Path]:
        """Return the next `batch_size` files whose sequence number >= start_at.

        Re-scans the directory on every call to pick up newly added files and
        updates `start_at` to the next file in the directory list.
        """
        if self.end_at and self.start_at > self.end_at:
            return []

        matched = sorted(p for p in self.directory.iterdir() if p.is_file() and self.file_regex.match(p.name))
        if not matched:
            logger.warning("No matching files found in %s", str(self.directory))
            return []

        eligible = [
            p
            for p in matched
            if self._get_sequence_number(p) >= self.start_at
            and (not self.end_at or self._get_sequence_number(p) <= self.end_at)
        ]

        batch = eligible[: self.batch_size]
        if batch:
            self.start_at = self._get_sequence_number(batch[-1]) + 1

        return batch


def get_file_batches(settings: CtsSettings) -> Generator[list[Path], Any]:
    """Yield successive batches of files to process, driven by a NumericFileSequenceBatcher.

    :param settings: pipeline config with input_dir
    :type  settings: CtsSettings
    :yield: batches of file paths
    :rtype: Generator[list[Path], Any]
    """
    batcher = NumericFileSequenceBatcher(directory=settings.input_dir, batch_size=DEFAULT_PIPELINE_BATCH_SIZE)
    while files := batcher.get_batch():
        logger.debug("Files to be processed:%s", "".join("- " + str(f) + "\n" for f in files))
        yield files
