"""Tests of file system-related utilities."""

import logging
from collections.abc import Callable
from copy import deepcopy
from pathlib import Path
from typing import Any
from unittest.mock import MagicMock

import pytest
from pydantic import ValidationError

import cdm_data_loaders.utils.batcher
from cdm_data_loaders.core.fields import DEFAULT_PIPELINE_BATCH_SIZE, MIN_START_AT
from cdm_data_loaders.utils.batcher import (
    FILE_NAME_REGEX,
    MIN_BATCH_SIZE,
    MIN_END_AT,
    NumericFileSequenceBatcher,
    get_file_batches,
)

# the maximum file number in the directory
MAX_FILE_NUMBER = 15
# what the range should be set to (max range is exclusive)
MAX_RANGE_VALUE = 16

EXPECTED: dict[int | None, dict[int | None, Any]] = {
    # batch_size
    1: {
        # start_at
        1: [[r] for r in range(1, MAX_RANGE_VALUE)],
        6: [[r] for r in range(6, MAX_RANGE_VALUE)],
        8: [[r] for r in range(8, MAX_RANGE_VALUE)],
        11: [[r] for r in range(11, MAX_RANGE_VALUE)],
    },
    5: {
        # start_at
        1: [range(1, 6), range(6, 11), range(11, MAX_RANGE_VALUE)],
        6: [range(6, 11), range(11, MAX_RANGE_VALUE)],
        8: [range(8, 13), range(13, MAX_RANGE_VALUE)],
        11: [range(11, MAX_RANGE_VALUE)],
    },
    8: {
        1: [range(1, 9), range(9, MAX_RANGE_VALUE)],
        6: [range(6, 14), range(14, MAX_RANGE_VALUE)],
        8: [range(8, MAX_RANGE_VALUE)],
        11: [range(11, MAX_RANGE_VALUE)],
    },
    15: {
        1: [range(1, MAX_RANGE_VALUE)],
        6: [range(6, MAX_RANGE_VALUE)],
        8: [range(8, MAX_RANGE_VALUE)],
        11: [range(11, MAX_RANGE_VALUE)],
    },
}
# batch_size is not specified
EXPECTED[None] = EXPECTED[1]
# batch_size greater than # of records
EXPECTED[20] = EXPECTED[15]
# add in results where start_at is not specified (i.e. uses the default value)
for ix, vals in EXPECTED.items():
    if ix is not None:
        EXPECTED[ix][None] = vals[1]

EXPECTED_END_AT = deepcopy(EXPECTED)


def make_files(directory: Path, names: list[str]) -> list[Path]:
    """Touch each filename in *directory* and return the sorted Path list."""
    paths = []
    for name in names:
        p = directory / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.touch()
        paths.append(p)
    return sorted(paths)


def make_file_names(prefix: str, ext: str, numbers: list[int] | range) -> list[str]:
    """Create the file names for a given set of numbers."""
    if not isinstance(numbers, list):
        # convert the range into a list
        numbers = list(numbers)
    return [f"{prefix}_{n:05}.{ext}" for n in numbers]


def make_sequence(directory: Path, prefix: str, ext: str, numbers: list[int] | range) -> list[Path]:
    """Create files for the given sequence numbers and return sorted Paths."""
    names = make_file_names(prefix, ext, numbers)
    return make_files(directory, names)


@pytest.fixture
def file_dir(tmp_path: Path) -> Path:
    """Directory pre-populated with files numbered 00001 - 00015."""
    make_sequence(tmp_path, "report", "csv", range(1, MAX_RANGE_VALUE))
    return tmp_path


@pytest.fixture(scope="module")
def non_matching_file_names() -> list[str]:
    """List of file names that do not match the pattern required by NumericFileSequenceBatcher."""
    return [
        # no numbers
        "README.md",
        # contains non-\w character
        ".hidden_00001.txt",
        # no extension
        "data_00001",
        # files in nested dirs -- will not be found
        "nested/data_00010.txt",
        "nested/dir1/data_00020.txt",
    ]


@pytest.fixture
def mixed_dir(tmp_path: Path, non_matching_file_names: list[str]) -> Path:
    """Directory containing valid files alongside files that should be ignored."""
    make_sequence(tmp_path, "data", "txt", range(1, 6))
    more_files = [
        "data_123.txt",
        "data_000001.txt",
        "data_000100.txt.gz",
        "data_000200.txt.tar.gz",
        "data_000400.csv.gz",
        "file_000300.txt.tar.gz",
        *non_matching_file_names,
    ]
    make_files(tmp_path, more_files)
    return tmp_path


@pytest.fixture
def fake_batcher_factory() -> Callable[..., tuple[type, list[dict[str, Any]]]]:
    """Return a factory producing a fake NumericFileSequenceBatcher class fed by a fixed batch queue."""

    def _make(batches: list[list[Path]]) -> tuple[type, list[dict[str, Any]]]:
        construction_calls: list[dict[str, Any]] = []
        batch_queue = list(batches)

        class _FakeBatcher:
            def __init__(self, **kwargs: Any) -> None:
                construction_calls.append(kwargs)

            def get_batch(self) -> list[Path]:
                return batch_queue.pop(0) if batch_queue else []

        return _FakeBatcher, construction_calls

    return _make


def test_init_batcher_defaults(tmp_path: Path) -> None:
    """Ensure defaults are set correctly."""
    bc = NumericFileSequenceBatcher(directory=str(tmp_path))  # type: ignore[reportArgumentType]
    assert bc.batch_size == MIN_BATCH_SIZE
    assert bc.start_at == MIN_START_AT
    assert bc.end_at == MIN_END_AT
    assert bc.file_regex == FILE_NAME_REGEX
    assert bc.directory == tmp_path


def test_init_batcher_values(tmp_path: Path) -> None:
    """Ensure defaults are set correctly."""
    start_at = 26
    batch_size = 10
    end_at = 26
    bc = NumericFileSequenceBatcher(directory=tmp_path, start_at=start_at, batch_size=batch_size, end_at=end_at)
    assert bc.start_at == start_at
    assert bc.batch_size == batch_size
    assert bc.end_at == end_at
    assert bc.file_regex == FILE_NAME_REGEX
    assert bc.directory == tmp_path


def test_init_batcher_no_directory() -> None:
    """Ensure the batcher throws an error if a directory is not specified."""
    with pytest.raises(ValidationError, match=r"directory\s+Field required") as err:
        NumericFileSequenceBatcher()  # type: ignore[reportCallIssue]

    assert err.match(r"directory\s+Field required")


VALID_INTEGER = r"\s+Input should be a valid integer"
GTE_ONE = r"\s+Input should be greater than or equal to 1"
GTE_ZERO = r"\s+Input should be greater than or equal to 0"

GTE_ONE_TEST_VALUES = [
    (None, VALID_INTEGER),
    (0, GTE_ONE),
    (-1, GTE_ONE),
    (1.0, VALID_INTEGER),
    (1.2345678, VALID_INTEGER),
    (-15, GTE_ONE),
    (-1234567890, GTE_ONE),
    ("something", VALID_INTEGER),
    ("50", VALID_INTEGER),
    ("None", VALID_INTEGER),
]


@pytest.mark.parametrize(("batch_size", "err_msg"), GTE_ONE_TEST_VALUES)
def test_init_batcher_invalid_batch_size(batch_size: float | str | None, err_msg: str) -> None:
    """Test invalid batch_size, start_at, end_at, and file_regex parameters."""
    with pytest.raises(ValidationError, match=f"batch_size{err_msg}"):
        NumericFileSequenceBatcher(directory=".", batch_size=batch_size)  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize(("start_at", "err_msg"), GTE_ONE_TEST_VALUES)
def test_init_batcher_invalid_start_at_params(start_at: float | None, err_msg: str) -> None:
    """Test invalid start_at parameters."""
    with pytest.raises(ValidationError, match=f"start_at{err_msg}"):
        NumericFileSequenceBatcher(directory=".", start_at=start_at)  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize(
    ("end_at", "err_msg"),
    [
        (None, VALID_INTEGER),
        (-1, GTE_ZERO),
        (1.0, VALID_INTEGER),
        (1.2345678, VALID_INTEGER),
        (-15, GTE_ZERO),
        (-1234567890, GTE_ZERO),
        ("something", VALID_INTEGER),
        ("50", VALID_INTEGER),
        ("None", VALID_INTEGER),
    ],
)
def test_init_batcher_invalid_end_at_params(end_at: str | float | None, err_msg: str) -> None:
    """Test invalid end_at parameters."""
    with pytest.raises(ValidationError, match=f"end_at{err_msg}"):
        NumericFileSequenceBatcher(directory=".", end_at=end_at)  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize(
    ("batch_size", "batch_size_err_msg"),
    [
        (None, VALID_INTEGER),
        (0, GTE_ONE),
        (-1, GTE_ONE),
        (1.0, VALID_INTEGER),
        ("50", VALID_INTEGER),
    ],
)
@pytest.mark.parametrize(
    ("start_at", "start_at_err_msg"),
    [
        (None, VALID_INTEGER),
        (0, GTE_ONE),
        (1.0, VALID_INTEGER),
        (-1234567890, GTE_ONE),
        ("50", VALID_INTEGER),
    ],
)
@pytest.mark.parametrize(
    ("end_at", "end_at_err_msg"),
    [
        (None, VALID_INTEGER),
        (-1, GTE_ZERO),
        (1.2345678, VALID_INTEGER),
        (-1234567890, GTE_ZERO),
        ("50", VALID_INTEGER),
    ],
)
def test_init_batcher_multiple_errors(
    batch_size: str | float | None,
    start_at: str | float | None,
    end_at: str | float | None,
    batch_size_err_msg: str,
    start_at_err_msg: str,
    end_at_err_msg: str,
) -> None:
    """Ensure that supplying the NumericFileSequenceBatcher with lots of dodgy values triggers an error."""
    with pytest.raises(ValidationError, match="3 validation errors for NumericFileSequenceBatcher") as exc:
        NumericFileSequenceBatcher(directory=Path(), start_at=start_at, end_at=end_at, batch_size=batch_size)  # type: ignore[reportArgumentType]
    assert exc.match(f"start_at{start_at_err_msg}")
    assert exc.match(f"end_at{end_at_err_msg}")
    assert exc.match(f"batch_size{batch_size_err_msg}")


@pytest.mark.parametrize(("start_at", "end_at"), [(5, 3), (123456789, 123456788)])
def test_init_batcher_invalid_start_vs_end_at_params(start_at: int, end_at: int) -> None:
    """Ensure that an error is thrown if start_at and end_at are not compatible."""
    with pytest.raises(ValidationError, match="end_at must be greater than start_at"):
        NumericFileSequenceBatcher(directory=".", start_at=start_at, end_at=end_at)  # pyright: ignore[reportArgumentType]


@pytest.mark.parametrize(("start_at", "end_at"), [(1, 0), (1, 1), (20, 20), (5, 0), (3, 7)])
def test_init_batcher_valid_start_vs_end_at_params(start_at: int, end_at: int) -> None:
    """Ensure that 0 is a valid end_at parameter, regardless of start_at value."""
    bc = NumericFileSequenceBatcher(directory=".", start_at=start_at, end_at=end_at)  # pyright: ignore[reportArgumentType]
    assert bc.start_at == start_at
    assert bc.end_at == end_at


def test_end_at_greater_than_start_at_during_iteration() -> None:
    """Ensure that if end_at is smaller than start_at during iteration, an empty list is returned."""
    bc = NumericFileSequenceBatcher(directory=".", start_at=MIN_START_AT, end_at=5)  # pyright: ignore[reportArgumentType]
    assert bc.end_at == 5  # noqa: PLR2004
    bc.start_at = 10
    assert bc.get_batch() == []


CUTOFF_VALUE = 12


# your basic batch
@pytest.mark.parametrize("end_at", [None, 0, CUTOFF_VALUE])
@pytest.mark.parametrize("start_at", [None, 1, 6, 8, 11])
@pytest.mark.parametrize("batch_size", EXPECTED.keys())
def test_get_batch_parametrized(
    file_dir: Path,
    batch_size: int | None,
    start_at: int | None,
    end_at: int | None,
) -> None:
    """Test retrieval of batches of files."""
    cursor_params = {}
    if batch_size is not None:
        cursor_params["batch_size"] = batch_size
    if start_at is not None:
        cursor_params["start_at"] = start_at
    if end_at is not None:
        cursor_params["end_at"] = end_at

    cursor = NumericFileSequenceBatcher(directory=file_dir, **cursor_params)

    # generate the expected files
    expected_files: list[list[Path]] = [
        [file_dir / fn for fn in make_file_names("report", "csv", numbers)]
        for numbers in EXPECTED[batch_size][start_at]
    ]
    if end_at:
        expected_files = []
        for numbers in EXPECTED[batch_size][start_at]:
            if cutoffless := [n for n in numbers if n <= end_at]:
                expected_files.append([file_dir / fn for fn in make_file_names("report", "csv", cutoffless)])

    output: list[list[Path]] = []
    while batch := cursor.get_batch():
        output.append(batch)
        if cursor.start_at >= MAX_RANGE_VALUE:
            break

    # check the number of batches is correct
    assert len(output) == len(expected_files)

    # results are sorted
    for batch in output:
        # results are all file paths
        assert all(isinstance(p, Path) for p in batch)
        assert sorted(batch) == batch
    assert output == expected_files

    # if end_at is defined, the start_at value will be one greater than the end_at value
    if end_at:
        assert cursor.start_at == CUTOFF_VALUE + 1


def test_get_batch_default_start_at_is_zero(file_dir: Path) -> None:
    """Ensure that the default start_at is 0."""
    cursor_default = NumericFileSequenceBatcher(directory=file_dir, batch_size=3)
    cursor_explicit = NumericFileSequenceBatcher(directory=file_dir, batch_size=3, start_at=MIN_START_AT)
    assert cursor_default.get_batch() == cursor_explicit.get_batch()


def test_get_batch_start_at_matches_sequence_number(file_dir: Path) -> None:
    """Ensure start_at value matches sequence number."""
    cursor = NumericFileSequenceBatcher(directory=file_dir, batch_size=5, start_at=15)
    result = cursor.get_batch()
    assert len(result) == 1
    assert result[0].name == "report_00015.csv"


# advancing the cursor
def test_get_batch_start_at_advances_after_get_batch(file_dir: Path) -> None:
    """Ensure that the start_at value changes after each successful get_batch operation."""
    batch_size = 5
    cursor = NumericFileSequenceBatcher(directory=file_dir, batch_size=batch_size, start_at=MIN_START_AT)
    assert cursor.start_at == MIN_START_AT
    batch_1 = cursor.get_batch()
    assert cursor.start_at == batch_size + 1  # next file is report_00006.csv
    batch_2 = cursor.get_batch()
    assert cursor.start_at == batch_size * 2 + 1  # report_00011.csv
    batch_3 = cursor.get_batch()
    assert cursor.start_at == batch_size * 3 + 1  # report_00016.csv (does not exist)
    # next call returns nothing
    assert cursor.get_batch() == []

    # all files should be the sequential list of existing files
    all_files = batch_1 + batch_2 + batch_3
    assert all_files == [file_dir / f"report_{n:05}.csv" for n in range(1, MAX_RANGE_VALUE)]


def test_get_batch_cursor_does_not_advance_on_empty_result(file_dir: Path) -> None:
    """Ensure that the cursor does not advance if the batch is empty."""
    start_at = 999
    cursor = NumericFileSequenceBatcher(directory=file_dir, batch_size=5, start_at=start_at)
    cursor.get_batch()
    assert cursor.start_at == start_at


def test_get_batch_partial_batch_advances_correctly(file_dir: Path) -> None:
    """Ensure that the cursor only advances as far as the last file in the batch."""
    # Only 3 files remain from 13 onward
    cursor = NumericFileSequenceBatcher(directory=file_dir, batch_size=5, start_at=13)
    result = cursor.get_batch()
    assert result == [file_dir / f"report_{n:05}.csv" for n in [13, 14, 15]]
    assert cursor.start_at == 16  # noqa: PLR2004


def test_get_batch_cursor_can_be_reset(file_dir: Path) -> None:
    """Ensure that the cursor can be reset."""
    batch_size = 5
    cursor = NumericFileSequenceBatcher(directory=file_dir, batch_size=batch_size)
    original_result = cursor.get_batch()
    assert cursor.start_at == batch_size + 1
    # set cursor to 0
    cursor.start_at = 1
    reset_result = cursor.get_batch()
    assert cursor.start_at == batch_size + 1
    assert original_result == reset_result
    assert reset_result[0].name == "report_00001.csv"


# Edge cases -- boundaries
def test_get_batch_start_at_beyond_end_returns_empty_list(file_dir: Path) -> None:
    """Ensure that nothing is returned if start_at is too high."""
    cursor = NumericFileSequenceBatcher(directory=file_dir, batch_size=5, start_at=999)
    assert cursor.get_batch() == []


def test_get_batch_empty_directory_returns_empty_list(tmp_path: Path) -> None:
    """Ensure that an empty dir returns nothing."""
    cursor = NumericFileSequenceBatcher(directory=tmp_path, batch_size=5)
    assert cursor.get_batch() == []


def test_get_batch_batch_size_larger_than_remaining_files(file_dir: Path) -> None:
    """Ensure that batches are sized correctly for partial batches."""
    cursor = NumericFileSequenceBatcher(directory=file_dir, batch_size=10, start_at=10)
    result = cursor.get_batch()
    # should have 00010 - 00015
    assert result[-1].name == "report_00015.csv"
    assert result == [file_dir / f"report_{n:05}.csv" for n in range(10, MAX_RANGE_VALUE)]
    assert cursor.start_at == 16  # noqa: PLR2004


# gaps in the sequence
def test_get_batch_start_at_skips_to_next_available_when_gap(tmp_path: Path) -> None:
    """Ensure that gaps in the sequence are dealt with correctly."""
    # Files exist for 1,2,3 then jump to 10,11,12 — no 4-9
    make_sequence(tmp_path, "data", "csv.gz", [1, 2, 3, 10, 11, 12])
    cursor = NumericFileSequenceBatcher(directory=tmp_path, batch_size=5, start_at=5)
    # retrieve 5 files, starting at 00005
    assert cursor.get_batch() == [tmp_path / f"data_{n:05}.csv.gz" for n in [10, 11, 12]]
    assert cursor.start_at == 13  # noqa: PLR2004
    assert cursor.get_batch() == []


def test__get_batchsequential_calls_across_gap(tmp_path: Path) -> None:
    """Ensure that files are correctly retrieved across gaps in the sequence."""
    make_sequence(tmp_path, "data", "csv.gz", [1, 2, 3, 10, 11, 12])
    cursor = NumericFileSequenceBatcher(directory=tmp_path, batch_size=2)
    assert cursor.get_batch() == [tmp_path / f"data_{n:05}.csv.gz" for n in [1, 2]]
    assert cursor.start_at == 3  # noqa: PLR2004

    assert cursor.get_batch() == [tmp_path / f"data_{n:05}.csv.gz" for n in [3, 10]]
    assert cursor.start_at == 11  # noqa: PLR2004

    assert cursor.get_batch() == [tmp_path / f"data_{n:05}.csv.gz" for n in [11, 12]]
    assert cursor.start_at == 13  # noqa: PLR2004

    assert cursor.get_batch() == []
    assert cursor.start_at == 13  # noqa: PLR2004


def test_get_batch_sequential_calls_across_gap_with_end_at(tmp_path: Path) -> None:
    """Ensure that files are correctly retrieved across gaps in the sequence when end_at is specified."""
    make_sequence(tmp_path, "data", "csv.gz", [1, 2, 3, 5, 8, 11, 15])
    cursor = NumericFileSequenceBatcher(directory=tmp_path, batch_size=2, end_at=10)
    assert cursor.get_batch() == [tmp_path / f"data_{n:05}.csv.gz" for n in [1, 2]]
    assert cursor.start_at == 3  # noqa: PLR2004

    assert cursor.get_batch() == [tmp_path / f"data_{n:05}.csv.gz" for n in [3, 5]]
    assert cursor.start_at == 6  # noqa: PLR2004

    assert cursor.get_batch() == [tmp_path / f"data_{n:05}.csv.gz" for n in [8]]
    assert cursor.start_at == 9  # noqa: PLR2004

    assert cursor.get_batch() == []
    assert cursor.start_at == 9  # noqa: PLR2004


# File-name pattern filtering
def test_get_batch_ignores_invalid_filenames(mixed_dir: Path) -> None:
    """Ensure that filenames are matched correctly."""
    cursor = NumericFileSequenceBatcher(directory=mixed_dir, batch_size=20)
    generated_file_names = [f"data_{n:05}.txt" for n in range(1, 6)]
    file_names = sorted(
        [
            *generated_file_names,
            "data_123.txt",
            "data_000001.txt",
            "data_000100.txt.gz",
            "data_000200.txt.tar.gz",
            "data_000400.csv.gz",
            "file_000300.txt.tar.gz",
        ]
    )
    assert cursor.get_batch() == [mixed_dir / fn for fn in file_names]


def test_get_batch_mixed_extensions_sorted_correctly(tmp_path: Path) -> None:
    """Ensure that files with a mix of extensions are sorted numerially."""
    names = ["data_00001.csv", "data_00001.tar.gz", "data_00002.tar.gz", "data_00003.txt"]
    make_files(tmp_path, names)
    cursor = NumericFileSequenceBatcher(directory=tmp_path, batch_size=10)
    assert [p.name for p in cursor.get_batch()] == names


def test_get_batch_exits_early_with_no_matching_files(
    tmp_path: Path, non_matching_file_names: list[str], caplog: pytest.LogCaptureFixture
) -> None:
    """Ensure that the batcher exits from get_batch early if there are no files that match the file_regex."""
    make_files(tmp_path, non_matching_file_names)
    cursor = NumericFileSequenceBatcher(directory=tmp_path)
    assert cursor.get_batch() == []
    assert len(caplog.records) == 1
    assert caplog.records[-1].levelno == logging.WARNING
    assert caplog.records[-1].message == f"No matching files found in {tmp_path!s}"


# Dynamic / live-directory behaviour
def test_get_batch_picks_up_newly_added_files(tmp_path: Path) -> None:
    """Ensure that adding files to a dir during batching picks up new files correctly."""
    # dir contains log_00001.txt -> log_00003.txt
    make_sequence(tmp_path, "log", "txt", range(1, 4))
    cursor = NumericFileSequenceBatcher(directory=tmp_path, batch_size=10)
    assert cursor.get_batch() == [tmp_path / f"log_{n:05}.txt" for n in [1, 2, 3]]

    (tmp_path / "log_00004.txt").touch()
    # Reset cursor to re-scan from the start
    cursor.start_at = MIN_START_AT
    assert cursor.get_batch() == [tmp_path / f"log_{n:05}.txt" for n in [1, 2, 3, 4]]


def test_get_batch_new_files_within_current_window_are_included(tmp_path: Path) -> None:
    """Ensure that all new files within the current batching params are included, regardless of sequence position."""
    # dir contains log_00001.txt -> log_00005.txt
    make_sequence(tmp_path, "log", "txt", range(1, 6))
    cursor = NumericFileSequenceBatcher(directory=tmp_path, batch_size=3, end_at=13)

    assert cursor.get_batch() == [tmp_path / f"log_{n:05}.txt" for n in [1, 2, 3]]
    assert cursor.start_at == 4  # noqa: PLR2004

    # New files added within the next window before the next call
    (tmp_path / "log_00006.txt").touch()
    (tmp_path / "log_00007.txt").touch()
    assert cursor.get_batch() == [tmp_path / f"log_{n:05}.txt" for n in [4, 5, 6]]
    assert cursor.start_at == 7  # noqa: PLR2004

    # New files added within the next window before the next call
    # Note: we are missing 00008 and 00009
    (tmp_path / "log_00010.txt").touch()
    (tmp_path / "log_00011.txt").touch()
    assert cursor.get_batch() == [tmp_path / f"log_{n:05}.txt" for n in [7, 10, 11]]
    assert cursor.start_at == 12  # noqa: PLR2004

    # add in missing files -- nothing is returned as start_at is at 12
    (tmp_path / "log_00008.txt").touch()
    (tmp_path / "log_00009.txt").touch()
    assert cursor.get_batch() == []
    assert cursor.start_at == 12  # noqa: PLR2004

    # add more files
    (tmp_path / "log_00012.txt").touch()
    (tmp_path / "log_00013.txt").touch()
    assert cursor.get_batch() == [tmp_path / f"log_{n:05}.txt" for n in [12, 13]]
    assert cursor.start_at == 14  # noqa: PLR2004

    # add more files beyond the end_at value
    (tmp_path / "log_00014.txt").touch()
    (tmp_path / "log_00015.txt").touch()
    assert cursor.get_batch() == []
    assert cursor.start_at == 14  # noqa: PLR2004


def test_get_file_batches_pass_yields_all_batches_in_order(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fake_batcher_factory: Any
) -> None:
    """Verify batches are yielded from the batcher in the order it produces them."""
    batch1 = [tmp_path / "0001.xml", tmp_path / "0002.xml"]
    batch2 = [tmp_path / "0003.xml"]
    fake_batcher_cls, _ = fake_batcher_factory([batch1, batch2])
    monkeypatch.setattr(cdm_data_loaders.utils.batcher, "NumericFileSequenceBatcher", fake_batcher_cls)

    settings = MagicMock(input_dir=tmp_path, start_at=None)
    results = list(get_file_batches(settings))

    assert results == [batch1, batch2]


def test_get_file_batches_pass_no_batches_yields_nothing(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fake_batcher_factory: Any
) -> None:
    """Verify a batcher that immediately returns no files produces an empty generator."""
    fake_batcher_cls, _ = fake_batcher_factory([])
    monkeypatch.setattr(cdm_data_loaders.utils.batcher, "NumericFileSequenceBatcher", fake_batcher_cls)

    settings = MagicMock(input_dir=tmp_path, start_at=None)
    results = list(get_file_batches(settings))

    assert results == []


def test_get_file_batches_pass_batcher_constructed_with_input_dir_and_default_batch_size(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fake_batcher_factory: Any
) -> None:
    """Verify the batcher is constructed with the settings' input_dir and the module's default batch size."""
    fake_batcher_cls, construction_calls = fake_batcher_factory([])
    monkeypatch.setattr(cdm_data_loaders.utils.batcher, "NumericFileSequenceBatcher", fake_batcher_cls)

    settings = MagicMock(input_dir=tmp_path, start_at=None)
    list(get_file_batches(settings))

    assert construction_calls[0]["directory"] == tmp_path
    assert construction_calls[0]["batch_size"] == DEFAULT_PIPELINE_BATCH_SIZE


@pytest.mark.parametrize(
    ("start_at", "expect_start_at_key"),
    [
        (None, False),
        ("", False),
        (0, False),  # documents current (possibly unintended) behavior: falsy start_at is dropped
        ("20230101", True),
        (5, True),
    ],
)
def test_get_file_batches_pass_falsy_start_at_omitted_truthy_start_at_forwarded(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    fake_batcher_factory: Any,
    start_at: Any,
    expect_start_at_key: bool,
) -> None:
    """Verify `start_at` is forwarded to the batcher only when truthy, per the `if settings.start_at:` check."""
    fake_batcher_cls, construction_calls = fake_batcher_factory([])
    monkeypatch.setattr(cdm_data_loaders.utils.batcher, "NumericFileSequenceBatcher", fake_batcher_cls)

    settings = MagicMock(input_dir=tmp_path, start_at=start_at)
    list(get_file_batches(settings))

    assert ("start_at" in construction_calls[0]) == expect_start_at_key
    if expect_start_at_key:
        assert construction_calls[0]["start_at"] == start_at


def test_get_file_batches_pass_logs_debug_message_listing_files_per_batch(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, caplog: pytest.LogCaptureFixture, fake_batcher_factory: Any
) -> None:
    """Verify a debug-level log listing the batch's files is emitted once per batch, not once per file."""
    caplog.set_level(logging.DEBUG)
    batch1 = [tmp_path / "0001.xml", tmp_path / "0002.xml"]
    batch2 = [tmp_path / "0003.xml"]
    fake_batcher_cls, _ = fake_batcher_factory([batch1, batch2])
    monkeypatch.setattr(cdm_data_loaders.utils.batcher, "NumericFileSequenceBatcher", fake_batcher_cls)

    settings = MagicMock(input_dir=tmp_path, start_at=None)
    list(get_file_batches(settings))

    batch_logs = [r for r in caplog.records if r.msg == "Files to be processed:%s"]
    assert len(batch_logs) == 2
