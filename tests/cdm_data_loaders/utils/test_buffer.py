"""Tests for the page buffers in cdm_data_loaders.utils.buffer."""

import re
from collections.abc import Iterator
from typing import Any, Final

import dlt
import pyarrow as pa
import pytest
from pydantic import ValidationError

from cdm_data_loaders.utils.buffer import (
    DEFAULT_BUFFER_MAX_ITEMS,
    ArrowDictBuffer,
    ArrowListBuffer,
    BufferCore,
    DictBuffer,
    ListBuffer,
)

PEOPLE = "people"
EMAILS = "emails"
PEOPLE_SCHEMA: Final[pa.Schema] = pa.schema([pa.field("name", pa.string())])
EMAILS_SCHEMA: Final[pa.Schema] = pa.schema([pa.field("n", pa.int64())])

SCHEMA_MAP: Final[dict[str, pa.Schema]] = {
    PEOPLE: pa.schema([pa.field("name", pa.string())]),
    EMAILS: pa.schema([pa.field("n", pa.int64())]),
    "widget": pa.schema([pa.field("n", pa.int64())]),
}


def name_row(n: int) -> dict[str, Any]:
    """Build a row dict matching the people schema."""
    return {"name": str(n)}


def number_row(n: int) -> dict[str, Any]:
    """Build a row dict matching the emails schema."""
    return {"n": n}


def tagged(pages: list[Any]) -> list[tuple[str, list[Any]]]:
    """Unwrap yielded pages into (table_name, rows) pairs."""
    for page in pages:
        if isinstance(page.data, pa.Table):
            assert page.data.schema == SCHEMA_MAP[page.meta.table_name]

    return [
        (page.meta.table_name, page.data.to_pylist() if isinstance(page.data, pa.Table) else page.data)
        for page in pages
    ]


def arrow_tagged(pages: list[Any]) -> list[tuple[str, pa.Table]]:
    """Unwrap yielded pages into (table_name, pa.Table) pairs."""
    return [(page.meta.table_name, page.data) for page in pages]


# object instantiation and errors


@pytest.mark.parametrize("cls", [DictBuffer, ListBuffer, ArrowDictBuffer, ArrowListBuffer])
def test_buffer_classes_pass_default_max_items(
    cls: type[DictBuffer] | type[ListBuffer] | type[ArrowListBuffer] | type[ArrowDictBuffer],
) -> None:
    """max_items defaults to DEFAULT_BUFFER_MAX_ITEMS."""
    args = {}
    if cls == ListBuffer:
        args = {"table_name": "widget"}
    elif cls in (ArrowDictBuffer, ArrowListBuffer):
        args = {"table_to_schema_map": {PEOPLE: PEOPLE_SCHEMA}}

    buffer = cls(**args)  # pyright: ignore[reportArgumentType]

    assert buffer.max_items == DEFAULT_BUFFER_MAX_ITEMS


def test_dict_buffer_requires_no_args() -> None:
    """Test that the DictBuffer class has no required arguments."""
    buffer = DictBuffer()  # pyright: ignore[reportCallIssue]
    assert isinstance(buffer, BufferCore)


@pytest.mark.parametrize("cls", [ListBuffer, ArrowDictBuffer, ArrowListBuffer])
def test_buffer_classes_required_args(
    cls: type[ListBuffer] | type[ArrowListBuffer] | type[ArrowDictBuffer],
) -> None:
    """Ensure that the required args for each class are flagged in error messages."""
    with pytest.raises(ValidationError, match="Field required") as e:
        cls()  # pyright: ignore[reportCallIssue]

    assert "1 validation error for " in str(e)

    if cls == ListBuffer:
        # table_name
        pattern = re.compile("table_name\n\\s+Field required")
        assert re.search(pattern, str(e))

    if cls in (ArrowListBuffer, ArrowDictBuffer):
        pattern = re.compile("table_to_schema_map\n\\s+Field required")
        assert re.search(pattern, str(e))


@pytest.mark.parametrize("cls", [DictBuffer, ListBuffer, ArrowDictBuffer, ArrowListBuffer])
@pytest.mark.parametrize("max_items", [0, -1], ids=["zero", "negative"])
def test_buffer_classes_fail_non_positive_max_items_raises_error(
    max_items: int,
    cls: type[DictBuffer] | type[ListBuffer] | type[ArrowListBuffer] | type[ArrowDictBuffer],
) -> None:
    """A non-positive max_items raises a pydantic ValidationError.

    Note: some classes will also raise other validation errors, but all will throw the max_items error.
    """
    with pytest.raises(ValidationError, match="Input should be greater than 0"):
        cls(max_items=max_items)  # pyright: ignore[reportCallIssue]


@pytest.mark.parametrize("cls", [DictBuffer, ListBuffer, ArrowDictBuffer, ArrowListBuffer])
def test_buffer_classes_fail_unrecognised_args_throw_error(
    cls: type[DictBuffer] | type[ListBuffer] | type[ArrowListBuffer] | type[ArrowDictBuffer],
) -> None:
    """Unrecognised arguments throw a big fat error."""
    with pytest.raises(ValidationError, match="Extra inputs are not permitted"):
        cls(_buffers="some buffery thing")  # pyright: ignore[reportCallIssue]


@pytest.mark.parametrize(
    ("args", "err_msg"),
    [
        ({}, "Dictionary should have at least 1 item after validation"),
        ({"": PEOPLE_SCHEMA}, "String should have at least 1 character"),
        ({"foo": "bar"}, "Input should be an instance of Schema"),
        ({"pop": PEOPLE_SCHEMA, "pip": {"this": "that"}}, "Input should be an instance of Schema"),
    ],
)
@pytest.mark.parametrize("cls", [ArrowListBuffer, ArrowDictBuffer])
def test_arrow_buffer_classes_fail_invalid_table_to_schema_map(
    cls: type[ArrowListBuffer] | type[ArrowDictBuffer], args: Any, err_msg: str
) -> None:
    """Arrow buffer classes throw an error on invalid input args."""
    with pytest.raises(ValidationError, match=err_msg):
        cls(table_to_schema_map=args)  # pyright: ignore[reportCallIssue]


@pytest.mark.parametrize("cls", [ArrowListBuffer, ArrowDictBuffer])
def test_arrow_buffers_accept_empty_schemas(cls: type[ArrowListBuffer] | type[ArrowDictBuffer]) -> None:
    """Arrow buffers accept an empty schema."""
    buffer = cls(table_to_schema_map={"pip": pa.schema([])})  # pyright: ignore[reportCallIssue]
    assert isinstance(buffer, cls)


def test_arrow_list_buffer_fail_too_many_mappings() -> None:
    """Only one table to schema mapping allowed in an ArrowListBuffer class."""
    with pytest.raises(ValidationError, match="Dictionary should have at most 1 item after validation"):
        ArrowListBuffer(table_to_schema_map={"pop": PEOPLE_SCHEMA, "pip": PEOPLE_SCHEMA})  # pyright: ignore[reportCallIssue]


# buffer functions
@pytest.mark.parametrize("cls", [DictBuffer, ArrowDictBuffer])
def test_dict_buffer_add_item_pass_yields_page_at_max_items(cls: DictBuffer | ArrowDictBuffer) -> None:
    """A single item addition yields a page when max_items is reached."""
    args = {}
    if cls == ArrowDictBuffer:
        args["table_to_schema_map"] = SCHEMA_MAP
    buffer = cls(max_items=3, **args)  # pyright: ignore[reportCallIssue]

    # Add 2 items - no yield
    assert list(buffer.add_item(PEOPLE, name_row(1))) == []
    assert list(buffer.add_item(PEOPLE, name_row(2))) == []

    # Add 3rd item - yields page
    pages = list(buffer.add_item(PEOPLE, name_row(3)))
    assert len(pages) == 1
    assert tagged(pages) == [(PEOPLE, [name_row(1), name_row(2), name_row(3)])]

    # Buffer should be empty
    assert list(buffer.flush()) == []


@pytest.mark.parametrize("cls", [DictBuffer, ArrowDictBuffer])
def test_dict_buffer_core_item_pass_routes_to_different_tables(cls: DictBuffer | ArrowDictBuffer) -> None:
    """add_item correctly separates rows into different table buffers."""
    args = {}
    if cls == ArrowDictBuffer:
        args["table_to_schema_map"] = SCHEMA_MAP
    buffer = cls(max_items=2, **args)  # pyright: ignore[reportCallIssue]

    assert list(buffer.add_item(PEOPLE, name_row(1))) == []
    assert list(buffer.add_item(EMAILS, number_row(1))) == []

    # Fill people buffer
    pages_people = list(buffer.add_item(PEOPLE, name_row(2)))
    assert tagged(pages_people) == [(PEOPLE, [name_row(1), name_row(2)])]

    # Fill emails buffer
    pages_emails = list(buffer.add_item(EMAILS, number_row(2)))
    assert tagged(pages_emails) == [(EMAILS, [number_row(1), number_row(2)])]


def test_list_buffer_add_item_pass_delegates_to_core() -> None:
    """ListBuffer.add_item uses the core logic with its fixed table_name."""
    buffer = ListBuffer(table_name=PEOPLE, max_items=2)

    assert list(buffer.add_item(name_row(1))) == []
    pages = list(buffer.add_item(name_row(2)))
    assert tagged(pages) == [(PEOPLE, [name_row(1), name_row(2)])]


def test_arrow_dict_buffer_pass_flush_yields_remaining_rows() -> None:
    """flush() yields the rows still held per table as pa.Tables and empties the buffers."""
    buffer = ArrowDictBuffer(max_items=10, table_to_schema_map=SCHEMA_MAP)
    list(buffer.add_items({PEOPLE: [name_row(1)], EMAILS: [number_row(2)]}))

    people_page, email_page = arrow_tagged(list(buffer.flush()))
    assert people_page == (PEOPLE, pa.Table.from_pylist([name_row(1)], schema=SCHEMA_MAP[PEOPLE]))
    assert email_page == (EMAILS, pa.Table.from_pylist([number_row(2)], schema=SCHEMA_MAP[EMAILS]))
    assert list(buffer.flush()) == []


def test_arrow_dict_buffer_pass_empty_table_rows_are_ignored() -> None:
    """Empty row lists contribute nothing and yield no pages."""
    buffer = ArrowDictBuffer(max_items=1, table_to_schema_map=SCHEMA_MAP)

    assert list(buffer.add_items({PEOPLE: []})) == []
    assert list(buffer.flush()) == []


@pytest.mark.parametrize("cls", [DictBuffer, ArrowDictBuffer])
def test_dict_buffer_pass_rows_accumulate_until_max_items(cls: type[DictBuffer] | type[ArrowDictBuffer]) -> None:
    """No pages are yielded while every table's buffer is below max_items."""
    args: dict[str, Any] = {"max_items": 5}
    if cls == ArrowDictBuffer:
        args["table_to_schema_map"] = SCHEMA_MAP
    buffer = cls(**args)  # pyright: ignore[reportCallIssue]

    assert list(buffer.add_items({PEOPLE: [name_row(1), name_row(2)], EMAILS: [number_row(3)]})) == []


@pytest.mark.parametrize("cls", [DictBuffer, ArrowDictBuffer])
def test_dict_buffer_pass_full_table_yields_page_and_empties(cls: type[DictBuffer] | type[ArrowDictBuffer]) -> None:
    """A table reaching max_items yields its rows as one page and the buffer resets."""
    args: dict[str, Any] = {}
    if cls == ArrowDictBuffer:
        args["table_to_schema_map"] = SCHEMA_MAP
    buffer = cls(max_items=2, **args)  # pyright: ignore[reportCallIssue]

    first = tagged(list(buffer.add_items({PEOPLE: [name_row(1), name_row(2)]})))
    assert first == [(PEOPLE, [name_row(1), name_row(2)])]

    assert list(buffer.add_items({PEOPLE: [name_row(3)]})) == []

    second = tagged(list(buffer.add_items({PEOPLE: [name_row(4)]})))
    assert second == [(PEOPLE, [name_row(3), name_row(4)])]
    assert list(buffer.add_items({PEOPLE: [name_row(5)]})) == []


@pytest.mark.parametrize("cls", [DictBuffer, ArrowDictBuffer])
def test_dict_buffer_pass_flush_yields_remaining_rows(cls: type[DictBuffer] | type[ArrowDictBuffer]) -> None:
    """flush() yields the rows still held per table and empties the buffers."""
    args = {}
    if cls == ArrowDictBuffer:
        args["table_to_schema_map"] = SCHEMA_MAP
    buffer = cls(max_items=10, **args)  # type: ignore[reportCallIssue]
    list(buffer.add_items({PEOPLE: [name_row(1)], EMAILS: [number_row(2)]}))

    assert tagged(list(buffer.flush())) == [(PEOPLE, [name_row(1)]), (EMAILS, [number_row(2)])]
    assert list(buffer.flush()) == []


@pytest.mark.parametrize("cls", [DictBuffer, ArrowDictBuffer])
def test_dict_buffer_pass_empty_table_rows_are_ignored(cls: type[DictBuffer] | type[ArrowDictBuffer]) -> None:
    """Empty row lists contribute nothing and yield no pages."""
    args = {}
    if cls == ArrowDictBuffer:
        args["table_to_schema_map"] = SCHEMA_MAP
    buffer = cls(max_items=1, **args)  # type: ignore[reportCallIssue]

    assert list(buffer.add_items({PEOPLE: []})) == []
    assert list(buffer.flush()) == []


def test_dict_buffer_pass_unknown_table_is_created_on_demand() -> None:
    """A table name not seen at construction is buffered like any other."""
    buffer = DictBuffer()  # type: ignore[reportCallIssue]

    assert list(buffer.add_items({"surprise": [{"wtf?": "boo!"}]})) == []
    assert tagged(list(buffer.flush())) == [("surprise", [{"wtf?": "boo!"}])]


def test_arrow_dict_buffer_fail_unknown_table_throws_error() -> None:
    """A table name not seen at construction triggers an error in an ArrowDictBuffer."""
    buffer = ArrowDictBuffer(table_to_schema_map=SCHEMA_MAP)  # pyright: ignore[reportCallIssue]

    with pytest.raises(ValueError, match="Cannot add data for surprise table: no schema supplied"):
        # n.b. error does not show up until you try to use the generator
        list(buffer.add_items({"surprise": [{"wtf?": "boo!"}]}))


@pytest.mark.parametrize("cls", [DictBuffer, ArrowDictBuffer])
def test_dict_buffer_pass_multiple_tables_page_independently(cls: type[DictBuffer] | type[ArrowDictBuffer]) -> None:
    """Each table pages on its own count; a full table does not flush other tables."""
    args = {}
    if cls == ArrowDictBuffer:
        args["table_to_schema_map"] = SCHEMA_MAP
    buffer = cls(max_items=2, **args)  # type: ignore[reportCallIssue]
    first = tagged(list(buffer.add_items({PEOPLE: [name_row(1), name_row(2)], EMAILS: [number_row(3)]})))
    assert first == [(PEOPLE, [name_row(1), name_row(2)])]

    second = tagged(list(buffer.add_items({PEOPLE: [name_row(4)], EMAILS: [number_row(5)]})))
    assert second == [(EMAILS, [number_row(3), number_row(5)])]

    assert tagged(list(buffer.flush())) == [(PEOPLE, [name_row(4)])]
    assert list(buffer.flush()) == []


@pytest.mark.parametrize("cls", [ListBuffer, ArrowListBuffer])
def test_list_buffer_pass_single_item_adds_and_pages(cls: type[ListBuffer] | type[ArrowListBuffer]) -> None:
    """add_item routes to the buffer's table and pages when max_items is reached."""
    args = {"table_name": "widget"} if cls == ListBuffer else {"table_to_schema_map": {"widget": SCHEMA_MAP["widget"]}}
    buffer = cls(**args, max_items=3)  # type: ignore[reportCallIssue]
    assert list(buffer.add_item(number_row(1))) == []
    assert list(buffer.add_item(number_row(2))) == []

    pages = tagged(list(buffer.add_item(number_row(3))))

    assert pages == [("widget", [number_row(1), number_row(2), number_row(3)])]
    assert list(buffer.add_item(number_row(4))) == []


@pytest.mark.parametrize("cls", [ListBuffer, ArrowListBuffer])
def test_list_buffer_pass_flush_yields_remaining_rows(cls: type[ListBuffer] | type[ArrowListBuffer]) -> None:
    """flush() yields rows still held for the buffer's table."""
    args = {"table_name": "widget"} if cls == ListBuffer else {"table_to_schema_map": {"widget": SCHEMA_MAP["widget"]}}
    buffer = cls(**args, max_items=10)  # type: ignore[reportCallIssue]
    list(buffer.add_item(number_row(1)))

    assert tagged(list(buffer.flush())) == [("widget", [number_row(1)])]
    assert list(buffer.flush()) == []


@pytest.mark.parametrize("cls", [ListBuffer, ArrowListBuffer])
def test_arrow_list_buffer_pass_buffer_scope_is_call_local(cls: type[ListBuffer] | type[ArrowListBuffer]) -> None:
    """A list buffer created inside one transformer call yields every row added to it."""
    max_items = 3
    args = {"table_name": PEOPLE} if cls == ListBuffer else {"table_to_schema_map": {PEOPLE: PEOPLE_SCHEMA}}
    buffer = cls(**args, max_items=max_items)  # type: ignore[reportCallIssue]
    emitted: list[Any] = []
    for n in range(10):
        emitted.extend(buffer.add_item(name_row(n)))
    emitted.extend(buffer.flush())

    values = [page.data for page in emitted]
    if cls == ArrowListBuffer:
        assert [value for table in values for value in table.column("name").to_pylist()] == [str(n) for n in range(10)]
        assert [table.num_rows for table in values] == [3, 3, 3, 1]
    else:
        assert [value for table in values for value in table] == [{"name": str(n)} for n in range(10)]
        assert [len(row) for row in values] == [3, 3, 3, 1]


class _NaiveSharedBuffer:
    """Recreates the unsafe pattern rejected during research.

    One mutable buffer instance is shared across every transformer call, with the
    check-then-append-then-replace sequence split across call boundaries. When dlt
    runs a parallelized transformer, calls execute concurrently in worker threads:
    two calls interleave their read of `self.rows` and the later append replaces the
    list captured by the earlier call, silently dropping that call's row.
    """

    def __init__(self, max_items: int) -> None:
        self.rows: list[dict[str, int]] = []
        self.max_items = max_items

    def add(self, item: dict[str, int]) -> list[dict[str, int]]:
        if len(self.rows) >= self.max_items:
            page = self.rows
            self.rows = []
            return page
        self.rows.append(item)
        return []


def _run_naive_shared_buffer_transformer(total: int, max_items: int) -> list[int]:
    """Run a parallelized transformer over one shared _NaiveSharedBuffer. Return emitted values."""
    buffer = _NaiveSharedBuffer(max_items)

    @dlt.transformer(parallelized=True)
    def naive(item: int) -> Iterator[int]:
        return iter(row["v"] for row in buffer.add({"v": item}))

    @dlt.resource
    def source() -> Iterator[int]:
        yield from range(total)

    return [o["v"] for o in (source() | naive) if isinstance(o, dict)]


def _run_page_per_call_transformer(total: int, max_items: int, chunk_size: int) -> list[int]:
    """Run the shipped page-per-call pattern under a parallelized transformer. Return emitted values."""

    @dlt.transformer(parallelized=True)
    def safe(page: list[int]) -> Iterator[dict[str, int]]:
        buffer = ListBuffer(table_name="safe_table", max_items=max_items)  # pyright: ignore[reportCallIssue]
        for item in page:
            yield from buffer.add_item({"v": item})  # pyright: ignore[reportReturnType]
        yield from buffer.flush()  # pyright: ignore[reportReturnType]

    @dlt.resource
    def source() -> Iterator[list[int]]:
        for start in range(0, total, chunk_size):
            yield list(range(start, min(start + chunk_size, total)))

    return [o["v"] for o in (source() | safe) if isinstance(o, dict)]


def test_naive_shared_buffer_fail_loses_rows_under_parallelized_transformer() -> None:
    """Recreates the research race: a buffer shared across parallelized transformer calls loses rows.

    dlt submits each transformer call to a thread pool; the calls interleave their
    check-then-append-then-replace sequences on the shared buffer, so rows vanish.
    """
    total = 50
    values = _run_naive_shared_buffer_transformer(total=total, max_items=5)

    assert len(values) < total, "race did not trigger; increase total or retry"
    assert sorted(values) != list(range(total))


def test_list_buffer_pass_loses_nothing_under_parallelized_transformer() -> None:
    """The shipped pattern buffers inside one transformer call, so no rows are lost.

    The buffer is created per call and dies with the call: no state is shared across
    the worker threads dlt uses for a parallelized transformer.
    """
    total = 50
    values = _run_page_per_call_transformer(total=total, max_items=5, chunk_size=7)

    assert sorted(values) == list(range(total))


def test_list_buffer_pass_buffer_scope_is_call_local() -> None:
    """A ListBuffer created inside one transformer call yields every row added to it.

    Demonstrates the invariant that makes the pattern safe under parallelism: the
    buffer never escapes the call that created it, so concurrent calls cannot
    interleave reads and writes.
    """
    max_items = 3
    buffer = ListBuffer(table_name="widget", max_items=max_items)  # pyright: ignore[reportCallIssue]
    emitted: list[Any] = []
    for n in range(10):
        emitted.extend(buffer.add_item(name_row(n)))
    emitted.extend(buffer.flush())

    values = [page.data for page in emitted]
    assert [v for page in values for v in page] == [name_row(n) for n in range(10)]
    assert [len(page) for page in values] == [3, 3, 3, 1]
