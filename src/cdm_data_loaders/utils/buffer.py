"""Buffer classes for accumulating rows per table prior to yielding them to dlt.

dlt loads most efficiently when transformers and resources yield pages of rows
rather than individual items. These buffers accumulate parsed rows per table and
hand back full pages when a table's buffer reaches its size limit.
"""

from collections import defaultdict
from collections.abc import Generator
from typing import Annotated, Any, Self

import dlt
import pyarrow as pa
from dlt.extract.items import DataItemWithMeta
from pydantic import BaseModel, ConfigDict, Field, PositiveInt, PrivateAttr, model_validator

from cdm_data_loaders.core.fields import NonEmptyStr

DEFAULT_BUFFER_MAX_ITEMS = 1000


type Buffers = Annotated[dict[NonEmptyStr, list[Any]], Field(description="Buffer implementation")]

type MaxItems = Annotated[PositiveInt, Field(description="Maximum number of items per buffer")]

type TableNames = Annotated[list[NonEmptyStr], Field(description="Names of the tables to be stored in the buffer")]

type TableToSchemaMap = Annotated[
    dict[NonEmptyStr, pa.Schema], Field(min_length=1, description="Mapping of table names to schemas")
]


class BufferCore(BaseModel):
    """Accumulates rows per table and yields them as pages to dlt."""

    model_config = ConfigDict(extra="forbid")

    _table_names: TableNames = PrivateAttr(default_factory=list)
    _buffers: Buffers = PrivateAttr(default_factory=dict)

    max_items: MaxItems = Field(default=DEFAULT_BUFFER_MAX_ITEMS)

    @model_validator(mode="after")
    def set_up_buffers(self) -> Self:
        """Set up the buffers after successfully instantiating the Buffer instance."""
        self._init_buffers()
        return self

    def _init_buffers(self) -> None:
        """Initialise the buffers."""
        if self._table_names:
            self._buffers: dict[str, list[Any]] = {name: [] for name in self._table_names}
        else:
            self._buffers = defaultdict(list)

    def _add_item(self, table_name: str, item_to_add: Any) -> Generator[DataItemWithMeta, Any]:  # noqa: ANN401
        """Add an item to a table. Yield a page if the table reaches max_items.

        :param table_name: name of the table to add the item to
        :type  table_name: str
        :param item_to_add: the item
        :type  item_to_add: Any
        :yield: generator that yields collected rows to the appropriate dlt table
        :rtype: Generator[DataItemWithMeta, Any]
        """
        self._buffers[table_name].append(item_to_add)
        if len(self._buffers[table_name]) >= self.max_items:
            yield dlt.mark.with_table_name(self._buffers[table_name], table_name)
            self._buffers[table_name] = []

    def _add_items(self, items_to_add: dict[str, list[Any]]) -> Generator[DataItemWithMeta, Any]:
        """Add items to per-table buffers. Yield a page for any table that reaches max_items.

        :param items_to_add: dictionary of lists of items to add to per-table buffers, indexed by table name
        :type  items_to_add: dict[str, list[Any]]
        :yield: generator that yields collected rows to the appropriate dlt table(s)
        :rtype: Generator[DataItemWithMeta, Any]
        """
        for table_name, contents in items_to_add.items():
            if not contents:
                continue
            self._buffers[table_name].extend(contents)
            if len(self._buffers[table_name]) >= self.max_items:
                yield dlt.mark.with_table_name(self._buffers[table_name], table_name)
                self._buffers[table_name] = []

    def flush(self) -> Generator[DataItemWithMeta, Any]:
        """Yield any remaining buffered rows as final pages and empty the buffers."""
        for table_name, rows in self._buffers.items():
            if rows:
                yield dlt.mark.with_table_name(rows, table_name)
                self._buffers[table_name] = []


class DictBuffer(BufferCore):
    """A buffer for parsed results keyed by table name, e.g. from a parse_fn returning dict[str, rows]."""

    def add_item(self, table_name: str, item_to_add: Any) -> Generator[DataItemWithMeta, Any]:  # noqa: ANN401
        """Add an item to a table. Yield a page if the table reaches max_items.

        :param table_name: name of the table to add the item to
        :type  table_name: str
        :param item_to_add: the item
        :type  item_to_add: Any
        :yield: generator that yields collected rows to the appropriate dlt table
        :rtype: Generator[DataItemWithMeta, Any]
        """
        yield from super()._add_item(table_name, item_to_add)

    def add_items(self, items_to_add: dict[str, list[Any]]) -> Generator[DataItemWithMeta, Any]:
        """Add items to per-table buffers. Yield a page for any table that reaches max_items.

        :param items_to_add: dictionary of lists of items to add to per-table buffers, indexed by table name
        :type  items_to_add: dict[str, list[Any]]
        :yield: generator that yields collected rows to the appropriate dlt table(s)
        :rtype: Generator[DataItemWithMeta, Any]
        """
        yield from super()._add_items(items_to_add)


class ListBuffer(BufferCore):
    """A single-table buffer for routing one row at a time, e.g. validation results."""

    table_name: Annotated[NonEmptyStr, Field(..., description="Table to save the data to")]

    @model_validator(mode="after")
    def set_up_table_names(self) -> Self:
        """Set up the buffers after successfully instantiating the Buffer instance."""
        self._init_buffers()
        return self

    def add_item(self, item_to_add: Any) -> Generator[DataItemWithMeta, Any]:  # noqa: ANN401
        """Add an item to the table. Yield a page if the table reaches max_items.

        :param item_to_add: the item
        :type  item_to_add: Any
        :yield: generator that yields collected rows to the appropriate dlt table
        :rtype: Generator[DataItemWithMeta, Any]
        """
        yield from super()._add_item(self.table_name, item_to_add)


class ArrowBufferCore(BaseModel):
    """Accumulates rows per table and yields them as pyarrow tables to dlt."""

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    _buffers: Buffers = PrivateAttr(default_factory=dict)

    max_items: MaxItems = Field(default=DEFAULT_BUFFER_MAX_ITEMS)
    table_to_schema_map: TableToSchemaMap

    @model_validator(mode="after")
    def set_up_buffers(self) -> Self:
        """Set up the buffers after successfully instantiating the Buffer instance."""
        self._init_buffers()
        return self

    def _init_buffers(self) -> None:
        """Initialise the buffers."""
        if not self.table_to_schema_map:
            err_msg = "Must supply table-to-schema mappings to use ArrowBuffers"
            raise ValueError(err_msg)

        self._buffers: dict[str, list[Any]] = {name: [] for name in self.table_to_schema_map}

    def _add_item(self, table_name: str, item_to_add: Any) -> Generator[DataItemWithMeta, Any]:  # noqa: ANN401
        """Add an item to a table. Yield a page if the table reaches max_items.

        :param table_name: name of the table to add the item to
        :type  table_name: str
        :param item_to_add: the item
        :type  item_to_add: Any
        :yield: generator that yields collected rows to the appropriate dlt table
        :rtype: Generator[DataItemWithMeta, Any]
        """
        if table_name not in self._buffers:
            err_msg = f"Cannot add data for {table_name} table: no schema supplied"
            raise ValueError(err_msg)

        self._buffers[table_name].append(item_to_add)
        if len(self._buffers[table_name]) >= self.max_items:
            yield dlt.mark.with_table_name(
                pa.Table.from_pylist(self._buffers[table_name], schema=self.table_to_schema_map[table_name]),
                table_name,
            )
            self._buffers[table_name] = []

    def _add_items(self, items_to_add: dict[str, list[Any]]) -> Generator[DataItemWithMeta, Any]:
        """Add items to per-table buffers. Yield a page for any table that reaches max_items."""
        for table_name, contents in items_to_add.items():
            if not contents:
                continue
            if table_name not in self._buffers:
                err_msg = f"Cannot add data for {table_name} table: no schema supplied"
                raise ValueError(err_msg)

            self._buffers[table_name].extend(contents)
            if len(self._buffers[table_name]) >= self.max_items:
                yield dlt.mark.with_table_name(
                    pa.Table.from_pylist(self._buffers[table_name], schema=self.table_to_schema_map[table_name]),
                    table_name,
                )
                self._buffers[table_name] = []

    def flush(self) -> Generator[DataItemWithMeta, Any]:
        """Yield any remaining buffered rows as final pages and empty the buffers."""
        for table_name, rows in self._buffers.items():
            if rows:
                yield dlt.mark.with_table_name(
                    pa.Table.from_pylist(rows, schema=self.table_to_schema_map[table_name]), table_name
                )
                self._buffers[table_name] = []


class ArrowDictBuffer(ArrowBufferCore):
    """A pyarrow buffer for parsed results keyed by table name. Results are yielded as pyarrow tables."""

    def add_item(self, table_name: str, item_to_add: Any) -> Generator[DataItemWithMeta, Any]:  # noqa: ANN401
        """Add an item to a table. Yield a page if the table reaches max_items.

        :param table_name: name of the table to add the item to
        :type  table_name: str
        :param item_to_add: the item
        :type  item_to_add: Any
        :yield: generator that yields collected rows to the appropriate dlt table
        :rtype: Generator[DataItemWithMeta, Any]
        """
        yield from super()._add_item(table_name, item_to_add)

    def add_items(self, items_to_add: dict[str, list[Any]]) -> Generator[DataItemWithMeta, Any]:
        """Add items to per-table buffers. Yield a page for any table that reaches max_items.

        :param items_to_add: dictionary of lists of items to add to per-table buffers, indexed by table name
        :type  items_to_add: dict[str, list[Any]]
        :yield: generator that yields collected rows to the appropriate dlt table(s)
        :rtype: Generator[DataItemWithMeta, Any]
        """
        yield from super()._add_items(items_to_add)


class ArrowListBuffer(ArrowBufferCore):
    """A single-table pyarrow buffer for routing one row at a time. Results are yielded as pyarrow tables."""

    _table_name: NonEmptyStr = PrivateAttr()
    table_to_schema_map: TableToSchemaMap = Field(min_length=1, max_length=1)

    @model_validator(mode="after")
    def set_up_table_name(self) -> Self:
        """Set up the buffers after successfully instantiating the Buffer instance."""
        if len(self.table_to_schema_map) != 1:
            err_msg = "ArrowListBuffers can only be used for a single table and schema"
            raise ValueError(err_msg)
        self._table_name = next(iter(self.table_to_schema_map.keys()))
        return self

    def add_item(self, item_to_add: Any) -> Generator[DataItemWithMeta, Any]:  # noqa: ANN401
        """Add an item to the table. Yield a page if the table reaches max_items.

        :param item_to_add: the item
        :type  item_to_add: Any
        :yield: generator that yields collected rows to the appropriate dlt table
        :rtype: Generator[DataItemWithMeta, Any]
        """
        yield from super()._add_items({self._table_name: [item_to_add]})
