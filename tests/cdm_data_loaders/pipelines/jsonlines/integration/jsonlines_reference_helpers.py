"""Reconstruction and canonicalization helpers for the dataset_report jsonl reference tests.

The pipeline flattens nested dicts into ``__``-separated columns, turns lists of
scalars into child tables holding a single ``value`` column, coerces naive
timestamp strings into TZ-aware UTC datetimes, and represents some scalar-or-object
fields (e.g. ``gc_percent``) as ``{"v_double": ...}`` when a float is present. The
helpers here rebuild the original nested entries from the pipeline tables and
canonicalize both sides so they can be compared exactly.
"""

from contextlib import suppress
from datetime import UTC, date, datetime
from typing import Any, Final

from tests.cdm_data_loaders.pipelines.helpers import reconstruct_entries

EXPECTED_ENTRY_COUNT: Final[int] = 52

NAIVE_TIMESTAMP_KEYS: Final[frozenset[str]] = frozenset(
    {"when", "submission_date", "last_updated", "publication_date", "release_date_first_public"}
)

DATE_KEYS: Final[frozenset[str]] = frozenset({"release_date"})

DOUBLE_KEYS: Final[frozenset[str]] = frozenset({"gc_percent"})


def reconstructed_entries(
    tables: dict[str, list[dict[str, Any]]], root_table: str | None = None
) -> list[dict[str, Any]]:
    """Rebuild nested entries from the pipeline's flattened output tables.

    The main table is named after the dataset (e.g. ``dataset``); child tables are
    prefixed with it. The xmltodict reconstruction walk expects a top-level
    ``entry`` table, so the shared root prefix is remapped before delegating.
    """
    tables = {
        name: rows
        for name, rows in tables.items()
        if not name.startswith("_dlt")
        and (root_table is None or name == root_table or name.startswith(f"{root_table}__"))
    }
    if not tables:
        return []
    top_tables = [name for name in tables if "__" not in name]
    if len(top_tables) != 1:
        err_msg = f"Expected exactly one top-level table, found: {sorted(top_tables)}"
        raise ValueError(err_msg)
    root = top_tables[0]
    renamed = {("entry" if name == root else f"entry{name.removeprefix(root)}"): rows for name, rows in tables.items()}
    return reconstruct_entries(renamed)


def canonicalize_reference_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Canonicalize an entry so pipeline output and source JSON compare equal.

    Reverses the transformations dlt's normalization applies, and normalizes the
    two storage representations (parquet datetimes vs jsonl timestamp strings):

    - drops keys whose value is None (dlt stores missing values as NULL columns
      that reconstruct to None);
    - drops empty dicts (e.g. ``"owner": {}``), which dlt stores as no columns;
    - wraps scalar list items in ``{"value": ...}`` (dlt stores lists of scalars
      in a child table with a single ``value`` column);
    - unwraps ``{"v_double": ...}`` on both sides (dlt's representation of a
      float scalar stored alongside a non-float variant of the same field;
    - converts naive ISO timestamp strings for known keys into TZ-aware UTC
      datetimes, matching dlt's timestamp coercion;
    - converts TZ-aware ISO timestamp strings for known keys into datetimes so
      jsonl output (which stores coerced timestamps as strings) matches.
    """
    return _canonicalize_value(entry)


def canonicalize_pipeline_entry(entry: dict[str, Any]) -> dict[str, Any]:
    """Canonicalize a reconstructed pipeline entry for comparison.

    Preserved JSON lists and reconstructed child rows share the same scalar-list
    representation. Null fields and empty objects disappear during flattening.
    Timestamps are normalized for both JSONL strings and Parquet datetimes.
    """
    return _canonicalize_value(entry)


def _canonicalize_value(value: Any, key: str = "") -> Any:
    """Normalize only known storage differences, retaining list order and duplicates."""
    if isinstance(value, dict):
        if key in DOUBLE_KEYS and set(value) == {"v_double"}:
            return _canonicalize_value(value["v_double"], key)
        canonical: dict[str, Any] = {}
        for child_key, child_value in value.items():
            if child_value is None:
                continue
            canonical_child = _canonicalize_value(child_value, child_key)
            if canonical_child == {}:
                continue
            canonical[child_key] = canonical_child
        return canonical
    if isinstance(value, list):
        return [
            (
                {"value": _canonicalize_value(item, key)}
                if not isinstance(item, (dict, list))
                else _canonicalize_value(item, key)
            )
            for item in value
        ]
    if key in NAIVE_TIMESTAMP_KEYS and isinstance(value, str):
        value = _parse_timestamp(value)
    if key in DATE_KEYS and isinstance(value, str):
        with suppress(ValueError):
            value = date.fromisoformat(value)
    if key in DOUBLE_KEYS and isinstance(value, int) and not isinstance(value, bool):
        value = float(value)
    return value


def _parse_timestamp(value: str) -> Any:
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        return value
    if parsed.tzinfo is None:
        return parsed.replace(tzinfo=UTC)
    return parsed


def reference_entries_from_duckdb(connection: Any) -> list[dict[str, Any]]:
    """Extract canonical nested entries from the DuckDB reference dataset.

    Restore inferred dates to JSON strings and remove null STRUCT members added
    for absent fields, so model defaults behave as they do for the source JSONL.
    """
    rows = connection.execute("SELECT * FROM reference").fetchall()
    columns = [description[0] for description in connection.description]
    entries = [_restore_reference_json(dict(zip(columns, row, strict=True))) for row in rows]
    return sorted(entries, key=lambda entry: entry.get("accession") or "")


def _restore_reference_json(value: Any) -> Any:
    """Undo DuckDB's date inference and padding of absent STRUCT fields."""
    if isinstance(value, dict):
        return {key: _restore_reference_json(item) for key, item in value.items() if item is not None}
    if isinstance(value, list):
        return [_restore_reference_json(item) for item in value]
    if isinstance(value, (date, datetime)):
        return value.isoformat()
    return value
