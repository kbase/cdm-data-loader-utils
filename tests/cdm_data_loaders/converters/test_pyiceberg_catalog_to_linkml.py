"""End-to-end tests for pyiceberg_to_linkml.

Chains a real pyiceberg SQL catalog (file-backed SQLite) through
dump_catalog_schemas and validates the emitted documents.
"""

import json
import logging
from collections.abc import Generator
from datetime import UTC, datetime
from pathlib import Path

import pytest
from linkml.linter.linter import Linter
from linkml_runtime.linkml_model.meta import ClassDefinition, SchemaDefinition, SlotDefinition
from linkml_runtime.loaders import yaml_loader
from pyiceberg import catalog as pyiceberg_catalog
from pyiceberg.catalog import Catalog, load_catalog
from pyiceberg.schema import Schema
from pyiceberg.types import (
    DateType,
    DecimalType,
    DoubleType,
    IntegerType,
    ListType,
    LongType,
    MapType,
    NestedField,
    StringType,
    StructType,
    TimestamptzType,
)
from pyiceberg.utils.config import Config

from cdm_data_loaders.converters.pyiceberg_catalog_to_linkml import (
    IcebergToLinkMLSettings,
    dump_catalog_schemas,
)

CATALOG_NAME = "testcatlinkml"


@pytest.fixture
def catalog_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Point the testcat catalog at a file-backed SQLite DB under tmp_path."""
    db_path = tmp_path / "catalog_linkml.db"
    warehouse = tmp_path / "warehouse_linkml"
    warehouse.mkdir()
    monkeypatch.setenv(f"PYICEBERG_CATALOG__{CATALOG_NAME.upper()}__TYPE", "sql")
    monkeypatch.setenv(f"PYICEBERG_CATALOG__{CATALOG_NAME.upper()}__URI", f"sqlite:///{db_path}")
    monkeypatch.setenv(f"PYICEBERG_CATALOG__{CATALOG_NAME.upper()}__WAREHOUSE", f"file://{warehouse}")
    monkeypatch.setattr(pyiceberg_catalog, "_ENV_CONFIG", Config())
    return tmp_path


@pytest.fixture
def catalog(catalog_env: Path) -> Generator[Catalog]:  # noqa: ARG001
    """Load the test catalog, closing its SQLite-backed engine on teardown."""
    with load_catalog(CATALOG_NAME) as loaded_catalog:
        yield loaded_catalog


@pytest.fixture
def wide_schema() -> Schema:
    """A schema exercising scalars, nested structs, lists, and maps."""
    return Schema(
        NestedField(1, "id", LongType(), required=True, doc="primary key"),
        NestedField(2, "name", StringType(), required=True),
        NestedField(3, "count", IntegerType(), required=False),
        NestedField(4, "ratio", DoubleType(), required=False),
        NestedField(5, "as_of", TimestamptzType(), required=False),
        NestedField(6, "born_on", DateType(), required=False),
        NestedField(7, "rate", DecimalType(12, 4), required=False),
        NestedField(8, "tags", ListType(element_id=9, element=StringType(), element_required=False), required=False),
        NestedField(
            10,
            "attrs",
            MapType(key_id=11, key_type=StringType(), value_id=12, value_type=DoubleType(), value_required=False),
            required=False,
        ),
        NestedField(
            13,
            "address",
            StructType(
                NestedField(14, "street", StringType(), required=True),
                NestedField(15, "city", StringType(), required=False),
            ),
            required=False,
        ),
    )


@pytest.mark.parametrize("group_by_namespace", [False, True], ids=["per-table", "per-namespace"])
def test_dump_catalog_schemas_pass_yaml_models(
    catalog: Catalog, catalog_env: Path, simple_schema: Schema, wide_schema: Schema, group_by_namespace: bool
) -> None:
    """Write YAML schemas with valid LinkML models and a shared timezone-aware timestamp."""
    catalog.create_namespace("ns")
    catalog.create_table(("ns", "simple"), schema=simple_schema)
    catalog.create_table(("ns", "wide"), schema=wide_schema)

    out_dir = catalog_env / "schemas"
    settings = IcebergToLinkMLSettings(  # pyright: ignore[reportCallIssue]
        catalog=CATALOG_NAME, output_dir=str(out_dir), group_by_namespace=group_by_namespace
    )
    dump_catalog_schemas(settings)

    files = sorted(out_dir.iterdir())
    assert [f.name for f in files] == (
        ["ns.schema.yaml"] if group_by_namespace else ["ns.simple.schema.yaml", "ns.wide.schema.yaml"]
    )

    simple_doc = yaml_loader.load(str(files[0]), target_class=SchemaDefinition)
    assert simple_doc.name == ("ns" if group_by_namespace else "simple")
    assert simple_doc.id == ("urn:iceberg:ns" if group_by_namespace else "urn:iceberg:ns.simple")
    assert {"simple": simple_doc.classes["simple"]} == {
        "simple": ClassDefinition(
            name="simple",
            attributes={
                "id": SlotDefinition(name="id", range="integer", required=True, description="primary key"),
                "name": SlotDefinition(name="name", range="string", description="display name"),
            },
        )
    }

    wide_doc = simple_doc if group_by_namespace else yaml_loader.load(str(files[1]), target_class=SchemaDefinition)
    assert set(wide_doc.classes) == (
        {"simple", "wide", "wide_attrs_map", "wide_address"}
        if group_by_namespace
        else {"wide", "wide_attrs_map", "wide_address"}
    )
    assert wide_doc.classes["wide"].attributes == {
        "id": SlotDefinition(name="id", range="integer", required=True, description="primary key"),
        "name": SlotDefinition(name="name", range="string", required=True),
        "count": SlotDefinition(name="count", range="integer"),
        "ratio": SlotDefinition(name="ratio", range="float"),
        "as_of": SlotDefinition(name="as_of", range="datetime"),
        "born_on": SlotDefinition(name="born_on", range="date"),
        "rate": SlotDefinition(name="rate", range="decimal"),
        "tags": SlotDefinition(name="tags", range="string", multivalued=True),
        "attrs": SlotDefinition(name="attrs", range="wide_attrs_map", multivalued=True),
        "address": SlotDefinition(name="address", range="wide_address"),
    }
    assert wide_doc.classes["wide_attrs_map"].attributes == {
        "key": SlotDefinition(name="key", range="string", required=True),
        "value": SlotDefinition(name="value", range="float"),
    }
    assert wide_doc.classes["wide_address"].attributes == {
        "street": SlotDefinition(name="street", range="string", required=True),
        "city": SlotDefinition(name="city", range="string"),
    }
    generated_at = simple_doc.annotations["generated_at"].value
    assert wide_doc.annotations["generated_at"].value == generated_at
    assert datetime.fromisoformat(generated_at).tzinfo == UTC

    for schema_path in files:
        with pytest.raises(json.JSONDecodeError):
            json.loads(schema_path.read_text(encoding="utf-8"))
        assert list(Linter.validate_schema(str(schema_path))) == []


@pytest.mark.parametrize("group_by_namespace", [False, True], ids=["per-table", "per-namespace"])
@pytest.mark.parametrize("empty_namespace", [False, True], ids=["empty-catalog", "empty-namespace"])
def test_dump_catalog_schemas_pass_empty_catalog(
    catalog: Catalog, catalog_env: Path, group_by_namespace: bool, empty_namespace: bool
) -> None:
    """Empty catalogs and namespaces write no files."""
    if empty_namespace:
        catalog.create_namespace("empty")
    out_dir = catalog_env / "schemas"
    settings = IcebergToLinkMLSettings(  # pyright: ignore[reportCallIssue]
        catalog=CATALOG_NAME, output_dir=str(out_dir), group_by_namespace=group_by_namespace
    )
    dump_catalog_schemas(settings)
    assert list(out_dir.iterdir()) == []


@pytest.mark.parametrize("group_by_namespace", [False, True], ids=["per-table", "per-namespace"])
@pytest.mark.parametrize("include_good", [False, True], ids=["all-fail", "partial-failure"])
def test_dump_catalog_schemas_pass_skips_failing_table_and_logs(  # noqa: PLR0917
    catalog: Catalog,
    catalog_env: Path,
    simple_schema: Schema,
    caplog: pytest.LogCaptureFixture,
    group_by_namespace: bool,
    include_good: bool,
) -> None:
    """A table that fails to convert is logged and skipped; other tables still succeed."""
    catalog.create_namespace("ns")
    bad_table = catalog.create_table(("ns", "bad"), schema=simple_schema)
    bad_table.io.delete(bad_table.metadata_location)
    if include_good:
        catalog.create_table(("ns", "good"), schema=simple_schema)

    out_dir = catalog_env / "schemas"
    settings = IcebergToLinkMLSettings(  # pyright: ignore[reportCallIssue]
        catalog=CATALOG_NAME, output_dir=str(out_dir), group_by_namespace=group_by_namespace
    )
    with caplog.at_level(logging.WARNING):
        dump_catalog_schemas(settings)

    expected_files = ["ns.schema.yaml" if group_by_namespace else "ns.good.schema.yaml"] if include_good else []
    assert [f.name for f in sorted(out_dir.iterdir())] == expected_files
    if include_good:
        written = yaml_loader.load(str(out_dir / expected_files[0]), target_class=SchemaDefinition)
        assert set(written.classes) == {"good"}
    assert "Failed to convert table ('ns', 'bad'); skipping" in caplog.text
    assert "Completed with 1 table(s) skipped due to conversion failures" in caplog.text


@pytest.mark.parametrize("group_by_namespace", [False, True], ids=["per-table", "per-namespace"])
def test_dump_catalog_schemas_pass_overwrites_existing_file_with_warning(
    catalog: Catalog,
    catalog_env: Path,
    simple_schema: Schema,
    caplog: pytest.LogCaptureFixture,
    group_by_namespace: bool,
) -> None:
    """A pre-existing output file for a table is overwritten, with a warning logged first."""
    catalog.create_namespace("ns")
    catalog.create_table(("ns", "simple"), schema=simple_schema)

    out_dir = catalog_env / "schemas"
    out_dir.mkdir(parents=True)
    stale_file = out_dir / ("ns.schema.yaml" if group_by_namespace else "ns.simple.schema.yaml")
    stale_file.write_text("stale")

    settings = IcebergToLinkMLSettings(  # pyright: ignore[reportCallIssue]
        catalog=CATALOG_NAME, output_dir=str(out_dir), group_by_namespace=group_by_namespace
    )
    with caplog.at_level(logging.WARNING):
        dump_catalog_schemas(settings)

    assert "Overwriting existing schema file" in caplog.text
    written = yaml_loader.load(str(stale_file), target_class=SchemaDefinition)
    assert written.name == ("ns" if group_by_namespace else "simple")


def test_dump_catalog_schemas_pass_namespace_isolation_and_collisions(
    catalog: Catalog, catalog_env: Path, simple_schema: Schema
) -> None:
    """Keep namespaces separate and retain colliding table and nested class definitions."""
    nested_schema = Schema(NestedField(1, "child", StructType(NestedField(2, "value", StringType()))))
    for namespace in ("first", "second"):
        catalog.create_namespace(namespace)
        catalog.create_table((namespace, "record"), schema=nested_schema, properties={"comment": "Table description."})
    catalog.create_table(("first", "record_child"), schema=simple_schema)
    catalog.create_table(("first", "a-b"), schema=simple_schema)
    catalog.create_table(("first", "a_b"), schema=nested_schema)
    out_dir = catalog_env / "schemas"
    settings = IcebergToLinkMLSettings(  # pyright: ignore[reportCallIssue]
        catalog=CATALOG_NAME, output_dir=str(out_dir), group_by_namespace=True
    )
    dump_catalog_schemas(settings)
    assert sorted(file.name for file in out_dir.iterdir()) == ["first.schema.yaml", "second.schema.yaml"]
    first = yaml_loader.load(str(out_dir / "first.schema.yaml"), target_class=SchemaDefinition)
    second = yaml_loader.load(str(out_dir / "second.schema.yaml"), target_class=SchemaDefinition)
    assert set(first.classes) == {"a_b", "a_b_2", "a_b_2_child", "record", "record_child", "record_child_2"}
    assert set(second.classes) == {"record", "record_child"}
    assert first.classes["a_b"].attributes == first.classes["record_child_2"].attributes
    assert first.classes["a_b_2"].attributes["child"].range == "a_b_2_child"
    assert first.classes["a_b_2_child"].attributes == first.classes["record_child"].attributes
    assert first.classes["record_child"].attributes == {
        "value": SlotDefinition(name="value", range="string"),
    }
    for namespace, document in (("first", first), ("second", second)):
        assert document.id == f"urn:iceberg:{namespace}"
        assert document.name == namespace
        assert document.classes["record"].description == "Table description."
        assert document.classes["record"].attributes == {"child": SlotDefinition(name="child", range="record_child")}
        assert list(Linter.validate_schema(str(out_dir / f"{namespace}.schema.yaml"))) == []
    assert first.annotations["generated_at"] == second.annotations["generated_at"]


@pytest.mark.parametrize(
    ("arguments", "expected"),
    [([], False), (["--group-by-namespace", "true"], True), (["--group-by-namespace", "false"], False)],
    ids=["default", "enabled", "disabled"],
)
def test_iceberg_to_linkml_settings_pass_grouping_cli(tmp_path: Path, arguments: list[str], expected: bool) -> None:
    """Parse namespace grouping from the CLI while retaining the per-table default."""
    settings = IcebergToLinkMLSettings(  # pyright: ignore[reportCallIssue]
        _cli_parse_args=["--catalog", CATALOG_NAME, "--output-dir", str(tmp_path), *arguments]
    )
    assert settings.group_by_namespace is expected
