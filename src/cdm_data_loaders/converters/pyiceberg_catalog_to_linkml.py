"""CLI interface for creating a dump of table schemas from an Iceberg catalog."""

from datetime import UTC, datetime
from logging import Logger, getLogger
from pathlib import Path
from typing import Annotated, Final, cast

from linkml_runtime.dumpers import yaml_dumper
from linkml_runtime.linkml_model.meta import Annotation, Prefix, SchemaDefinition
from pydantic import Field
from pydantic_settings import SettingsConfigDict
from pyiceberg.catalog import load_catalog

from cdm_data_loaders.converters.emitters.linkml import LINKML_PREFIX
from cdm_data_loaders.converters.pyiceberg_to_linkml import table_to_linkml
from cdm_data_loaders.core.fields import (
    OUTPUT_DIR,
    NonEmptyStr,
)
from cdm_data_loaders.core.path_fields import LocalOutputDir
from cdm_data_loaders.core.settings import (
    LoggerSettings,
    default_settings_with_shortcuts,
)

logger: Logger = getLogger(__name__)

PIPELINE_NAME: Final[str] = "iceberg_catalog_linkml_converter"


class IcebergToLinkMLSettings(LoggerSettings):
    """Settings for generating a set of LinkML schemas from an Iceberg data catalog."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME, cli_shortcuts={OUTPUT_DIR: "o", "catalog": "c", "group_by_namespace": "g"}
    )

    catalog: Annotated[NonEmptyStr, Field(description="Name of the catalog to retrieve schemas from")]
    output_dir: LocalOutputDir
    group_by_namespace: bool = Field(default=False, description="Export one LinkML schema per namespace")


def _write_schema(schema: SchemaDefinition, out_file: Path, generated_at: str) -> None:
    """Write a timestamped schema, warning before replacing an existing file."""
    cast("dict[str, Annotation]", schema.annotations)["generated_at"] = Annotation(
        tag="generated_at", value=generated_at
    )
    schema_yaml = yaml_dumper.dumps(schema)
    if out_file.exists():
        logger.warning("Overwriting existing schema file %s", out_file)
    out_file.write_text(schema_yaml, encoding="utf-8")
    logger.info("Wrote schema to %s", out_file)


def dump_catalog_schemas(settings: IcebergToLinkMLSettings) -> None:
    """Iterate through a catalog and dump out the table schemas as LinkML.

    A table that fails to load or convert is logged and skipped; it does not
    abort the dump of the remaining tables. `generated_at` is stamped
    identically on every table in a single run (it reflects the run, not each
    table's individual read time). An existing output file for a table is
    overwritten, with a warning logged first. With ``group_by_namespace``,
    each nonempty namespace produces one file containing its table classes.

    :param settings: Catalog name and output directory.
    :returns: None.
    """
    generated_at = datetime.now(UTC).isoformat()
    failures = 0

    with load_catalog(settings.catalog) as catalog:
        namespaces = catalog.list_namespaces()
        logger.info("Dumping schemas for %d namespace(s) from catalog %r", len(namespaces), settings.catalog)
        for namespace in namespaces:
            namespace_name = ".".join(namespace)
            namespace_schema = (
                SchemaDefinition(
                    id=f"urn:iceberg:{namespace_name}",
                    name=namespace_name,
                    prefixes=[
                        Prefix(prefix_prefix="iceberg", prefix_reference=f"urn:iceberg:{namespace_name}#"),
                        Prefix(prefix_prefix="linkml", prefix_reference=LINKML_PREFIX),
                    ],
                    default_prefix="iceberg",
                    imports=["linkml:types"],
                )
                if settings.group_by_namespace
                else None
            )
            for identifier in sorted(catalog.list_tables(namespace)):
                try:
                    table = catalog.load_table(identifier)
                    schema_doc = table_to_linkml(table, identifier, schema=namespace_schema)
                except Exception:
                    failures += 1
                    logger.exception("Failed to convert table %s; skipping", identifier)
                    continue

                if namespace_schema is None:
                    _write_schema(
                        schema_doc, settings.output_dir / (".".join(identifier) + ".schema.yaml"), generated_at
                    )

            if namespace_schema is not None and namespace_schema.classes:
                _write_schema(namespace_schema, settings.output_dir / (namespace_name + ".schema.yaml"), generated_at)

    if failures:
        logger.warning("Completed with %d table(s) skipped due to conversion failures", failures)


if __name__ == "__main__":
    settings = IcebergToLinkMLSettings()  # pyright: ignore[reportCallIssue]
    dump_catalog_schemas(settings)
