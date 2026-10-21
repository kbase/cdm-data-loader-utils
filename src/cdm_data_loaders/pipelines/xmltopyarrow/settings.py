"""Settings class for the XML to PyArrow pipeline for the KBase CTS."""

from pathlib import Path
from typing import Annotated, Final, Self

import pyarrow as pa
import xmlschema
from pydantic import Field, FilePath, PositiveInt, PrivateAttr, model_validator
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import (
    DEFAULT_XML_FILE_GLOB,
    FILE_GLOB,
    DatasetName,
    FileGlob,
    NonEmptyStr,
    TableName,
)
from cdm_data_loaders.core.settings import (
    FileFormatIngestSettings,
    default_settings_with_shortcuts,
)
from cdm_data_loaders.readers.xsd import (
    DEFAULT_MAX_RECURSION_DEPTH,
    generate_pyarrow_schema,
    load_schema,
)
from cdm_data_loaders.schemas.uniprot.schema import SCHEMA as UNIPROT_SCHEMA
from cdm_data_loaders.schemas.uniref.schema import SCHEMA as UNIREF_SCHEMA

PIPELINE_NAME: Final[str] = "xmltopyarrow"
UNIPROT_PIPELINE_NAME: Final[str] = "uniprot"
UNIREF_PIPELINE_NAME: Final[str] = "uniref"

ENTRY_XML_TAG: Final[str] = "entry"
ENTRY_TABLE_NAME: Final[str] = "entry"

SCHEMAS_DIR: Final[Path] = Path(__file__).resolve().parents[2] / "schemas"
UNIPROT_XSD_FILE: Final[Path] = SCHEMAS_DIR / "uniprot" / "uniprot.xsd"
UNIREF_XSD_FILE: Final[Path] = SCHEMAS_DIR / "uniref" / "uniref.xsd"


class XmlToPyArrowSettings(FileFormatIngestSettings):
    """Settings for the XML to PyArrow ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts={
            **FileFormatIngestSettings.model_config["cli_shortcuts"],
            FILE_GLOB: "g",
        },
    )

    _xsd_schema: xmlschema.XMLSchema = PrivateAttr()
    _pyarrow_schema: pa.Schema = PrivateAttr()
    file_glob: Annotated[
        FileGlob,
        Field(
            default=DEFAULT_XML_FILE_GLOB,
            description="Glob pattern for XML files inside each entity's input subdirectory.",
        ),
    ]
    table_name: TableName
    xml_tag: Annotated[
        NonEmptyStr,
        Field(
            description="XML tag to capture the contents of",
        ),
    ]
    xsd_file: Annotated[
        FilePath,
        Field(
            description="Path to the XSD file describing the XML documents to be imported. Including the schema will ensure that the appropriate relationships (1-to-1 vs 1-to-many) are created between parent and child elements, and that element text is stored consistently.",
        ),
    ]
    max_recursion_depth: Annotated[
        PositiveInt,
        Field(
            default=DEFAULT_MAX_RECURSION_DEPTH,
            description="Maximum recursion depth for self-referential XSD types before serializing deeper content as JSON.",
        ),
    ]

    @model_validator(mode="after")
    def check_xsd_file(self) -> Self:
        """Check whether the schema file is present and valid; load the XSD and resolve the PyArrow schema."""
        try:
            self._xsd_schema = load_schema(self.xsd_file)
            self._pyarrow_schema = self._resolve_pyarrow_schema()
        except Exception as e:
            err_msg = f"Could not generate PyArrow schema from XSD: {e!s}"
            raise ValueError(err_msg) from e

        return self

    def _resolve_pyarrow_schema(self) -> pa.Schema:
        """Return the PyArrow schema for this settings instance."""
        return generate_pyarrow_schema(
            self._xsd_schema,
            self.xml_tag,
            max_recursion_depth=self.max_recursion_depth,
        )

    @property
    def xsd_schema(self) -> xmlschema.XMLSchema:
        """Return the loaded XMLSchema instance."""
        return self._xsd_schema

    @property
    def pyarrow_schema(self) -> pa.Schema:
        """Return the generated PyArrow schema for xml_tag."""
        return self._pyarrow_schema


class UniProtXmlToPyArrowSettings(XmlToPyArrowSettings):
    """Settings for the UniProt XML to PyArrow ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=UNIPROT_PIPELINE_NAME,
        cli_shortcuts={
            **FileFormatIngestSettings.model_config["cli_shortcuts"],
            FILE_GLOB: "g",
        },
    )

    dataset_name: Annotated[
        DatasetName,
        Field(
            default=UNIPROT_PIPELINE_NAME,
            description="The name of the dataset being produced",
        ),
    ]
    table_name: Annotated[
        TableName,
        Field(
            default=ENTRY_TABLE_NAME,
            description="The name of the table to export parsed data to",
        ),
    ]
    xml_tag: Annotated[
        NonEmptyStr,
        Field(
            default=ENTRY_XML_TAG,
            description="XML tag to capture the contents of",
        ),
    ]
    xsd_file: Annotated[
        FilePath,
        Field(
            default=UNIPROT_XSD_FILE,
            description="Path to the UniProt XSD file describing the XML documents to be imported.",
        ),
    ]

    def _resolve_pyarrow_schema(self) -> pa.Schema:
        """Return the precomputed UniProt PyArrow schema."""
        return UNIPROT_SCHEMA


class UniRefXmlToPyArrowSettings(XmlToPyArrowSettings):
    """Settings for the UniRef XML to PyArrow ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=UNIREF_PIPELINE_NAME,
        cli_shortcuts={
            **FileFormatIngestSettings.model_config["cli_shortcuts"],
            FILE_GLOB: "g",
        },
    )

    dataset_name: Annotated[
        DatasetName,
        Field(
            default=UNIREF_PIPELINE_NAME,
            description="The name of the dataset being produced",
        ),
    ]
    table_name: Annotated[
        TableName,
        Field(
            default=ENTRY_TABLE_NAME,
            description="The name of the table to export parsed data to",
        ),
    ]
    xml_tag: Annotated[
        NonEmptyStr,
        Field(
            default=ENTRY_XML_TAG,
            description="XML tag to capture the contents of",
        ),
    ]
    xsd_file: Annotated[
        FilePath,
        Field(
            default=UNIREF_XSD_FILE,
            description="Path to the UniRef XSD file describing the XML documents to be imported.",
        ),
    ]

    def _resolve_pyarrow_schema(self) -> pa.Schema:
        """Return the precomputed UniRef PyArrow schema."""
        return UNIREF_SCHEMA
