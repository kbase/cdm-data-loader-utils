"""Settings class for the XML to PyArrow pipeline for the KBase CTS."""

from pathlib import Path
from typing import Annotated, ClassVar, Final, Self

import pyarrow as pa
import xmlschema
from pydantic import Field, FilePath, PositiveInt, PrivateAttr, model_validator
from pydantic_settings import SettingsConfigDict

from cdm_data_loaders.core.fields import (
    DatasetName,
    NonEmptyStr,
    TableName,
)
from cdm_data_loaders.core.settings import (
    XML_CLI_SHORTCUTS,
    XmlIngestSettings,
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


class XmlToPyArrowSettings(XmlIngestSettings):
    """Settings for the XML to PyArrow ingestion pipeline."""

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts=XML_CLI_SHORTCUTS,
    )

    _xsd_schema: xmlschema.XMLSchema = PrivateAttr()
    _pyarrow_schema: pa.Schema = PrivateAttr()
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


class PresetXmlToPyArrowSettings(XmlToPyArrowSettings):
    """Base settings for XML-to-PyArrow pipelines backed by a bundled XSD and precomputed PyArrow schema."""

    _PRECOMPUTED_SCHEMA: ClassVar[pa.Schema]

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

    def _resolve_pyarrow_schema(self) -> pa.Schema:
        """Return the precomputed PyArrow schema for this preset pipeline."""
        return self._PRECOMPUTED_SCHEMA


class UniProtXmlToPyArrowSettings(PresetXmlToPyArrowSettings):
    """Settings for the UniProt XML to PyArrow ingestion pipeline."""

    _PRECOMPUTED_SCHEMA: ClassVar[pa.Schema] = UNIPROT_SCHEMA

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=UNIPROT_PIPELINE_NAME,
        cli_shortcuts=XML_CLI_SHORTCUTS,
    )

    dataset_name: Annotated[
        DatasetName,
        Field(
            default=UNIPROT_PIPELINE_NAME,
            description="The name of the dataset being produced",
        ),
    ]
    xsd_file: Annotated[
        FilePath,
        Field(
            default=UNIPROT_XSD_FILE,
            description="Path to the UniProt XSD file describing the XML documents to be imported.",
        ),
    ]


class UniRefXmlToPyArrowSettings(PresetXmlToPyArrowSettings):
    """Settings for the UniRef XML to PyArrow ingestion pipeline."""

    _PRECOMPUTED_SCHEMA: ClassVar[pa.Schema] = UNIREF_SCHEMA

    model_config: SettingsConfigDict = default_settings_with_shortcuts(
        cli_prog_name=UNIREF_PIPELINE_NAME,
        cli_shortcuts=XML_CLI_SHORTCUTS,
    )

    dataset_name: Annotated[
        DatasetName,
        Field(
            default=UNIREF_PIPELINE_NAME,
            description="The name of the dataset being produced",
        ),
    ]
    xsd_file: Annotated[
        FilePath,
        Field(
            default=UNIREF_XSD_FILE,
            description="Path to the UniRef XSD file describing the XML documents to be imported.",
        ),
    ]
