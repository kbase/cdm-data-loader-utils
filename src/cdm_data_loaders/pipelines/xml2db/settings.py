"""Settings class for the xml2db-based XML ingestion pipeline for the KBase CTS."""

from pathlib import Path
from typing import Annotated, Final

from pydantic import Field, PositiveInt, field_validator

from cdm_data_loaders.core.fields import (
    DEFAULT_XML_FILE_GLOB,
    FILE_GLOB,
    BufferSize,
    DatasetName,
    FileGlob,
    LoaderFileFormat,
    LogInterval,
    NonEmptyStr,
)
from cdm_data_loaders.core.settings import CLI_SHORTCUTS, CtsSettings, default_settings_with_shortcuts

PIPELINE_NAME: Final[str] = "xml2db"
DEFAULT_XML2DB_SHORT_NAME: Final[str] = "xml2db"
DEFAULT_XML2DB_CHUNK_SIZE: Final[int] = 1000


class Xml2DbSettings(CtsSettings):
    """Settings for the xml2db-based XML ingestion pipeline.

    Unlike `XmlToDictSettings`, there is no `xml_tag`/`table_name` to configure: xml2db derives
    the full set of output tables (and their names) from the XSD itself, so the only schema-related
    input required is the path to that XSD.

    xml2db loads an entire XML document into memory before it can be converted into rows, which
    risks OOM for very large files. If the files being imported follow the shape of a single root
    element directly followed by many repeated elements of one type (e.g. `<UniRef50><entry/>
    ...<entry/></UniRef50>`), set `xml_tag` to that repeated element's tag: the file will
    then be streamed and parsed in bounded-size pieces of `chunk_size` elements at a time instead
    of all at once. See `cdm_data_loaders.readers.xml2db_doc.iter_xml2db_chunk_fragments` for the
    details and constraints of this mode. Note that xml2db's content-hash deduplication of
    "reused" tables only happens within a single parsed document, so chunking (like parsing
    multiple separate files) trades whole-run deduplication for bounded memory: rows in tables
    shared across many elements (e.g. UniRef's `property` table) may be duplicated once per
    document (file or chunk) they occur in, rather than collapsed into a single row.
    `compact_reused_tables` cleans this up after the fact, regardless of chunking.
    """

    model_config = default_settings_with_shortcuts(
        cli_prog_name=PIPELINE_NAME,
        cli_shortcuts={
            **CLI_SHORTCUTS,
            FILE_GLOB: "g",
        },
    )

    buffer_size: BufferSize
    xml_tag: Annotated[
        NonEmptyStr | None,
        Field(
            default=None,
            description=(
                "Repeated top-level child element to stream and parse in bounded-size pieces "
                "(e.g. 'entry', or a namespaced tag like '{http://uniprot.org/uniref}entry'). "
                "If unset (the default), each file is parsed with xml2db in a single pass, which "
                "loads the whole document into memory. Only local names are matched; the "
                "namespace portion of a namespaced tag, if given, is ignored."
            ),
        ),
    ]
    chunk_size: Annotated[
        PositiveInt,
        Field(
            default=DEFAULT_XML2DB_CHUNK_SIZE,
            description=(
                "Maximum number of xml_tag elements parsed by xml2db at a time when xml_tag is set. Ignored otherwise."
            ),
        ),
    ]
    compact_reused_tables: Annotated[
        bool,
        Field(
            default=True,
            description=(
                "After a run, deduplicate literal duplicate rows out of every xml2db 'reused' "
                "table (see cdm_data_loaders.pipelines.xml2db.compaction for why parsing more "
                "than one file and/or chunk can produce them, and why this is always a safe, "
                "reference-preserving row removal). Only takes effect when use_destination is "
                "'local_fs'; ignored otherwise."
            ),
        ),
    ]
    dataset_name: DatasetName
    file_glob: Annotated[
        FileGlob,
        Field(
            default=DEFAULT_XML_FILE_GLOB,
            description="Glob pattern for XML files inside settings.input_dir.",
        ),
    ]
    loader_file_format: LoaderFileFormat
    log_interval: LogInterval
    short_name: Annotated[
        NonEmptyStr,
        Field(
            default=DEFAULT_XML2DB_SHORT_NAME,
            description=(
                "Short identifier for the xml2db data model; used to name the virtual root table "
                "when the XSD declares more than one top-level element."
            ),
        ),
    ]
    skip_xml_validation: Annotated[
        bool,
        Field(
            default=True,
            description=(
                "Skip validating each XML document against the XSD before parsing it. "
                "Validation is slower but will raise on schema-invalid documents."
            ),
        ),
    ]
    xml2db_config_file: Annotated[
        NonEmptyStr | None,
        Field(
            default=None,
            description=(
                "Optional path to a YAML file with xml2db model configuration overrides "
                "(e.g. per-table `reuse` settings). See xml2db's ModelConfig documentation."
            ),
        ),
    ]
    xsd_file: Annotated[
        NonEmptyStr,
        Field(
            description="Path to the XSD file describing the XML documents to be imported.",
        ),
    ]

    @field_validator("xsd_file", "xml2db_config_file")
    @classmethod
    def validate_file_exists(cls, v: str | None) -> str | None:
        """Validate that a file path setting, if supplied, points to an existing file.

        :param v: the file path to validate
        :type v: str | None
        :raises ValueError: if the path is supplied but does not point to an existing file.
        :return: the validated file path
        :rtype: str | None
        """
        if v is not None and not Path(v).is_file():
            err_msg = f"File does not exist: {v}"
            raise ValueError(err_msg)
        return v
