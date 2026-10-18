"""Settings class for the xml2db-based XML ingestion pipeline for the KBase CTS."""

from typing import Annotated, Final

from pydantic import Field, FilePath, PositiveInt

from cdm_data_loaders.core.fields import (
    DEFAULT_XML_FILE_GLOB,
    FILE_GLOB,
    FileGlob,
    NonEmptyStr,
)
from cdm_data_loaders.core.settings import FileFormatIngestSettings, default_settings_with_shortcuts

PIPELINE_NAME: Final[str] = "xml2db"
DEFAULT_XML2DB_SHORT_NAME: Final[str] = "xml2db"
DEFAULT_XML2DB_CHUNK_SIZE: Final[int] = 1000


class Xml2DbSettings(FileFormatIngestSettings):
    """Settings for the xml2db-based XML ingestion pipeline.

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
            **FileFormatIngestSettings.model_config["cli_shortcuts"],
            FILE_GLOB: "g",
        },
    )

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
            description=("Maximum number of xml_tag elements parsed at a time. Only relevant if xml_tag is set."),
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
    file_glob: Annotated[
        FileGlob,
        Field(
            default=DEFAULT_XML_FILE_GLOB,
            description="Glob pattern for XML files inside settings.input_dir.",
        ),
    ]
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
        FilePath | None,
        Field(
            default=None,
            description=(
                "Optional path to a YAML file with xml2db model configuration overrides "
                "(e.g. per-table `reuse` settings). See xml2db's ModelConfig documentation."
            ),
        ),
    ]
    xsd_file: Annotated[
        FilePath,
        Field(
            description="Path to the XSD file describing the XML documents to be imported.",
        ),
    ]
