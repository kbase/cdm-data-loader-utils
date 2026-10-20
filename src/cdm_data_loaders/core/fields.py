"""Common fields and defaults for data loading and wrangling pipelines."""

from enum import StrEnum
from typing import Annotated, Final

from frozendict import frozendict
from pydantic import Field, PositiveInt, StringConstraints

INPUT_MOUNT: Final[str] = "/input_dir"
OUTPUT_MOUNT: Final[str] = "/output_dir"

DEFAULT_JSONL_FILE_GLOB: Final[str] = "*.jsonl*"
DEFAULT_XML_FILE_GLOB: Final[str] = "*.xml*"
GZIP_SUFFIX: Final[str] = ".gz"

# output file formats
JSONL: Final[str] = "jsonl"
PARQUET: Final[str] = "parquet"


class LoaderFileFormatEnum(StrEnum):
    """Valid values for loader file format."""

    JSONL = JSONL
    PARQUET = PARQUET


# destination config block names
CTS: Final[str] = "cts"

# field names
BATCH_SIZE: Final[str] = "batch_size"
BUFFER_SIZE: Final[str] = "buffer_size"
DATASET_NAME: Final[str] = "dataset_name"
DLT_DEV_MODE: Final[str] = "dlt_dev_mode"
DISABLE_OUTPUT_COMPRESSION: Final[str] = "disable_output_compression"
FILE_GLOB: Final[str] = "file_glob"
INPUT_DIR: Final[str] = "input_dir"
LOADER_FILE_FORMAT: Final[str] = "loader_file_format"
LOG_CONFIG_FILE: Final[str] = "log_config_file"
LOG_INTERVAL: Final[str] = "log_interval"
MAX_TABLE_NESTING: Final[str] = "max_table_nesting"
OUTPUT_DIR: Final[str] = "output_dir"
PRESERVE_TABLE_NESTING: Final[str] = "preserve_table_nesting"
SAVE_RAW_RESPONSES: Final[str] = "save_raw_responses"
TABLE_NAME: Final[str] = "table_name"
USE_DESTINATION: Final[str] = "use_destination"
USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: Final[str] = "use_output_dir_for_pipeline_metadata"


MIN_START_AT: Final[int] = 1

# Default values for the common fields
DEFAULTS = frozendict(
    {
        BATCH_SIZE: 1000,
        BUFFER_SIZE: 100,
        DLT_DEV_MODE: False,
        DISABLE_OUTPUT_COMPRESSION: False,
        FILE_GLOB: "*",
        LOADER_FILE_FORMAT: LoaderFileFormatEnum.PARQUET,
        LOG_CONFIG_FILE: None,
        LOG_INTERVAL: 1000,
        MAX_TABLE_NESTING: 0,
        PRESERVE_TABLE_NESTING: False,
        SAVE_RAW_RESPONSES: False,
        USE_DESTINATION: CTS,
        USE_OUTPUT_DIR_FOR_PIPELINE_METADATA: False,
    }
)

NonEmptyStr = Annotated[str, StringConstraints(min_length=1)]

BatchSize = Annotated[
    PositiveInt,
    Field(
        default=DEFAULTS[BATCH_SIZE],
        description="Number of items per batch",
    ),
]
BufferSize = Annotated[
    PositiveInt,
    Field(
        default=DEFAULTS[BUFFER_SIZE],
        description="Number of rows to buffer per table before yielding a batch to the destination. Must be a positive integer.",
    ),
]
DatasetName = Annotated[NonEmptyStr, Field(description="The name of the dataset being produced")]
DltDevMode = Annotated[
    bool,
    Field(
        default=DEFAULTS[DLT_DEV_MODE],
        description="Enables DLT's dev mode. In dev mode, dlt resets the pipeline state on each run and creates a new dataset with a unique timestamped name for each run.",
    ),
]
DisableOutputCompression = Annotated[
    bool,
    Field(
        default=DEFAULTS[DISABLE_OUTPUT_COMPRESSION],
        description="Disables output compression for the pipeline. Output files are compressed by default; setting this to False turns off compression.",
    ),
]
FileGlob = Annotated[str, Field(default=DEFAULTS[FILE_GLOB], description="File glob for input files")]
InputDir = Annotated[
    NonEmptyStr,
    Field(
        # no default - value must be supplied
        description="Location of directory containing file(s) to import",
    ),
]
LoaderFileFormat = Annotated[
    LoaderFileFormatEnum,
    Field(
        default=DEFAULTS[LOADER_FILE_FORMAT],
        description=f"Format of the output files. Choices: {[member.value for member in LoaderFileFormatEnum]}",
    ),
]
LogConfigFile = Annotated[
    NonEmptyStr | None,
    Field(
        default=DEFAULTS[LOG_CONFIG_FILE],
        description="Location of configuration file for the logger",
    ),
]
LogInterval = Annotated[
    PositiveInt,
    Field(
        default=DEFAULTS[LOG_INTERVAL],
        description="How often (in number of processed entries) to emit a progress log message.",
    ),
]
MaxTableNesting = Annotated[
    int,
    Field(
        default=DEFAULTS[MAX_TABLE_NESTING],
        description="Maximum level of nesting of output datasets. For infinite nesting, set to 0.",
        ge=0,
    ),
]
OutputDir = Annotated[
    NonEmptyStr | None,
    Field(
        description="Location to save raw data files and, if use_output_dir_for_pipeline_metadata is set, pipeline metadata to. Has no effect on where the dlt destination loads data.",
    ),
]
PreserveTableNesting = Annotated[
    bool,
    Field(
        default=DEFAULTS[PRESERVE_TABLE_NESTING],
        description="Whether or not nested data should be flattened out into separate tables",
    ),
]
SaveRawResponses = Annotated[
    bool,
    Field(
        default=DEFAULTS[SAVE_RAW_RESPONSES],
        description="Whether or not to save the raw responses from the API to the output directory.",
    ),
]
TableName = Annotated[NonEmptyStr, Field(description="The name of the table to export parsed data to")]
UseDestination = Annotated[
    NonEmptyStr,
    Field(
        # no default - value must be supplied
        description="Name of the dlt destination to load data into. The dlt config must set destination.<name>.destination_type; the destination's location and credentials also come from the dlt config.",
    ),
]
UseOutputDirForPipelineMetadata = Annotated[
    bool,
    Field(
        default=DEFAULTS[USE_OUTPUT_DIR_FOR_PIPELINE_METADATA],
        description="Store pipeline metadata in `<output_dir>/.dlt_conf`. output_dir must be set and on the local filesystem.",
    ),
]
