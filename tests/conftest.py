"""Global configuration settings for tests."""

import datetime
import logging
import shutil
from collections.abc import Generator
from copy import deepcopy
from importlib.util import find_spec
from pathlib import Path
from typing import Any, Final
from unittest.mock import patch

import boto3
import dlt
import pytest
from frozendict import frozendict
from moto import mock_aws
from pyspark.conf import SparkConf
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql.types import (
    ArrayType,
    BooleanType,
    DateType,
    FloatType,
    IntegerType,
    StringType,
    StructField,
    StructType,
)
from types_boto3_s3.client import S3Client

import cdm_data_loaders.utils.file_transfer.s3.client as s3_client
from cdm_data_loaders.audit.schema import (
    NAMESPACE,
    PIPELINE,
    ROW_ERRORS,
    RUN_ID,
    SOURCE,
)
from cdm_data_loaders.core.pipeline_run import PipelineRun
from cdm_data_loaders.readers.dsv import INVALID_DATA_FIELD
from cdm_data_loaders.utils.file_transfer.s3.client import _client_config, reset_s3_client

SAVE_DIR: Final[str] = "spark.sql.warehouse.dir"

TEST_NS: Final[str] = "test_ns"
PIPELINE_RUN = frozendict({RUN_ID: "1234-5678-90", PIPELINE: "KeystoneXL", SOURCE: "/path/to/file"})
ALT_PIPELINE_RUN = frozendict({RUN_ID: "9876-5432-10", PIPELINE: "KeystoneXXXL", SOURCE: "/path/to/dir"})

BASE_DIR: Final[Path] = Path("tests").parent
TEST_DATA_DIR: Final[Path] = Path("tests") / "data"
CASSETTES_DIR: Final[Path] = Path("tests") / "cassettes"

DEFAULT_VCR_CONFIG = frozendict(
    {
        "cassette_library_dir": str(CASSETTES_DIR),
        "record_mode": "once",  # record on first run, replay thereafter
        "serializer": "yaml",
        "match_on": ["method", "scheme", "host", "path", "query"],
        # strip the NCBI API key from cassettes
        "filter_query_parameters": ["api_key"],
        "filter_headers": ["api_key"],
        "decode_compressed_response": True,
        "allow_playback_repeats": True,
    }
)

_notebook_utils_available: bool | None = None


def _find_notebook_utils() -> bool:
    """Check whether the current environment has the KBase Lakehouse notebook_utils lib available.

    :return: True if they are available
    :rtype: bool
    """
    global _notebook_utils_available  # noqa: PLW0603
    if _notebook_utils_available is not None:
        return _notebook_utils_available

    # set to false by default
    _notebook_utils_available = False
    try:
        sss = find_spec("berdl_notebook_utils.setup_spark_session")
        db = find_spec("berdl_notebook_utils.spark.database")
        if sss is not None and db is not None:
            _notebook_utils_available = True
    except Exception:
        logging.getLogger(__name__).exception("Notebook utils not available: requires_notebook_utils tests will fail")

    return _notebook_utils_available


def pytest_runtest_setup(item: pytest.Item) -> None:
    """Set up tests to skip anything where dependencies are missing.

    :param item: test item
    :type item: pytest.Item
    """
    if "requires_spark" not in item.keywords:
        return

    # notebook utils are available: we're fine
    if _find_notebook_utils():
        return

    # the raw string passed to -m
    markexpr = item.config.option.markexpr
    # tests matching any "not" markers will already have been filtered out
    # markexpr will contain "requires_spark" when -m requires_spark is used
    # N.b. this blatantly ignores more complex situations, like "not(marker or marker)"
    if markexpr and "requires_spark" in markexpr and "not requires_spark" not in markexpr:
        pytest.fail(
            "Test is marked requires_spark, but Spark is not available in this environment.",
            pytrace=False,
        )
    else:
        # no markexpr or unrelated -m expression: skip silently
        pytest.skip("Test is marked requires_spark, but Spark is not available in this environment.")


@pytest.fixture(autouse=True)
def logging_setup(caplog: pytest.LogCaptureFixture) -> None:
    """Fiddle with the loggers used in the tests for a better experience."""
    vcr_logger = logging.getLogger("vcr")
    vcr_logger.setLevel("ERROR")
    # turn on log propagation for the dlt logger
    dlt_logger = logging.getLogger("dlt")
    dlt_logger.propagate = True
    caplog.set_level(logging.INFO)
    caplog.clear()


@pytest.fixture
def spark(tmp_path: Path) -> Generator[SparkSession, Any]:
    """Generate a spark session with spark.sql.warehouse.dir set to the pytest temporary directory."""
    logger = logging.getLogger(__name__)
    try:
        from berdl_notebook_utils.setup_spark_session import generate_spark_conf  # pyright: ignore[reportMissingImports]  # noqa: I001, PLC0415
    except ModuleNotFoundError:
        logger.exception("berdl_notebook_utils not available: cannot create a spark session")
        raise

    config = generate_spark_conf("test_delta_app", local=True)
    test_config = {
        "spark.sql.shuffle.partitions": 5,
        "spark.default.parallelism": 9,
        # Disabling rdd and map output compression as data is already small for tests
        "spark.rdd.compress": False,
        "spark.shuffle.compress": False,
        # Disable Spark UI for tests
        "spark.ui.enabled": False,
        "spark.ui.showConsoleProgress": False,
        # Extra configs to optimize Delta internal operations on tests
        "spark.databricks.delta.snapshotPartitions": 2,
        "delta.log.cacheSize": 3,
    }
    config.update(test_config)
    config[SAVE_DIR] = str(tmp_path)
    spark_conf = SparkConf().setAll(list(config.items()))
    logger.info("starting spark session...")
    spark = SparkSession.builder.config(conf=spark_conf).getOrCreate()
    save_dir = spark.conf.get(SAVE_DIR).removeprefix("file:")  # pyright: ignore[reportOptionalMemberAccess]
    if save_dir != str(tmp_path):
        logger.error("spark dir: %s; save dir: %s", tmp_path, save_dir)
    yield spark
    spark.catalog.clearCache()
    spark.stop()
    logger.info("stopping spark session...")
    shutil.rmtree(save_dir)


@pytest.fixture(scope="session")
def test_data_dir() -> Path:
    """Test data directory."""
    return TEST_DATA_DIR


@pytest.fixture(scope="session")
def json_test_strings() -> dict[str, Any]:
    """A selection of JSON strings for testing."""
    return {
        "null": "null",
        "empty": "",
        "empty_str": '""',
        "str": "some random string",
        "quoted_str": '"some random string"',
        "ws": "\n\n\t\n\t   ",
        "quoted_ws": '"\n\n\t\n\t   "',
        "empty_object": "{}",
        "empty_array": "[]",
        "null_key": '{"key": {null: "value"}}',
        "empty_key": '{"key": {"": "value"}}',
        "unclosed_str": '{"key": "value}',
        "array_null": "[null]",
        "array_mixed": '[null, "", 0, 1, 1.2345, "string", ["that", "this"], {"this": "that"}]',
        "array_of_str": '["this", "that", "the", "other"]',
        "array_of_arrays": '[[1,2,3],["this","that"],["what?"]]',
        "array_of_objects": '[{"key": "value"}]',
        "object": '{"key": "value"}',
    }


CONFIG_BUCKET = frozendict({"local_fs": "/output_dir", "s3": "s3://some/s3/bucket"})

TEST_DLT_CONFIG = frozendict(
    {
        "destination.local_fs.bucket_url": CONFIG_BUCKET["local_fs"],
        "destination.local_fs.destination_type": "filesystem",
        "destination.s3.bucket_url": CONFIG_BUCKET["s3"],
        "destination.s3.destination_type": "filesystem",
        "normalize.data_writer.disable_compression": False,
    }
)


def _generate_dlt_config() -> dict[str, Any]:
    """Return a fresh DLT config dict (same shape as the conftest fixture)."""
    return {
        "destination": {
            "local_fs": {"bucket_url": CONFIG_BUCKET["local_fs"]},
            "s3": {"bucket_url": CONFIG_BUCKET["s3"]},
        },
        "destination.local_fs.bucket_url": CONFIG_BUCKET["local_fs"],
        "destination.s3.bucket_url": CONFIG_BUCKET["s3"],
        "normalize.data_writer.disable_compression": False,
    }


@pytest.fixture
def dlt_config() -> dict[str, Any]:
    """DLT config for testing purposes."""
    return _generate_dlt_config()


@pytest.fixture
def patch_dlt_config(dlt_config: dict[str, Any]) -> Generator[Any]:
    """Monkeypatch `dlt.config` to the test dlt config."""
    with pytest.MonkeyPatch.context() as m:
        m.setattr(dlt, "config", dlt_config)
        yield


@pytest.fixture
def empty_df_schema() -> list[StructField]:
    """List of fields corresponding to the empty dataframe."""
    return [
        StructField("name", StringType(), nullable=True),
        StructField("id", IntegerType(), nullable=False),
    ]


@pytest.fixture
def empty_df(spark: SparkSession, empty_df_schema: list[StructField]) -> Generator[DataFrame, Any]:
    """Empty dataframe for testing usage."""
    df = spark.createDataFrame([], schema=StructType(empty_df_schema))
    yield df
    assert df.schema == StructType(empty_df_schema)


# Various CSV permutations for testing
VALID = "valid_csv"
MISSING_REQUIRED = "invalid_csv_missing_required"
TYPE_MISMATCH = "invalid_csv_type_mismatch"
TOO_FEW_COLS = "invalid_csv_too_few_cols"
TOO_MANY_COLS = "invalid_csv_too_many_cols"
ALL_LINES = "all_lines"


@pytest.fixture
def csv_schema() -> list[StructField]:
    """List of fields for parsing the various CSV snippets."""
    return [
        StructField("col1", IntegerType(), nullable=False),
        StructField("col2", DateType(), nullable=False),
        StructField("col3", FloatType(), nullable=False),
        StructField("col4", BooleanType(), nullable=False),
        StructField("col5", StringType(), nullable=False),
    ]


@pytest.fixture(scope="session")
def valid_csv(test_data_dir: Path) -> Path:
    """Valid CSV data.

    1,20250301,1.2345,true,EcoCyc:EG10986-MONOMER
    2,20250201,0.2,false,MetaCyc:EG10986-MONOMER
    3,20250801,23,True,4261555
    4,00010101,.1234,False,col5

    """
    return test_data_dir / "dsv" / "valid.csv"


@pytest.fixture(scope="session")
def invalid_csv_missing_required(test_data_dir: Path) -> Path:
    """CSV data with required fields missing.

    # correct number of cols, but some cols are empty
    1,,,,col5
    # missing leading cols
    ,,2.345,True,col5
    # missing trailing cols
    3,20250531,23.45,,
    # all missing
    ,,,,
    """
    return test_data_dir / "dsv" / "missing_required.csv"


@pytest.fixture
def invalid_csv_missing_required_annots() -> list[list[str]]:
    """Generate the expected error annotations for the lines in invalid_csv_missing_required.

    :return: list of list of error strings
    :rtype: list[list[str]]
    """
    valid_invalid_fields = [[1, 0, 0, 0, 1], [0, 0, 1, 1, 1], [1, 1, 1, 0, 0], [0, 0, 0, 0, 0]]
    return [[f"missing_required: col{n + 1}" for n in range(5) if not row[n]] for row in valid_invalid_fields]


@pytest.fixture(scope="session")
def invalid_csv_type_mismatch(test_data_dir: Path) -> Path:
    """CSV data with incorrect data types.

    # Y, N, Y, N, Y
    1,2,3,4,5
    # N, N, Y, N, Y
    1.234,2.3456,3.45,4.5,5
    # N, N, N, Y, Y
    true,false,true,false,true
    # Y, Y, N, N, Y
    00200202,00200202,00200202,00200202,00200202
    """
    return test_data_dir / "dsv" / "type_mismatch.csv"


@pytest.fixture(scope="session")
def invalid_csv_too_few_cols(test_data_dir: Path) -> Path:
    """CSV data containing rows with too few columns.

    # too few cols
    1
    2,20250502,0.2345
    3,,23.56,False
    ,,,
    """
    return test_data_dir / "dsv" / "too_few_cols.csv"


@pytest.fixture(scope="session")
def invalid_csv_too_many_cols(test_data_dir: Path) -> Path:
    """CSV data containing rows with too many columns.

    # too many cols - all have 6 cols
    ,,,,,
    2,20250710,col3,True,,col6
    # empty trailing
    3,20250101,,,,
    # empty leading
    ,,,,,col6
    """
    return test_data_dir / "dsv" / "too_many_cols.csv"


@pytest.fixture(scope="session")
def all_lines(test_data_dir: Path) -> Path:
    """All the CSV lines in a single fixture!"""
    return test_data_dir / "dsv" / "all_lines.csv"


@pytest.fixture(scope="session")
def all_lines_tsv(test_data_dir: Path) -> Path:
    """All the CSV lines in a single fixture!"""
    return test_data_dir / "dsv" / "all_lines.tsv"


@pytest.fixture
def annotated_df_schema(csv_schema: list[StructField]) -> StructType:
    """The schema for the annotated dataframe produced by validating one of the CSV files above."""
    actual_csv_schema = list(csv_schema)
    for r in actual_csv_schema:
        r.nullable = True

    return StructType([*actual_csv_schema, INVALID_DATA_FIELD, StructField(ROW_ERRORS, ArrayType(StringType()))])


@pytest.fixture(scope="session")
def annotated_df_data() -> list[dict[str, Any]]:
    """The output of running all_lines (above) through the dsv parser and df_nullable_fields validator."""
    return deepcopy(
        [
            {
                "col1": 1,
                "col2": datetime.date(2025, 3, 1),
                "col3": 1.2345000505447388,
                "col4": True,
                "col5": "EcoCyc:EG10986-MONOMER",
                "__invalid_data__": None,
                "errors_in_record": [],
            },
            {
                "col1": 2,
                "col2": datetime.date(2025, 2, 1),
                "col3": 0.20000000298023224,
                "col4": False,
                "col5": "MetaCyc:EG10986-MONOMER",
                "__invalid_data__": None,
                "errors_in_record": [],
            },
            {
                "col1": 3,
                "col2": datetime.date(2025, 8, 1),
                "col3": 23.0,
                "col4": True,
                "col5": "4261555",
                "__invalid_data__": None,
                "errors_in_record": [],
            },
            {
                "col1": 4,
                "col2": datetime.date(1, 1, 1),
                "col3": 0.1234000027179718,
                "col4": False,
                "col5": "col5",
                "__invalid_data__": None,
                "errors_in_record": [],
            },
            {
                "col1": 1,
                "col2": None,
                "col3": None,
                "col4": None,
                "col5": "col5",
                "__invalid_data__": None,
                "errors_in_record": [
                    "missing_required: col2",
                    "missing_required: col3",
                    "missing_required: col4",
                ],
            },
            {
                "col1": None,
                "col2": None,
                "col3": 2.3450000286102295,
                "col4": True,
                "col5": "col5",
                "__invalid_data__": None,
                "errors_in_record": ["missing_required: col1", "missing_required: col2"],
            },
            {
                "col1": 3,
                "col2": datetime.date(2025, 5, 31),
                "col3": 23.450000762939453,
                "col4": None,
                "col5": None,
                "__invalid_data__": None,
                "errors_in_record": ["missing_required: col4", "missing_required: col5"],
            },
            {
                "col1": None,
                "col2": None,
                "col3": None,
                "col4": None,
                "col5": None,
                "__invalid_data__": None,
                "errors_in_record": [
                    "missing_required: col1",
                    "missing_required: col2",
                    "missing_required: col3",
                    "missing_required: col4",
                    "missing_required: col5",
                ],
            },
            {
                "col1": None,
                "col2": None,
                "col3": None,
                "col4": None,
                "col5": None,
                "__invalid_data__": ",,,,,",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": 2,
                "col2": datetime.date(2025, 7, 10),
                "col3": None,
                "col4": True,
                "col5": None,
                "__invalid_data__": "2,20250710,col3,True,,col6",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": 3,
                "col2": datetime.date(2025, 1, 1),
                "col3": None,
                "col4": None,
                "col5": None,
                "__invalid_data__": "3,20250101,,,,",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": None,
                "col2": None,
                "col3": None,
                "col4": None,
                "col5": None,
                "__invalid_data__": ",,,,,col6",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": 1,
                "col2": None,
                "col3": None,
                "col4": None,
                "col5": None,
                "__invalid_data__": "1",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": 2,
                "col2": datetime.date(2025, 5, 2),
                "col3": 0.2345000058412552,
                "col4": None,
                "col5": None,
                "__invalid_data__": "2,20250502,0.2345",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": 3,
                "col2": None,
                "col3": 23.559999465942383,
                "col4": False,
                "col5": None,
                "__invalid_data__": "3,,23.56,False",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": None,
                "col2": None,
                "col3": None,
                "col4": None,
                "col5": None,
                "__invalid_data__": ",,,",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": 1,
                "col2": None,
                "col3": 3.0,
                "col4": None,
                "col5": "5",
                "__invalid_data__": "1,2,3,4,5",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": None,
                "col2": None,
                "col3": 3.450000047683716,
                "col4": None,
                "col5": "5",
                "__invalid_data__": "1.234,2.3456,3.45,4.5,5",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": None,
                "col2": None,
                "col3": None,
                "col4": False,
                "col5": "true",
                "__invalid_data__": "true,false,true,false,true",
                "errors_in_record": ["parse_error"],
            },
            {
                "col1": 200202,
                "col2": datetime.date(20, 2, 2),
                "col3": 200202.0,
                "col4": None,
                "col5": "00200202",
                "__invalid_data__": "00200202,00200202,00200202,00200202,00200202",
                "errors_in_record": ["parse_error"],
            },
        ]
    )


@pytest.fixture(scope="session")
def annotated_df_errors(annotated_df_data: list[dict[str, Any]]) -> set[str]:
    """Unique errors found in annotated_df_data."""
    list_of_errors = []
    for r in annotated_df_data:
        list_of_errors.extend(r[ROW_ERRORS])

    return set(list_of_errors)


# Audit-related stuff
@pytest.fixture(scope="package")
def pipeline_run() -> PipelineRun:
    """Generate a pipeline run."""
    return PipelineRun(**{**PIPELINE_RUN, NAMESPACE: TEST_NS})


@pytest.fixture(scope="package")
def alt_pipeline_run() -> PipelineRun:
    """Generate a different pipeline run."""
    return PipelineRun(**{**ALT_PIPELINE_RUN, NAMESPACE: TEST_NS})


"""S3 Client mocks"""


TEST_BUCKET: Final[str] = "test_bucket"
ALT_BUCKET: Final[str] = "alt_bucket"
BUCKETS = [TEST_BUCKET, ALT_BUCKET]


@pytest.fixture
def mock_s3_client(monkeypatch: pytest.MonkeyPatch) -> Generator[S3Client, Any]:
    """Yield a mocked S3 client with both valid buckets created.

    The function get_s3_client() is patched to ensure that all module functions use this client.

    Resets the cached client before and after to prevent state leaking between tests.
    """
    # Remove any real endpoint/credential env vars so moto intercepts all HTTP calls.
    monkeypatch.setenv("AWS_ACCESS_KEY_ID", "testing")
    monkeypatch.setenv("AWS_SECRET_ACCESS_KEY", "testing")
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-east-1")
    monkeypatch.delenv("AWS_ENDPOINT_URL", raising=False)
    monkeypatch.delenv("AWS_ENDPOINT_URL_S3", raising=False)
    boto3.DEFAULT_SESSION = None

    with mock_aws():
        reset_s3_client()
        client: S3Client = boto3.client("s3", config=_client_config())

        for bucket in BUCKETS:
            client.create_bucket(Bucket=bucket)

        # delete any existing client
        reset_s3_client()
        assert s3_client._s3_client is None  # noqa: SLF001

        # patch in the client that we have just created
        with patch.object(s3_client, "get_s3_client", return_value=client):
            yield client

        reset_s3_client()
        assert s3_client._s3_client is None  # noqa: SLF001
