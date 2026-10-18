"""Global configuration settings for tests."""

import gzip
import logging
import os
import shutil
import sys
from collections.abc import Callable, Generator
from copy import deepcopy
from dataclasses import dataclass
from importlib.util import find_spec
from pathlib import Path
from typing import Any, Final, Self
from unittest.mock import patch

import boto3
import dlt
import pytest
from frozendict import frozendict
from moto import mock_aws
from pyspark.sql import SparkSession
from types_boto3_s3.client import S3Client

import cdm_data_loaders.utils.file_transfer.s3.client as s3_client
from cdm_data_loaders.utils import cdm_logger
from cdm_data_loaders.utils.file_transfer.s3.client import _client_config, reset_s3_client
from tests.dlt_config_isolation import isolated_dlt_config

pytest_plugins = ["tests.jsonlines_fixtures", "tests.jsonschema_fixtures"]

LOCAL_FS: Final[str] = "local_fs"
S3: Final[str] = "s3"

SAVE_DIR: Final[str] = "spark.sql.warehouse.dir"

BASE_DIR: Final[Path] = Path("tests").parent
TEST_DATA_DIR: Final[Path] = Path("tests") / "data"
CASSETTES_DIR: Final[Path] = Path("tests") / "cassettes"

REFERENCE_XML_NS: Final[str] = "http://uniprot.org/uniref"
REFERENCE_XML_TAG: Final[str] = f"{{{REFERENCE_XML_NS}}}entry"
REFERENCE_XML_FIXTURE_DIR: Final[Path] = TEST_DATA_DIR / "uniprot" / "uniref"
REFERENCE_XSD: Final[Path] = REFERENCE_XML_FIXTURE_DIR / "uniref.xsd"
N_REFERENCE_XML_ENTRIES: Final[int] = 100

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
    vcr_logger.setLevel(logging.ERROR)
    # turn on log propagation for the dlt logger
    dlt_logger = logging.getLogger("dlt")
    dlt_logger.propagate = True
    caplog.set_level(logging.INFO)
    caplog.clear()


@pytest.fixture(autouse=True)
def _isolated_cli_and_env(monkeypatch: pytest.MonkeyPatch) -> Generator[None]:
    """Clear sys.argv and the environment to prevent test-to-test pollution."""
    monkeypatch.setattr(sys, "argv", ["pytest"])

    current_env = deepcopy(os.environ)
    with patch.dict(os.environ, {k: v for k, v in current_env.items() if not k.lower().startswith("cdl_")}, clear=True):
        yield


@pytest.fixture
def spark(tmp_path: Path) -> Generator[SparkSession, Any]:
    """Generate a local spark session with spark.sql.warehouse.dir set to the pytest temporary directory."""
    logger = logging.getLogger(__name__)
    try:
        from berdl_notebook_utils.setup_spark_session import get_spark_session  # pyright: ignore[reportMissingImports]  # noqa: PLC0415
    except ModuleNotFoundError:
        logger.exception("berdl_notebook_utils not available: cannot create a spark session")
        raise

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
        SAVE_DIR: str(tmp_path),
    }
    # local=True builds the session without any BERDL environment (no get_settings(),
    # no Spark Connect server); the override applies the test-tuning configs after.
    logger.info("starting spark session...")
    spark = get_spark_session("test_app", local=True, override=test_config)
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


@pytest.fixture
def write_gzip_file() -> Callable[[Path, str, str], Path]:
    """Return a function that writes text content to a gzip-compressed file."""

    def _write(directory: Path, filename: str, content: str) -> Path:
        directory.mkdir(parents=True, exist_ok=True)
        file_path = directory / filename
        with gzip.open(file_path, "wb") as f:
            f.write(content.encode("utf-8"))
        return file_path

    return _write


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


CONFIG_BUCKET = frozendict({LOCAL_FS: "/output_dir", S3: "s3://some/s3/bucket"})

TEST_DLT_CONFIG = frozendict(
    {
        "destination.local_fs.bucket_url": CONFIG_BUCKET[LOCAL_FS],
        "destination.local_fs.destination_type": "filesystem",
        "destination.s3.bucket_url": CONFIG_BUCKET[S3],
        "destination.s3.destination_type": "filesystem",
        "normalize.data_writer.disable_compression": False,
    }
)


def _generate_dlt_config() -> dict[str, Any]:
    """Return a fresh DLT config dict (same shape as the conftest fixture)."""
    return {
        "destination": {
            LOCAL_FS: {"bucket_url": CONFIG_BUCKET[LOCAL_FS], "destination_type": "filesystem"},
            S3: {"bucket_url": CONFIG_BUCKET[S3], "destination_type": "filesystem"},
        },
        "destination.local_fs.bucket_url": CONFIG_BUCKET[LOCAL_FS],
        "destination.local_fs.destination_type": "filesystem",
        "destination.s3.bucket_url": CONFIG_BUCKET[S3],
        "destination.s3.destination_type": "filesystem",
        "normalize.data_writer.disable_compression": False,
    }


@pytest.fixture
def dlt_config() -> dict[str, Any]:
    """DLT config for testing purposes."""
    return _generate_dlt_config()


@pytest.fixture(autouse=True)
def _isolated_dlt_config() -> Generator[None]:
    """Isolate every test's dlt config/secrets provider chain, seeded with the standard test config.

    Prevents mutation of the global dlt.config during tests. dlt.config can be set explicitly by
    passing dlt_config=... to a Settings object or using `isolated_dlt_config(...)`.
    """
    with isolated_dlt_config(_generate_dlt_config()):
        yield


@pytest.fixture
def dlt_destination_config(tmp_path: Path) -> Generator[str]:
    """Set the 'local_fs' destination to tmp_path in dlt.config."""
    bucket = tmp_path / "bucket"
    bucket.mkdir()
    with dlt.config.values(
        {
            "destination.local_fs.destination_type": "filesystem",
            "destination.local_fs.bucket_url": str(bucket),
        }
    ):
        yield LOCAL_FS


"""S3 Client mocks"""


TEST_BUCKET: Final[str] = "test_bucket"
ALT_BUCKET: Final[str] = "alt_bucket"
BUCKETS = [TEST_BUCKET, ALT_BUCKET]


@pytest.fixture
def mock_s3_client(monkeypatch: pytest.MonkeyPatch) -> Generator[S3Client]:
    """Yield a mocked S3 client with two valid buckets created (TEST_BUCKET, ALT_BUCKET).

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
        client: S3Client = boto3.client(S3, config=_client_config())  # pyright: ignore[reportCallIssue]

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


@dataclass(frozen=True, slots=True)
class LoggerState:
    """Settings of a single logger that tests may change."""

    level: int = logging.NOTSET
    propagate: bool = True
    disabled: bool = False
    filters: tuple[Any, ...] = ()
    handlers: tuple[logging.Handler, ...] = ()
    formatters: tuple[logging.Formatter | None, ...] = ()

    @classmethod
    def capture(cls, lg: logging.Logger, exclude: logging.Handler | None = None) -> Self:
        """Record the current state of ``lg``, leaving out the ``exclude`` handler."""
        handlers = tuple(h for h in lg.handlers if h is not exclude)
        return cls(
            level=lg.level,
            propagate=lg.propagate,
            disabled=lg.disabled,
            filters=tuple(lg.filters),
            handlers=handlers,
            formatters=tuple(h.formatter for h in handlers),
        )

    def restore(self, lg: logging.Logger, *extra_handlers: logging.Handler) -> None:
        """Put ``lg`` back into the recorded state and attach ``extra_handlers``."""
        lg.setLevel(self.level)
        lg.propagate = self.propagate
        lg.disabled = self.disabled
        lg.filters[:] = self.filters
        lg.handlers[:] = [*self.handlers, *extra_handlers]
        for handler, formatter in zip(self.handlers, self.formatters, strict=True):
            handler.setFormatter(formatter)


@pytest.fixture
def isolated_logging(monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture) -> Generator[None]:
    """Reset init_logger's record of the applied config and undo any logging changes made by the test."""
    root = logging.getLogger()
    manager = logging.Logger.manager
    # pytest attaches a new caplog handler to the root logger for each test phase, so this one must not be restored
    root_state = LoggerState.capture(root, exclude=caplog.handler)
    saved = {name: LoggerState.capture(lg) for name, lg in manager.loggerDict.items() if isinstance(lg, logging.Logger)}
    monkeypatch.setattr(cdm_logger, "applied_key", None)

    yield

    loggers = {name: lg for name, lg in manager.loggerDict.items() if isinstance(lg, logging.Logger)}
    known = {*root_state.handlers, caplog.handler, *(h for state in saved.values() for h in state.handlers)}
    for lg in [root, *loggers.values()]:
        for handler in lg.handlers:
            if handler not in known:
                handler.close()
    root_state.restore(root, caplog.handler)
    for name, lg in loggers.items():
        saved.get(name, LoggerState()).restore(lg)
