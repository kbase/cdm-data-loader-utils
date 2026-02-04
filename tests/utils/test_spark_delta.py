"""Tests of the Spark / Delta utilities."""

import logging
from collections.abc import Generator
from pathlib import Path
from typing import Any

import pytest
from pyspark.sql import DataFrame, DataFrameWriter, Row, SparkSession
from tests.helpers import assertDataFrameEqual

from cdm_data_loader_utils.utils import spark_delta
from cdm_data_loader_utils.utils.spark_delta import (
    APPEND,
    DEFAULT_APP_NAME,
    DEFAULT_NAMESPACE,
    ERROR,
    ERROR_IF_EXISTS,
    IGNORE,
    OVERWRITE,
    WRITE_MODE,
    get_spark,
    preview_or_skip,
    write_delta,
)

original_set_up_ws_fn = spark_delta.set_up_workspace

SAVE_DIR = "spark.sql.warehouse.dir"
DEFAULT_WRITE_MODE = ERROR
DEFAULT_SAMPLE_DATA = {"a": "A1", "b": "B1"}
TENANT_NAME = "The_Breakers"


def gen_ns_save_dir(current_save_dir: str, namespace: str, tenant_name: str | None) -> tuple[str, str]:
    """Generate the projected namespace and save directory, given a file path, a namespace, and a tenant name."""
    db_location = f"tenant/{tenant_name}/{namespace}.db" if tenant_name else f"user/some_user/{namespace}.db"
    namespace = db_location.replace("/", "__").replace(".db", "")
    save_dir = f"{current_save_dir.replace('file:', '')}/{db_location}"
    return (namespace, save_dir)


def fake_create_namespace_if_not_exists(
    spark: SparkSession,
    namespace: str = "default",
    append_target: bool = True,
    tenant_name: str | None = None,
    **kwargs,
) -> str:
    """Mock create_namespace_if_not_exists without external calls."""
    current_save_dir = spark.conf.get(SAVE_DIR)
    if not current_save_dir:
        msg = f"Error setting up fixtures: {SAVE_DIR} not set"
        raise ValueError(msg)

    if append_target:
        delta_ns, db_location = gen_ns_save_dir(current_save_dir, namespace, tenant_name)
        spark.sql(f"CREATE DATABASE IF NOT EXISTS {delta_ns} LOCATION '{db_location}'")
        print(f"Namespace {delta_ns} is ready to use at location {db_location}.")
    else:
        delta_ns = namespace
        spark.sql(f"CREATE DATABASE IF NOT EXISTS {delta_ns}")
        print(f"Namespace {delta_ns} is ready to use.")

    return delta_ns


@pytest.fixture
def spark_db(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Generator[tuple[SparkSession, str, str], Any]:
    """Provide a Spark session with a per-test warehouse dir and patched workspace setup."""
    # patch the create_namespace_if_not_exists function
    monkeypatch.setattr(
        "cdm_data_loader_utils.utils.spark_delta.create_namespace_if_not_exists",
        fake_create_namespace_if_not_exists,
    )

    def set_up_test_workspace(*args: str, **kwargs: str | bool) -> tuple[SparkSession, str]:  # noqa: ARG001
        """Local override of set_up_workspace."""
        return original_set_up_ws_fn(*args, local=True, delta_lake=True, override={SAVE_DIR: str(tmp_path)})

    monkeypatch.setattr(spark_delta, "set_up_workspace", set_up_test_workspace)

    _, save_dir = gen_ns_save_dir(str(tmp_path), DEFAULT_NAMESPACE, TENANT_NAME)
    (spark, delta_ns) = spark_delta.set_up_workspace("test_delta_app", DEFAULT_NAMESPACE, TENANT_NAME)

    yield (spark, delta_ns, save_dir)
    spark.stop()


@pytest.mark.parametrize("app_name", [None, "", "my_fave_app"])
def test_get_spark(app_name: str | None, monkeypatch: pytest.MonkeyPatch) -> None:
    """Test of the get_spark utility's ability to fill in the app name if not provided."""

    def fake_get_spark_session(*args: str, **kwargs: str | bool) -> str:  # noqa: ARG001
        if app_name == "my_fave_app":
            assert args[0] == app_name
        else:
            assert args[0] == DEFAULT_APP_NAME
        return "fake spark session"

    monkeypatch.setattr(
        "cdm_data_loader_utils.utils.spark_delta.get_spark_session",
        fake_get_spark_session,
    )

    spark = get_spark(app_name)
    assert spark == "fake spark session"


@pytest.mark.requires_spark
@pytest.mark.parametrize("app_name", [None, "", "my_fave_app"])
def test_get_spark_live(app_name: str | None) -> None:
    """Test of the get_spark utility's ability to fill in the app name if not provided.

    Runs against the live spark engine.
    """
    spark = get_spark(app_name, local=True)
    assert isinstance(spark, SparkSession)
    assert spark.conf.get("spark.app.name") == "my_fave_app" if app_name == "my_fave_app" else DEFAULT_APP_NAME
    spark.stop()


@pytest.mark.parametrize("app_name", [None, "", "my_fave_app"])
@pytest.mark.parametrize("tenant_name", [None, "", "some_tenant"])
@pytest.mark.parametrize("namespace", [None, "", "some_namespace"])
@pytest.mark.parametrize("data_dir", [None, "", "path/to/ws"])
def test_set_up_workspace_defaults(
    app_name: str | None,
    tenant_name: str | None,
    namespace: str | None,
    data_dir: str | None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Check the default values when setting up a workspace."""

    def fake_get_spark_session(*args: str, **kwargs: str | bool) -> str:  # noqa: ARG001
        assert args[0] == app_name if app_name else DEFAULT_APP_NAME
        return "spark session"

    def fake_create_ns(*args: str, **kwargs: str | bool) -> str:
        assert args[0] == "spark session"
        if namespace:
            assert args[1] == namespace
        else:
            assert args[1] == DEFAULT_NAMESPACE
        assert kwargs["tenant_name"] == tenant_name
        return "delta namespace"

    monkeypatch.setattr(
        "cdm_data_loader_utils.utils.spark_delta.get_spark_session",
        fake_get_spark_session,
    )

    monkeypatch.setattr(
        "cdm_data_loader_utils.utils.spark_delta.create_namespace_if_not_exists",
        fake_create_ns,
    )

    if data_dir:
        with pytest.raises(NotImplementedError, match="The data_dir parameter has not been implemented\\."):
            spark_delta.set_up_workspace(app_name, namespace, tenant_name, data_dir)
        return

    spark, delta_ns = spark_delta.set_up_workspace(app_name, namespace, tenant_name, data_dir)
    assert spark == "spark session"
    assert delta_ns == "delta namespace"


@pytest.mark.requires_spark
@pytest.mark.parametrize("tenant_name", [pytest.param(_, id=f"tenant_{_}") for _ in [None, "some_tenant"]])
@pytest.mark.parametrize("namespace", [pytest.param(_, id=f"ns_{_}") for _ in [None, "some_namespace"]])
def test_set_up_workspace_creates_database(
    spark: SparkSession,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    tenant_name: str | None,
    namespace: str | None,
) -> None:
    """Test that setting up a workspace creates the appropriate namespace.

    Mimics the functionality of BERDL's `create_namespace_if_not_exists`.
    """
    app_name = "test_app"
    # expected delta_ns, according to the namespace and tenant_name arguments
    delta_ns = {
        # namespace
        None: {
            # tenant
            None: f"user__some_user__{DEFAULT_NAMESPACE}",
            "some_tenant": f"tenant__some_tenant__{DEFAULT_NAMESPACE}",
        },
        "some_namespace": {
            None: "user__some_user__some_namespace",
            "some_tenant": "tenant__some_tenant__some_namespace",
        },
    }

    expected = delta_ns[namespace][tenant_name]
    assert not spark.catalog.databaseExists(expected)

    def fake_create_namespace_if_not_exists(
        spark: SparkSession,
        namespace: str = "default",
        append_target: bool = True,  # noqa: ARG001
        tenant_name: str | None = None,
    ) -> str:
        """Mock create_namespace_if_not_exists without external calls."""
        namespace = f"tenant__{tenant_name}__{namespace}" if tenant_name else f"user__some_user__{namespace}"
        assert not spark.catalog.databaseExists(namespace)
        spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")
        return namespace

    # patch the create_namespace_if_not_exists function
    monkeypatch.setattr(
        "cdm_data_loader_utils.utils.spark_delta.create_namespace_if_not_exists",
        fake_create_namespace_if_not_exists,
    )

    def set_up_test_workspace(*args: str, **kwargs: str | bool) -> tuple[SparkSession, str]:  # noqa: ARG001
        """Local override of set_up_workspace."""
        return original_set_up_ws_fn(*args, local=True, delta_lake=True, override={SAVE_DIR: str(tmp_path)})

    # patch the set_up_workspace function to add in the various extra kwargs for local use
    monkeypatch.setattr(spark_delta, "set_up_workspace", set_up_test_workspace)

    # create a spark session and ensure that the appropriate database has been created
    (spark, delta_ns) = spark_delta.set_up_workspace(app_name, namespace, tenant_name)
    assert expected == delta_ns
    assert spark.catalog.databaseExists(expected)


@pytest.mark.requires_spark
@pytest.mark.parametrize("dataframe", [None, [1, 2, 3], {}, True])
def test_write_delta_no_data(
    spark: SparkSession,
    dataframe: DataFrame | bool | list | dict | None,
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Ensure that the appropriate message is logged if there is no data to save."""
    if isinstance(dataframe, bool):
        dataframe = spark.createDataFrame([], "name: string, age: int").show()

    output = write_delta(spark, dataframe, "what", "ever", DEFAULT_WRITE_MODE)  # type: ignore
    assert output is None
    assert len(caplog.records) == 1
    for record in caplog.records:
        assert record.levelno == logging.WARNING
        assert record.message == "No data to write to what.ever"


@pytest.mark.requires_spark
@pytest.mark.parametrize("mode", ["some", "mode", 123, None, "whatever"])
def test_write_delta_invalid_write_mode(spark: SparkSession, mode: str, caplog: pytest.LogCaptureFixture) -> None:
    """Ensure that an error is logged if an invalid write mode is supplied."""
    error_msg = f"Invalid mode supplied for writing delta table: {mode}"
    with pytest.raises(ValueError, match=error_msg):
        write_delta(spark, {}, "what", "ever", mode)

    assert len(caplog.records) == 1
    for record in caplog.records:
        assert record.levelno == logging.ERROR
        assert record.message == error_msg


def check_query_output(spark: SparkSession, db_table: str, expected: list[dict[str, Any]]) -> None:
    """Check that the query output matches the expected output."""
    # ensure that the table exists
    assert spark.catalog.tableExists(db_table)
    # run the query
    results = spark.sql(f"SELECT * FROM {db_table}").collect()
    expected_df = spark.createDataFrame(expected)
    assertDataFrameEqual(results, expected_df)


def check_logger_output_successful_write(records: list[logging.LogRecord], db_table: str, mode: str, rows: int) -> None:
    """Check that the logger has emitted the appropriate messages on a successful db write."""
    first_message = records[0]
    assert f"Writing table {db_table} in mode {mode} (rows={rows})" in first_message.message
    assert first_message.levelno == logging.INFO
    last_message = records[-1]
    assert f"Saved managed table {db_table} (rows={rows})" in last_message.message
    assert last_message.levelno == logging.INFO


def check_logger_output_successful_location_write(
    records: list[logging.LogRecord], db_table: str, mode: str, rows: int
) -> None:
    """Check that the logger has emitted the appropriate messages on a successful db write."""
    first_message = records[0]
    assert f"Writing table {db_table} in mode {mode} (rows={rows})" in first_message.message
    assert first_message.levelno == logging.INFO
    last_message = records[-1]
    assert f"Saved external table {db_table} (rows={rows}) to " in last_message.message
    assert last_message.levelno == logging.INFO


def _check_saved_files(parquet_dir: Path) -> None:
    # the directory where table data is expected to be stored
    assert parquet_dir.is_dir()

    # use `sorted` to shortcut the iterator
    parquet_files = sorted(parquet_dir.glob("*.parquet"))
    assert parquet_files

    # the directory where delta table logs are expected to be stored
    delta_log_dir = parquet_dir / "_delta_log"
    assert delta_log_dir.is_dir()
    log_files = sorted(delta_log_dir.glob("*.json"))
    assert log_files


def check_saved_files(ns_save_dir: str | Path, table: str) -> None:
    """Check that the file save operation has saved files in the expected location.

    :param ns_save_dir: save directory for a given namespace
    :type ns_save_dir: str | Path
    :param table: table name
    :type table: str
    """
    parquet_dir = Path(ns_save_dir) / table
    _check_saved_files(parquet_dir)


def populate_db(
    spark: SparkSession, caplog: pytest.LogCaptureFixture, delta_ns: str, table: str, ns_save_dir: str
) -> None:
    """Populate a database, save it as a delta table, and register it with Hive."""
    db_table = f"{delta_ns}.{table}"
    # save a very boring dataframe to a new db_table
    write_delta(
        spark=spark,
        sdf=spark.createDataFrame([DEFAULT_SAMPLE_DATA]),
        delta_ns=delta_ns,
        table=table,
        mode=DEFAULT_WRITE_MODE,
    )
    assert spark.catalog.databaseExists(delta_ns)
    assert spark.catalog.tableExists(db_table)
    # check the db contents are as expected
    check_query_output(spark, db_table, [DEFAULT_SAMPLE_DATA])
    # check there are saved files
    check_saved_files(ns_save_dir, table)
    # check the logger output
    check_logger_output_successful_write(caplog.records, db_table, DEFAULT_WRITE_MODE, 1)


@pytest.mark.requires_spark
@pytest.mark.parametrize("mode", WRITE_MODE)
def test_write_delta_managed_table(
    mode: str,
    spark_db: tuple[SparkSession, str, str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test that a delta table is correctly written and registered in the Hive metastore.

    All valid write modes are tested.
    """
    spark, delta_ns, ns_save_dir = spark_db
    table = f"{mode}_example"
    db_table = f"{delta_ns}.{table}"

    df = spark.createDataFrame([DEFAULT_SAMPLE_DATA])
    write_delta(
        spark=spark,
        sdf=df,
        delta_ns=delta_ns,
        table=table,
        mode=mode,
    )
    check_query_output(spark, db_table, [DEFAULT_SAMPLE_DATA])
    assert len(caplog.records) > 1
    check_logger_output_successful_write(caplog.records, db_table, mode, 1)
    check_saved_files(ns_save_dir, table)


@pytest.mark.requires_spark
def test_write_delta_append_schema_merge(
    spark_db: tuple[SparkSession, str, str], caplog: pytest.LogCaptureFixture
) -> None:
    """Test adding data to an existing db using 'append' mode.

    Append mode should merge schemas, adding new columns without dropping existing data.
    """
    mode = APPEND
    spark, delta_ns, ns_save_dir = spark_db
    table = f"{mode}_test"
    db_table = f"{delta_ns}.{table}"

    populate_db(spark, caplog, delta_ns, table, ns_save_dir)
    caplog.clear()

    # second write - two rows with three columns (SCHEMA CHANGE ALERT!)
    new_rows = [{"a": "A2", "b": "B2", "c": "C2"}, {"a": "A3", "b": "B3", "c": "C3"}]
    write_delta(
        spark=spark,
        sdf=spark.createDataFrame(new_rows),
        delta_ns=delta_ns,
        table=table,
        mode=mode,
    )
    check_saved_files(ns_save_dir, table)
    check_logger_output_successful_write(caplog.records, db_table, mode, 2)
    check_query_output(spark, db_table, [{**DEFAULT_SAMPLE_DATA, "c": None}, *new_rows])


@pytest.mark.requires_spark
def test_write_delta_overwrite_schema(
    spark_db: tuple[SparkSession, str, str], caplog: pytest.LogCaptureFixture
) -> None:
    """Test adding data to an existing db using 'overwrite' mode.

    Overwrite mode should overwrite the original schema and replace any existing data.
    """
    mode = OVERWRITE
    spark, delta_ns, ns_save_dir = spark_db
    table = f"{mode}_test"
    db_table = f"{delta_ns}.{table}"

    populate_db(spark, caplog, delta_ns, table, ns_save_dir)
    caplog.clear()

    # second write - two rows with three columns (SCHEMA CHANGE ALERT!)
    write_delta(
        spark=spark,
        sdf=spark.createDataFrame([Row(x="X2", y="Y2", z="Z2"), Row(x="X3", y="Y3", z="Z3")]),
        delta_ns=delta_ns,
        table=table,
        mode=mode,
    )
    check_saved_files(ns_save_dir, table)
    check_logger_output_successful_write(caplog.records, db_table, mode, 2)
    check_query_output(spark, db_table, [{"x": f"X{n}", "y": f"Y{n}", "z": f"Z{n}"} for n in [2, 3]])


@pytest.mark.requires_spark
@pytest.mark.parametrize("mode", [IGNORE, ERROR, ERROR_IF_EXISTS])
def test_write_delta_ignore_error(
    spark_db: tuple[SparkSession, str, str], caplog: pytest.LogCaptureFixture, mode: str
) -> None:
    """Test adding data to an existing db using 'ignore' or either of the 'error' modes.

    "ignore" will write data if no table exists, but will not write anything if the table already exists.
    "error" and "error_if_exists" would throw an error if the table already exists, but `write_delta` exits early.
    """
    spark, delta_ns, ns_save_dir = spark_db
    table = f"{mode}_test"
    db_table = f"{delta_ns}.{table}"

    populate_db(spark, caplog, delta_ns, table, ns_save_dir)
    caplog.clear()

    # second write - two rows with three columns (SCHEMA CHANGE ALERT!)
    write_delta(
        spark=spark,
        sdf=spark.createDataFrame([Row(x="X2", y="Y2", z="Z2"), Row(x="X3", y="Y3", z="Z3")]),
        delta_ns=delta_ns,
        table=table,
        mode=mode,
    )
    check_saved_files(ns_save_dir, table)
    last_logger_message = caplog.records[-1]
    assert last_logger_message.levelno == logging.WARNING
    assert (
        last_logger_message.message
        == f"Database table {db_table} already exists and writer is set to {mode} mode, so no data would be written. Aborting."
    )
    # check the db contents
    check_query_output(spark, db_table, [DEFAULT_SAMPLE_DATA])


@pytest.mark.requires_spark
def test_write_delta_raise_error(
    spark_db: tuple[SparkSession, str, str], caplog: pytest.LogCaptureFixture, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Ensure that errors are handled gracefully if something terrible happens during saveAsFile."""
    spark, delta_ns, _ = spark_db
    table = "error_handling"
    db_table = f"{delta_ns}.{table}"

    def save_as_oh_crap(*args: str, **kwargs: str | bool) -> None:  # noqa: ARG001
        """Local override of set_up_workspace."""
        msg = "Oh crap!"
        raise RuntimeError(msg)

    monkeypatch.setattr(DataFrameWriter, "saveAsTable", save_as_oh_crap)

    with pytest.raises(Exception, match="Oh crap!"):
        write_delta(
            spark=spark,
            sdf=spark.createDataFrame([Row(x=2, y=3)]),
            delta_ns=delta_ns,
            table=table,
            mode=DEFAULT_WRITE_MODE,
        )
    last_log_record = caplog.records[-1]
    assert last_log_record.levelno == logging.ERROR
    assert last_log_record.message == f"Error writing managed table {db_table}"


@pytest.mark.requires_spark
@pytest.mark.parametrize("mode", WRITE_MODE)
def test_write_delta_uninited_namespace(
    mode: str,
    spark_db: tuple[SparkSession, str, str],
    caplog: pytest.LogCaptureFixture,
) -> None:
    """Test that a namespace that has not been registered throws an error."""
    spark, _delta_ns, _ns_save_dir = spark_db
    table = f"{mode}_example"
    err_msg = "Could not find an appropriate base directory for saving data."
    df = spark.createDataFrame([DEFAULT_SAMPLE_DATA])
    with pytest.raises(RuntimeError, match=err_msg):
        write_delta(
            spark=spark,
            sdf=df,
            delta_ns="namespace_I_just_made_up",
            table=table,
            mode=mode,
        )

    assert caplog.records[-1].levelno == logging.ERROR
    assert caplog.records[-1].message.startswith(err_msg)


@pytest.mark.skip("Not yet implemented")
@pytest.mark.requires_spark
@pytest.mark.parametrize("mode", [APPEND, OVERWRITE])
def test_write_delta_existing_proposed_path_warning(
    mode: str, spark_db: tuple[SparkSession, str, str], caplog: pytest.LogCaptureFixture, tmp_path: Path
) -> None:
    """Test that a warning is emitted if there already exists data saved in another location."""
    spark, delta_ns, _ns_save_dir = spark_db
    table = f"{mode}_example"
    db_table = f"{delta_ns}.{table}"
    err_msg = "Existing path does not match the projected base path for the table. Data written to this directory must be tracked manually."
    save_dir = tmp_path / "save" / "some" / "data" / "here"

    # set up a save directory for the table
    spark.sql(f"CREATE TABLE IF NOT EXISTS {db_table} USING DELTA LOCATION '{save_dir!s}'")

    df = spark.createDataFrame([DEFAULT_SAMPLE_DATA])
    write_delta(
        spark=spark,
        sdf=df,
        delta_ns=delta_ns,
        table=table,
        mode=mode,
    )
    assert caplog.records[0].levelno == logging.WARNING
    assert caplog.records[0].message.startswith(err_msg)


# END write_delta tests. PHEW!


@pytest.mark.requires_spark
def test_preview_or_skip_existing(
    spark_db: tuple[SparkSession, str, str], caplog: pytest.LogCaptureFixture, capsys: pytest.CaptureFixture
) -> None:
    """Test the preview or skip function with an extant db."""
    spark, delta_ns, ns_save_dir = spark_db
    table = "preview_test"
    db_table = f"{delta_ns}.{table}"
    populate_db(spark, caplog, delta_ns, table, ns_save_dir)
    caplog.clear()

    preview_or_skip(spark, delta_ns, table)

    assert caplog.records[0].message == f"Preview for {db_table}:"
    captured = capsys.readouterr().out
    # N.b. this may be fragile if formatting of "show" statements changes
    for k, v in DEFAULT_SAMPLE_DATA.items():
        assert f"|{k} " in captured
        assert f"|{v} " in captured


@pytest.mark.requires_spark
def test_preview_or_skip_missing(spark: SparkSession, caplog: pytest.LogCaptureFixture) -> None:
    """Test the preview or skip function with a missing db."""
    db = "missing"
    table = "not_found"
    preview_or_skip(spark, db, table)

    last_log_message = caplog.records[-1]
    assert last_log_message.message == f"Table {db}.{table} not found. Skipping preview."
