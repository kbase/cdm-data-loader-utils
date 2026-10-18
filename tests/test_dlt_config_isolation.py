"""Regression tests for the dlt config isolation helpers in tests/dlt_config_isolation.py."""

from pathlib import Path

import dlt
from dlt.common.configuration.container import Container
from dlt.common.pipeline import PipelineContext

from tests.conftest import LOCAL_FS, S3, _generate_dlt_config
from tests.dlt_config_isolation import deactivated_pipeline, dlt_config_unset, isolated_dlt_config

SEED_CONFIG_KEY: str = "destination.local_fs.bucket_url"


def test_isolated_dlt_config_scopes_seeded_values() -> None:
    """A nested isolated_dlt_config resolves its own seed and restores the outer chain on exit."""
    seed_bucket_url: str = _generate_dlt_config()[SEED_CONFIG_KEY]
    with isolated_dlt_config({"destination": {LOCAL_FS: {"bucket_url": "/isolated/nested"}}}):
        assert dlt.config.get(SEED_CONFIG_KEY) == "/isolated/nested"
    assert dlt.config.get(SEED_CONFIG_KEY) == seed_bucket_url


def test_isolated_dlt_config_restores_chain_after_direct_mutation() -> None:
    """A raw dlt.config[key] = value inside an isolated scope cannot leak past the scope."""
    seed_bucket_url: str = _generate_dlt_config()[SEED_CONFIG_KEY]
    with isolated_dlt_config():
        dlt.config[SEED_CONFIG_KEY] = "/isolated/leaky"
        assert dlt.config.get(SEED_CONFIG_KEY) == "/isolated/leaky"
    assert dlt.config.get(SEED_CONFIG_KEY) == seed_bucket_url


def test_isolated_dlt_config_seeds_and_restores_secrets() -> None:
    """The isolated chain resolves seeded secrets and drops them on exit."""
    with isolated_dlt_config(secrets={"destination": {S3: {"credentials": {"access_key_id": "test-key"}}}}):
        assert dlt.secrets.get("destination.s3.credentials.access_key_id") == "test-key"
    assert dlt.secrets.get("destination.s3.credentials.access_key_id") is None


def test_dlt_config_unset_makes_dlt_config_none() -> None:
    """dlt_config_unset simulates dlt.config being entirely absent, restoring it on exit."""
    with dlt_config_unset():
        assert dlt.config is None
    assert dlt.config is not None


def test_deactivated_pipeline_deactivates_pipeline_created_inside(tmp_path: Path) -> None:
    """A dlt.pipeline created inside the block is deactivated when the block exits."""
    with deactivated_pipeline():
        dlt.pipeline(pipeline_name="deactivated_pipeline_test", pipelines_dir=str(tmp_path / "pipelines"))
        assert Container()[PipelineContext].is_active()
    assert not Container()[PipelineContext].is_active()


def test_autouse_isolation_scopes_direct_mutation() -> None:
    """A raw dlt.config mutation in a test body is scoped to that test by the autouse fixture."""
    seed_bucket_url: str = _generate_dlt_config()[SEED_CONFIG_KEY]
    assert dlt.config.get(SEED_CONFIG_KEY) == seed_bucket_url
    dlt.config[SEED_CONFIG_KEY] = "/isolated/leaky"
    assert dlt.config.get(SEED_CONFIG_KEY) == "/isolated/leaky"


def test_autouse_isolation_chain_is_seed_clean_between_tests() -> None:
    """The chain a test starts with matches the autouse seed, proving prior mutations did not leak."""
    seed_bucket_url: str = _generate_dlt_config()[SEED_CONFIG_KEY]
    assert dlt.config.get(SEED_CONFIG_KEY) == seed_bucket_url
