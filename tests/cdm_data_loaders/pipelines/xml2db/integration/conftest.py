"""Shared fixtures for xml2db pipeline integration tests."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest
from dlt.common.pipeline import LoadInfo
from frozendict import frozendict

from cdm_data_loaders.pipelines.xml2db.pipeline import run_xml2db_pipeline
from cdm_data_loaders.pipelines.xml2db.settings import Xml2DbSettings
from tests.conftest import REFERENCE_XSD
from tests.helpers import make_cts_settings
from tests.xml_samples import LIBRARY_XSD

DEFAULT_XML2DB_SETTINGS: frozendict = frozendict(
    {
        "buffer_size": 10,
        "dlt_dev_mode": False,
        "file_glob": "*.xml*",
        "log_interval": 1000,
        "use_output_dir_for_pipeline_metadata": True,
    }
)


@pytest.fixture
def library_xsd(tmp_path: Path) -> Path:
    """Write the small library.xsd schema used for basic pipeline-mechanics tests."""
    xsd_path = tmp_path / "library.xsd"
    xsd_path.write_text(LIBRARY_XSD, encoding="utf-8")
    return xsd_path


@pytest.fixture(scope="session")
def reference_xsd() -> Path:
    """Path to the uniref.xsd schema."""
    return REFERENCE_XSD


@pytest.fixture
def settings_factory(tmp_path: Path, dlt_destination_config: str, library_xsd: Path) -> Callable[..., Xml2DbSettings]:
    """Return a factory that builds a valid Xml2DbSettings using the library.xsd schema."""

    def _factory(**overrides: Any) -> Xml2DbSettings:  # noqa: ANN401
        defaults: dict[str, Any] = {
            **DEFAULT_XML2DB_SETTINGS,
            "dataset_name": "xml2db_test_dataset",
            "output_dir": str(tmp_path / "bucket"),
            "use_destination": dlt_destination_config,
            "xsd_file": str(library_xsd),
        }
        defaults.update(overrides)
        return make_cts_settings(Xml2DbSettings, tmp_path, **defaults)

    return _factory


@pytest.fixture
def run_xml2db_pipeline_fxt(
    reference_xml_data_dir: Path,
    reference_xsd: Path,
    dlt_destination_config: str,
    isolated_run_factory: Callable[..., Callable[..., tuple[LoadInfo | None, Path]]],
) -> Callable[..., Any]:
    """Return a factory that runs the xml2db pipeline on a chunk dir with full isolation.

    Each run gets a unique output_dir and dataset_name, so parallel-safe
    repeated runs never collide or accumulate rows. `isolated_run_factory` is called once per
    test, not once per chunk_dir/override combination, so its output_dir counter is shared and
    keeps incrementing across every call this fixture's factory makes within one test -- calling
    it again per call would reset the counter and collide on the same output_dir every time.
    """
    run = isolated_run_factory(
        "reference_run",
        run_xml2db_pipeline,
        {
            "buffer_size": 10,
            "dlt_dev_mode": False,
            "file_glob": "*.xml*",
            "log_interval": 1000,
            "short_name": "uniref",
            "settings_cls": Xml2DbSettings,
            "use_destination": dlt_destination_config,
            "use_output_dir_for_pipeline_metadata": True,
            "xsd_file": str(reference_xsd),
        },
    )

    def _factory(chunk_dir: str, **overrides: Any) -> tuple[LoadInfo | None, Path]:  # noqa: ANN401
        return run(input_dir=str(reference_xml_data_dir / chunk_dir), **overrides)

    return _factory
