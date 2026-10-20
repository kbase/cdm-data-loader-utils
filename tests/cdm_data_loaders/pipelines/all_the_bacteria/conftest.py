from typing import Any

import pytest
from frozendict import frozendict

from cdm_data_loaders.core.fields import INPUT_DIR, OUTPUT_DIR, USE_DESTINATION
from cdm_data_loaders.pipelines.all_the_bacteria import PROJECT_PART_REGEX, AtbSettings, load_patterns
from tests.conftest import S3


@pytest.fixture(params=["default", "v1"])
def atb_settings_pair(
    request: pytest.FixtureRequest,
    test_cts_settings: frozendict[str, Any],
    test_cts_settings_reconciled: frozendict[str, Any],
) -> tuple[frozendict[str, Any], frozendict[str, Any]]:
    """Return (settings, reconciled) dicts for AtbSettings variants."""
    output_dir = test_cts_settings_reconciled[OUTPUT_DIR]
    if request.param == "default":
        settings = frozendict(**test_cts_settings, version="1.2.3", pattern_file=None)
        reconciled = frozendict(
            **test_cts_settings_reconciled,
            version="1.2.3",
            pattern_file=None,
            raw_data_dir=output_dir / "raw_data" / "1.2.3",
            pattern_matches=PROJECT_PART_REGEX,
        )
        return settings, reconciled

    pattern_rel = "some/path/or/other.txt"
    pattern_path = test_cts_settings_reconciled[INPUT_DIR] / pattern_rel
    pattern_path.parent.mkdir(parents=True, exist_ok=True)
    pattern_path.write_text("hello world\n")
    settings = frozendict(**test_cts_settings, version=1, pattern_file=pattern_rel)
    reconciled = frozendict(
        **test_cts_settings_reconciled,
        version="1",
        pattern_file=pattern_rel,
        raw_data_dir=output_dir / "raw_data" / "1",
        pattern_matches=load_patterns(pattern_path),
    )
    return settings, reconciled


@pytest.fixture
def test_settings(test_defaults: dict[str, Any]) -> AtbSettings:
    """Generate fake settings for testing."""
    return AtbSettings(**test_defaults)


@pytest.fixture
def test_s3_settings(test_defaults: dict[str, Any]) -> AtbSettings:
    """Generate fake settings that use s3."""
    return AtbSettings(**{**test_defaults, USE_DESTINATION: S3})
