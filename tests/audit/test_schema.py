"""Tests for the schema package in the audit section."""

import pytest

from cdm_data_loader_utils.audit.schema import PIPELINE, RUN_ID, SOURCE, current_run_expr, match_run
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from tests.audit.conftest import PIPELINE_RUN


def test_current_run_expr_no_args() -> None:
    """Test generation of the current run expression without any arguments."""
    assert current_run_expr() == f"t.{RUN_ID} = s.{RUN_ID} AND t.{SOURCE} = s.{SOURCE} AND t.{PIPELINE} = s.{PIPELINE}"


@pytest.mark.parametrize("target", ["", "some_target"])
@pytest.mark.parametrize("source", ["", "some_source"])
def test_current_run_expr(target: str, source: str) -> None:
    """Test generation of the current run expression with some arguments."""
    t_str = target if target else "t"
    s_str = source if source else "s"
    assert (
        current_run_expr(target, source)
        == f"{t_str}.{RUN_ID} = {s_str}.{RUN_ID} AND {t_str}.{SOURCE} = {s_str}.{SOURCE} AND {t_str}.{PIPELINE} = {s_str}.{PIPELINE}"
    )


def test_match_run(pipeline_run: PipelineRun) -> None:
    """Test generation of the current run expression without any arguments."""
    assert (
        match_run(pipeline_run)
        == f"{RUN_ID} = '{PIPELINE_RUN[RUN_ID]}' AND {PIPELINE} = '{PIPELINE_RUN[PIPELINE]}' AND {SOURCE} = '{PIPELINE_RUN[SOURCE]}'"
    )
