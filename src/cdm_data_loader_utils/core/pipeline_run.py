"""Pipeline run class."""

from dataclasses import dataclass


@dataclass(frozen=True)
class PipelineRun:
    """Dataclass for capturing ingestion run information."""

    run_id: str
    pipeline: str
    source_path: str
    namespace: str
