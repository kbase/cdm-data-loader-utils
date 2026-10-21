"""UniRef XML to PyArrow pipeline for the KBase CTS."""

from dlt.common.pipeline import LoadInfo

from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.pipelines.xmltopyarrow.pipeline import run_xmltopyarrow_pipeline
from cdm_data_loaders.pipelines.xmltopyarrow.settings import (
    UNIREF_PIPELINE_NAME as PIPELINE_NAME,
)
from cdm_data_loaders.pipelines.xmltopyarrow.settings import (
    UniRefXmlToPyArrowSettings,
)

__all__ = ["PIPELINE_NAME", "UniRefXmlToPyArrowSettings", "cli", "run_uniref_pipeline"]


def run_uniref_pipeline(settings: UniRefXmlToPyArrowSettings) -> LoadInfo | None:
    """Run the UniRef XML to PyArrow pipeline on matching files in settings.input_dir.

    :param settings: pipeline configuration
    :type  settings: UniRefXmlToPyArrowSettings
    :return: load information for the pipeline
    :rtype: LoadInfo | None
    """
    return run_xmltopyarrow_pipeline(settings)


cli = PipelineCli(UniRefXmlToPyArrowSettings, run_uniref_pipeline)


if __name__ == "__main__":
    cli()
