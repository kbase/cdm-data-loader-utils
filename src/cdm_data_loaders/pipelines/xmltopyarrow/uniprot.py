"""UniProt XML to PyArrow pipeline for the KBase CTS."""

from dlt.common.pipeline import LoadInfo

from cdm_data_loaders.pipelines.core import PipelineCli
from cdm_data_loaders.pipelines.xmltopyarrow.pipeline import run_xmltopyarrow_pipeline
from cdm_data_loaders.pipelines.xmltopyarrow.settings import (
    UNIPROT_PIPELINE_NAME as PIPELINE_NAME,
)
from cdm_data_loaders.pipelines.xmltopyarrow.settings import (
    UniProtXmlToPyArrowSettings,
)

__all__ = ["PIPELINE_NAME", "UniProtXmlToPyArrowSettings", "cli", "run_uniprot_pipeline"]


def run_uniprot_pipeline(settings: UniProtXmlToPyArrowSettings) -> LoadInfo | None:
    """Run the UniProt XML to PyArrow pipeline on matching files in settings.input_dir.

    :param settings: pipeline configuration
    :type  settings: UniProtXmlToPyArrowSettings
    :return: load information for the pipeline
    :rtype: LoadInfo | None
    """
    return run_xmltopyarrow_pipeline(settings)


cli = PipelineCli(UniProtXmlToPyArrowSettings, run_uniprot_pipeline)


if __name__ == "__main__":
    cli()

