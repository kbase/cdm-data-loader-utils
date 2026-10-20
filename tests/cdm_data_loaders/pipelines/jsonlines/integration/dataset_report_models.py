"""Entity registry for the NCBI assembly-report reference dataset."""

from pydantic import BaseModel

from cdm_data_loaders.parsers.ncbi.datasets_api.openapi_pydantic import V2reportsAssemblyDataReport

ENTITY_MODELS: dict[str, type[BaseModel]] = {"dataset": V2reportsAssemblyDataReport}
