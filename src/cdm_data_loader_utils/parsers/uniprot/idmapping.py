"""Parser for UniProt ID Mapping file.

UniProt provides comprehensive mappings between their protein IDs and many other databases. Mappings are extracted from the UniProt records where possible.

Source file: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz

Legacy mappings (pre-UniProt proteome redundancy reduction drive): https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.2015_03.gz

Docs: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/README

Metalink file: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/RELEASE.metalink

Retrieve the list of databases referenced from the UniProt API: https://rest.uniprot.org/database/stream?format=json&query=%28*%29

1) idmapping.dat
This file has three columns, delimited by tab:
1. UniProtKB-AC
2. ID_type
3. ID
where ID_type is the database name as appearing in UniProtKB cross-references,
and as supported by the ID mapping tool on the UniProt web site,
http://www.uniprot.org/mapping and where ID is the identifier in
that cross-referenced database.
"""

import datetime
from uuid import uuid4

import click
from pyspark.sql import DataFrame, SparkSession
from pyspark.sql import functions as sf
from pyspark.sql.types import StringType, StructField

from cdm_data_loader_utils.core.constants import INVALID_DATA_FIELD_NAME
from cdm_data_loader_utils.core.pipeline_run import PipelineRun
from cdm_data_loader_utils.readers.dsv import read
from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger
from cdm_data_loader_utils.utils.minio import list_remote_dir_contents
from cdm_data_loader_utils.utils.spark_delta import APPEND, set_up_workspace, write_delta
from cdm_data_loader_utils.validation.dataframe_validator import DataFrameValidator, Validator
from cdm_data_loader_utils.validation.df_nullable_fields import validate as check_nullable_fields

APP_NAME = "uniprot_idmapping"
NOW = datetime.datetime.now(tz=datetime.UTC)
DB = "db"
XREF = "xref"
ID = "id"
COLUMNS = [ID, DB, XREF]


logger = get_cdm_logger()

ID_MAPPING_SCHEMA = [StructField(n, StringType(), nullable=False) for n in COLUMNS]


def ingest(spark: SparkSession, run: PipelineRun, id_mapping_tsv: str) -> DataFrame:
    """Parse the ID mapping file and convert it to a dataframe.

    :param spark: spark sesh
    :type spark: SparkSession
    :param id_mapping_tsv: path to the ID mapping tsv file
    :type id_mapping_tsv: str
    :return: dataframe containing the ID mapping stuff
    :rtype: DataFrame
    """
    options = {
        "delimiter": "\t",
        "header": False,
        "ignoreLeadingWhiteSpace": True,
        "ignoreTrailingWhiteSpace": True,
        "enforceSchema": True,
        "inferSchema": False,
    }

    df = read(spark, id_mapping_tsv, ID_MAPPING_SCHEMA, options)
    id_map_parse_result = DataFrameValidator(spark).validate_dataframe(
        data_to_validate=df,
        schema=ID_MAPPING_SCHEMA,
        run=run,
        validator=Validator(check_nullable_fields, {"invalid_col": INVALID_DATA_FIELD_NAME}),
        invalid_col=INVALID_DATA_FIELD_NAME,
    )
    id_map_df = id_map_parse_result.valid_df

    # destination format:
    # "entity_id", "identifier", "description", "source", "relationship"
    return id_map_df.select(
        # prefix with UniProt
        sf.concat(sf.lit("UniProt:"), sf.col("id")).alias("uniprot_id"),
        sf.col(DB),
        sf.col(XREF),
        sf.lit(None).cast("string").alias("description"),
        sf.lit("UniProt ID mapping").alias("source"),
        sf.lit(None).cast("string").alias("relationship"),
    )


def read_and_write(spark: SparkSession, pipeline_run: PipelineRun, id_mapping_tsv: str) -> None:
    """Read in the UniProt ID mapping and write it out as a uniprot_identifier table.

    :param spark: spark sesh
    :type spark: SparkSession
    :param delta_ns: namespace to write to
    :type delta_ns: str
    :param id_mapping_tsv: path to the ID mapping file
    :type id_mapping_tsv: str
    :param mode: write mode (append or overwrite)
    :type mode: str
    """
    # get the metalink XML and retrieve data source info
    write_delta(
        spark, ingest(spark, pipeline_run, id_mapping_tsv), pipeline_run.namespace, "uniprot_identifier", APPEND
    )


@click.command()
@click.option(
    "--source",
    required=True,
    help="Full path to the source directory containing ID mapping file(s). S3 buckets should include the full path with s3a:// prefix; otherwise, the path will be assumed to be local.",
)
@click.option(
    "--namespace",
    default="uniprot",
    show_default=True,
    help="Delta Lake database name",
)
@click.option(
    "--tenant-name",
    default=None,
    help="Tenant warehouse to save processed data to; defaults to saving data to the user warehouse if a tenant is not specified",
)
def cli(source: str, namespace: str, tenant_name: str | None) -> None:
    """Run the UniProt ID Mapping importer.

    :param source: full path to the source directory containing ID mapping file(s)
    :type source: str
    :param namespace: Delta Lake database name
    :type namespace: str
    :param tenant_name: Tenant warehouse to save processed data to; defaults to saving data to the user warehouse if a tenant is not specified
    :type tenant_name: str | None
    """
    (spark, delta_ns) = set_up_workspace(APP_NAME, namespace, tenant_name)
    bucket_list = []
    if "://" in source and source.startswith("s3a://"):
        # we're golden
        bucket_list = list_remote_dir_contents(source)
    # TODO: other locations

    for file in bucket_list:
        # file names are in the 'Key' value
        # 'tenant-general-warehouse/kbase/datasets/uniprot/id_mapping/id_mapping_part_001.tsv.gz'
        pipeline_run = PipelineRun(str(uuid4()), APP_NAME, file["Key"], delta_ns)
        read_and_write(spark, pipeline_run, file["Key"])


if __name__ == "__main__":
    cli()
