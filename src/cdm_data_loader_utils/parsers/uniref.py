"""
UniRef XML Cluster ETL Pipeline

This script downloads a UniRef100 XML file, parses cluster and member information, and writes the extracted data into Delta Lake tables for downstream analysis.

Workflow:
1. Download the specified UniRef100 XML file if not already present locally.
2. Use streaming XML parser to iterate through each <entry>.
3. Extract the following CDM entities:
   - Cluster: Basic metadata for each UniRef cluster.
   - Entity: Tracks unique clusters and creation/update timestamps for idempotent ingest.
   - ClusterMember: Maps all protein members to clusters.
   - CrossReference: Captures additional database cross-references for each UniProt member.
4. Write each parsed table as a Delta Lake table to the specified output directory.
5. Print sample rows for each table after writing for inspection.

Usage:
cd ~/cdm-data-loader-utils

python src/parsers/uniref.py \
  --ftp-url https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.xml.gz \
  --output-dir cdm-data-loader-utils/output/uniref100_clusters \
  --batch-size 1000

python3 uniref.py \
  --ftp-url https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.xml.gz \
  --output-dir output_uniref \
  --batch-size 1000

Parameters:
- --ftp-url:       UniProt FTP URL to the UniRef100 gzipped XML file.
- --output-dir:    Output directory where Delta tables will be written.
- --batch-size:    Number of UniRef entries to process.

"""

import gzip

### ===== logging setup ===== ###
import logging
import os
import uuid
import xml.etree.ElementTree as ET
from datetime import datetime
from urllib.error import URLError
from datetime import timezone
from urllib.request import urlretrieve
import click
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.types import StringType, StructField, StructType
from pathlib import Path

from cdm_data_loader_utils.parsers.xml_utils import get_text, parse_properties


logger = logging.getLogger(__name__)


UNIREF_NS = {"ns": "http://uniprot.org/uniref"}
DATA_SOURCE = "UniRef 100/90/50"


PREFIX_TRANSLATION = {
    "UniProtKB ID": "UniProt",
    "UniProtKB accession": "UniProt",
    "UniParc ID": "UniParc",
    "UniRef90 ID": "UniRef90",
    "UniRef50 ID": "UniRef50",
    "UniRef100 ID": "UniRef100",
}


def generate_dbxref(db: str, acc: str) -> str:
    """Generate a database reference that uses BioRegistry prefixes."""
    return f"{PREFIX_TRANSLATION[db]}:{acc}"


def cdm_entity_id(value: str, prefix: str = "CDM:") -> str:
    """
    Deterministic UUIDv5-based CDM id generator.

    value must be non-empty.
    """
    if not value:
        raise ValueError("Value must be a non-empty string")

    return f"{prefix}{uuid.uuid5(uuid.NAMESPACE_OID, value)}"


# timestamp helper
def get_timestamps(
    uniref_id: str,
    existing_created: dict[str, str],
    now: datetime | None = None,
) -> tuple[str, str]:
    """
    Return (updated_time, created_time) for a given UniRef cluster ID.
    - All timestamps are UTC ISO8601 with timezone (e.g., 2026-01-05T12:34:56+00:00)
    - uniref_id must be non-empty (schema invariant)
    """
    if not uniref_id:
        raise ValueError("get_timestamps: uniref_id must be a non-empty string")

    now_dt = now or datetime.now(timezone.utc)
    updated_time = now_dt.isoformat(timespec="seconds")

    created_time = existing_created.get(uniref_id) or updated_time
    return updated_time, created_time


def download_file(url: str, local_path: str, overwrite: bool = False) -> str:
    """
    Download URL -> local_path.
    - Atomic: downloads to .part then os.replace
    - Idempotent: skips if exists unless overwrite=True
    Returns the final local_path.
    """
    dst = Path(local_path)
    dst.parent.mkdir(parents=True, exist_ok=True)

    if dst.exists() and not overwrite:
        logger.info("File already exists, skip download: %s", dst)
        return str(dst)

    tmp = dst.with_suffix(dst.suffix + ".part")

    try:
        if tmp.exists():
            tmp.unlink()
    except Exception:
        logger.exception("Failed to remove partial download: %s", tmp)

    logger.info("Downloading %s -> %s", url, dst)
    try:
        urlretrieve(url, str(tmp))
        os.replace(tmp, dst)
        logger.info("Download complete: %s", dst)
        return str(dst)
    except Exception:
        logger.exception("Failed to download %s", url)
        try:
            if tmp.exists():
                tmp.unlink()
        except Exception:
            logger.exception("Failed to cleanup tmp file: %s", tmp)
        raise


def load_existing_created(spark: SparkSession, entity_table: str | None) -> dict[str, str]:
    """
    Load mapping data_source_entity_id -> created timestamp from the Entity Delta table.
    Returns an empty dict if the table does not exist.
    """
    existing_created: dict[str, str] = {}
    if not entity_table:
        logger.warning("Entity table path not specified.")
        return existing_created

    try:
        df = spark.read.format("delta").load(entity_table).select("data_source_entity_id", "created")
        existing_created = {row["data_source_entity_id"]: row["created"] for row in df.collect()}
        logger.info(f"Loaded {len(existing_created)} existing created timestamps from {entity_table}.")
    except Exception as e:
        logger.warning(f"No existing Delta table found at {entity_table}. Starting fresh. ({e.__class__.__name__})")

    return existing_created


##### -------------- List utility function --------------- #####
def extract_cluster(
    elem: ET.Element,
    ns: dict[str, str],
    uniref_id: str,
) -> tuple[str, str]:
    """
    Extract a deterministic CDM cluster_id and the UniRef cluster name.
    """
    cluster_id = cdm_entity_id(value=uniref_id) or f"CDM:{uuid.uuid4()}"
    name_elem = elem.find("ns:name", ns)
    name = get_text(elem=name_elem, default="UNKNOWN") or "UNKNOWN"

    return cluster_id, name


def get_accession_and_seed(dbref: ET.Element | None, ns: dict[str, str]) -> tuple[str | None, bool]:
    """
    Extract UniProtKB accession and is_seed status from a dbReference element.
    """

    if dbref is None:
        return None, False

    props = parse_properties(dbref, ns)

    raw_acc = props.get("UniProtKB accession")
    if isinstance(raw_acc, list):
        accession = raw_acc[0] if raw_acc else None
    else:
        accession = raw_acc  # string or None

    raw_seed = props.get("isSeed")
    if isinstance(raw_seed, list):
        is_seed = bool(raw_seed) and raw_seed[0].lower() == "true"
    else:
        is_seed = raw_seed is not None and raw_seed.lower() == "true"

    return accession, is_seed


def add_cluster_members(
    cluster_id: str,
    repr_db: ET.Element | None,
    elem: ET.Element,
    cluster_member_rows: list[tuple[str, str, str, str, str]],
    ns: dict[str, str],
) -> None:
    """Populate cluster_member_rows with representative, member records."""
    dbrefs: list[tuple[ET.Element, bool]] = []
    if repr_db is not None:
        dbrefs.append((repr_db, True))
    for mem in elem.findall("ns:member/ns:dbReference", ns):
        dbrefs.append((mem, False))

    for dbref, is_representative in dbrefs:
        accession, is_seed = get_accession_and_seed(dbref, ns)
        if not accession:
            continue

        member_entity_id = cdm_entity_id(accession)
        if not member_entity_id:
            continue

        cluster_member_rows.append((
            cluster_id,
            member_entity_id,
            str(is_representative).lower(),
            str(is_seed).lower(),
            "1.0",  # score placeholder
        ))


def extract_cross_refs(dbref, cross_reference_data, ns) -> None:
    if dbref is None:
        return
    props = {p.attrib["type"]: p.attrib["value"] for p in dbref.findall("ns:property", ns)}
    entity_id = cdm_entity_id(dbref.attrib.get("id"))
    xref_types = ["UniRef90 ID", "UniRef50 ID", "UniParc ID"]
    for i in xref_types:
        if i in props:
            cross_reference_data.append((entity_id, i, props[i]))


def parse_uniref_entry(
    elem: ET.Element, existing_created: dict[str, str], ns: dict[str, str]
) -> dict[str, list[tuple]]:
    """
    Parse a single UniRef <entry> element into CDM-friendly row tuples.
    """
    cluster_rows: list[tuple[str, str, str, str | None, str]] = []
    entity_rows: list[tuple[str, str, str, str, str, str]] = []
    member_rows: list[tuple[str, str, str, str, str]] = []
    xref_rows: list[tuple[str, str, str]] = []

    uniref_id = elem.attrib.get("id") or ""

    cluster_id, name = extract_cluster(elem, ns, uniref_id)
    updated_time, created_time = get_timestamps(uniref_id, existing_created)

    # Cluster table
    cluster_rows.append((
        cluster_id,
        name,
        "protein",
        None,
        DATA_SOURCE,
    ))

    # Entity table
    entity_rows.append((
        cluster_id,
        uniref_id,
        "Cluster",
        DATA_SOURCE,
        updated_time,
        created_time,
    ))

    # Cross references from representative and members
    repr_db = elem.find("ns:representativeMember/ns:dbReference", ns)
    extract_cross_refs(repr_db, xref_rows, ns)

    for mem in elem.findall("ns:member/ns:dbReference", ns):
        extract_cross_refs(mem, xref_rows, ns)

    # Cluster members (representative + members)
    add_cluster_members(cluster_id, repr_db, elem, member_rows, ns)

    return {
        "cluster_data": cluster_rows,
        "entity_data": entity_rows,
        "cluster_member_data": member_rows,
        "cross_reference_data": xref_rows,
    }


##### -------------- Parse Uniref XML --------------- #####
def parse_uniref_xml(local_gz: str, batch_size: int, existing_created: dict[str, str]) -> dict[str, list[tuple]]:
    """
    Stream-parse UniRef XML (gzipped) and extract CDM-like row tuples.
    """
    ns = UNIREF_NS
    entry_count = 0

    cluster_data: list[tuple] = []
    entity_data: list[tuple] = []
    cluster_member_data: list[tuple] = []
    cross_reference_data: list[tuple] = []

    with gzip.open(local_gz, "rb") as f:
        context = ET.iterparse(f, events=("end",))
        for _, elem in context:
            if not elem.tag.endswith("entry"):
                continue

            parsed = parse_uniref_entry(elem, existing_created, ns)
            cluster_data.extend(parsed["cluster_data"])
            entity_data.extend(parsed["entity_data"])
            cluster_member_data.extend(parsed["cluster_member_data"])
            cross_reference_data.extend(parsed["cross_reference_data"])

            entry_count += 1
            if entry_count >= batch_size:
                break

            elem.clear()

    logger.info(f"Parsed {entry_count} clusters")
    return {
        "cluster_data": cluster_data,
        "entity_data": entity_data,
        "cluster_member_data": cluster_member_data,
        "cross_reference_data": cross_reference_data,
    }


##### -------------- Save delta table and print the preview --------------- #####
def save_delta_tables(spark, output_dir, data_dict):
    # Cluster
    cluster_schema = StructType([
        StructField("cluster_id", StringType(), False),
        StructField("name", StringType(), False),
        StructField("entity_type", StringType(), False),
        StructField("description", StringType(), True),
        StructField("protocol_id", StringType(), False),
    ])

    cluster_df = spark.createDataFrame(data_dict["cluster_data"], cluster_schema)
    cluster_df.write.format("delta").mode("overwrite").save(os.path.join(output_dir, "Cluster"))
    logger.info(f"Cluster Delta table written to: {os.path.join(output_dir, 'Cluster')}")

    # Entity
    entity_schema = StructType([
        StructField("entity_id", StringType(), False),
        StructField("data_source_entity_id", StringType(), False),
        StructField("entity_type", StringType(), False),
        StructField("data_source", StringType(), False),
        StructField("updated", StringType(), False),
        StructField("created", StringType(), False),
    ])

    entity_df = spark.createDataFrame(data_dict["entity_data"], entity_schema)
    entity_table_path = os.path.join(output_dir, "Entity")
    entity_df.write.format("delta").mode("overwrite").save(entity_table_path)
    logger.info(f"Entity Delta table written to: {entity_table_path}")

    # ClusterMember
    cluster_member_schema = StructType([
        StructField("cluster_id", StringType(), False),
        StructField("entity_id", StringType(), False),
        StructField("is_representative", StringType(), False),
        StructField("is_seed", StringType(), False),
        StructField("score", StringType(), False),
    ])

    cluster_member_df = spark.createDataFrame(data_dict["cluster_member_data"], cluster_member_schema)
    cluster_member_path = os.path.join(output_dir, "ClusterMember")
    cluster_member_df.write.format("delta").mode("overwrite").save(cluster_member_path)
    logger.info(f"ClusterMember Delta table written to: {cluster_member_path}")

    # CrossReference
    cross_reference_schema = StructType([
        StructField("entity_id", StringType(), False),
        StructField("xref_type", StringType(), False),
        StructField("xref_value", StringType(), False),
    ])

    cross_reference_df = spark.createDataFrame(data_dict["cross_reference_data"], cross_reference_schema)
    cross_reference_path = os.path.join(output_dir, "CrossReference")
    cross_reference_df.write.format("delta").mode("overwrite").save(cross_reference_path)
    logger.info(f"CrossReference Delta table written to: {cross_reference_path}")

    # Previews
    logger.info("Sample Clusters:")
    cluster_df.createOrReplaceTempView("Cluster")
    spark.sql("SELECT * FROM Cluster LIMIT 20").show(truncate=False)

    logger.info("Sample Entities:")
    entity_df.createOrReplaceTempView("Entity")
    spark.sql("SELECT * FROM Entity LIMIT 20").show(truncate=False)

    logger.info("Sample ClusterMembers:")
    cluster_member_df.createOrReplaceTempView("ClusterMember")
    spark.sql("SELECT * FROM ClusterMember LIMIT 20").show(truncate=False)

    logger.info("Sample CrossReferences:")
    cross_reference_df.createOrReplaceTempView("CrossReference")
    spark.sql("SELECT * FROM CrossReference LIMIT 20").show(truncate=False)


def build_spark_session():
    # Build a Spark session configured for Delta Lake support
    builder = (
        SparkSession.builder.appName("UniRef Cluster Extractor")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config(
            "spark.sql.catalog.spark_catalog",
            "org.apache.spark.sql.delta.catalog.DeltaCatalog",
        )
    )
    return configure_spark_with_delta_pip(builder).getOrCreate()


@click.command()
@click.option("--ftp-url", required=True, help="FTP URL to UniRef100 XML file")
@click.option("--output-dir", required=True, help="Output directory for Delta table")
@click.option("--batch-size", default=1000, help="Number of UniRef entries to parse (limit)")
def main(ftp_url, output_dir, batch_size):
    # set up logging in CLI context
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] (%(name)s:%(lineno)d %(message)s",
    )

    logger.info("Starting UniRef100/90/50 Import Pipeline")

    # Set local path for downloaded gzipped XML file
    local_gz = os.path.join("/tmp", os.path.basename(ftp_url))

    # Download UniRef XML if not already present
    try:
        download_file(ftp_url, local_gz)
    except URLError as e:
        logger.error(f"Error! Cannot download file: {e.reason}")
        return

    # Start Spark session with Delta Lake support
    logger.info("Building Spark session:")
    spark = build_spark_session()

    # Load existing entity creation timestamps
    try:
        entity_table_path = os.path.join(output_dir, "Entity")
        existing_created = load_existing_created(spark, entity_table_path)

        # Parse the UniRef XML and extract all CDM table data
        logger.info("Parsing UniRef XML:")
        data_dict = parse_uniref_xml(local_gz, batch_size, existing_created)

        # Write parsed data to Delta tables in output directory
        logger.info("Saving Delta tables:")
        save_delta_tables(spark, output_dir, data_dict)

        logger.info("UniRef100/90/50 Import Pipeline completed successfully.")

    finally:
        spark.stop()
        logger.info("Spark session stopped.")


if __name__ == "__main__":
    main()
