"""
UniRef XML Cluster ETL Pipeline.

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

**Parameters:**
- --ftp-url:       UniProt FTP URL to the UniRef100 gzipped XML file.
- --output-dir:    Output directory where Delta tables will be written.
- --batch-size:    Number of UniRef entries to process.

"""

import gzip
import os
import uuid
import xml.etree.ElementTree as ET
from datetime import datetime
from urllib.request import URLError, urlretrieve

import click
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.types import StringType, StructField, StructType


# Generate a unique CDM entity_id based on accession
def cdm_entity_id(accession) -> str | None:
    if not accession:
        return None
    uuid_part = uuid.uuid5(uuid.NAMESPACE_OID, accession)
    return f"CDM:{uuid_part}"


# Download a file from the specified URL to the local path if it does not already exist
def download_file(url, local_path) -> None:
    """
    If the file is already present at local, the function does nothing.
    If the download fails, any partially downloaded file will be removed.
    """
    if not os.path.exists(local_path):
        print(f"Downloading from URL link: {url}")
        try:
            urlretrieve(url, local_path)
            print("Download completed!")
        except Exception as e:
            print(f"Failed to download {url}: {e}")
            if os.path.exists(local_path):
                os.remove(local_path)
            raise
    else:
        print(f"File already exists: {local_path}")


# Load mapping from data_source_entity_id to created timestamp from Delta table
def load_existing_created(spark, entity_table):
    existing_created = {}
    if not entity_table:
        print("Entity table path not specified.")
        return existing_created

    try:
        df = spark.read.format("delta").load(entity_table).select("data_source_entity_id", "created")
        existing_created = {row["data_source_entity_id"]: row["created"] for row in df.collect()}
        print(f"Loaded {len(existing_created)} existing created timestamps.")
    except Exception as e:
        print(f"No existing Delta table found at {entity_table}. Starting fresh. ({e.__class__.__name__})")

    return existing_created


##### -------------- List utility function --------------- #####


# Helper function to extract basic cluster info from XML entry element
def extract_cluster(elem, ns):
    cluster_id = f"CDM:{uuid.uuid4()}"
    name_elem = elem.find("ns:name", ns)
    name = name_elem.text if name_elem is not None else "UNKNOWN"
    return cluster_id, name


# Returns tuple of (updated_time, created_time)
def get_timestamps(uniref_id, existing_created, now=None):
    now_dt = now or datetime.now()
    formatted_now = now_dt.strftime("%Y-%m-%dT%H:%M:%S")
    created = existing_created.get(uniref_id)
    created_time = (created.split(".")[0] if "." in created else created) if created else formatted_now
    return formatted_now, created_time


# Extract UniProtKB accession and is_seed status from a dbReference element
def get_accession_and_seed(dbref, ns):
    if dbref is None:
        return None, False
    prop_elems = dbref.findall("ns:property", ns)

    props = {}
    for prop in prop_elems:
        t = prop.attrib["type"]
        v = prop.attrib["value"]
        props[t] = v

    acc = props.get("UniProtKB accession") or dbref.attrib.get("id")
    is_seed = props.get("isSeed", "false").lower() == "true"
    return acc, is_seed


# Add both representative and other cluster members into cluster_member_data list
def add_cluster_members(cluster_id, repr_db, elem, cluster_member_data, ns) -> None:
    dbrefs = []
    if repr_db is not None:
        dbrefs.append((repr_db, True))
    for mem in elem.findall("ns:member/ns:dbReference", ns):
        dbrefs.append((mem, False))

    for dbref, is_representative in dbrefs:
        acc, is_seed = get_accession_and_seed(dbref, ns)
        if acc:
            member_entity_id = cdm_entity_id(acc)
            cluster_member_data.append(
                (cluster_id, member_entity_id, str(is_representative).lower(), str(is_seed).lower(), "1.0")
            )


# Extract cross-references (UniRef90/50/UniParc) from a dbReference element
def extract_cross_refs(dbref, cross_reference_data, ns) -> None:
    if dbref is None:
        return
    props = {p.attrib["type"]: p.attrib["value"] for p in dbref.findall("ns:property", ns)}
    entity_id = cdm_entity_id(dbref.attrib.get("id"))
    xref_types = ["UniRef90 ID", "UniRef50 ID", "UniParc ID"]
    for i in xref_types:
        if i in props:
            cross_reference_data.append((entity_id, i, props[i]))


##### -------------- Parse Uniref XML --------------- #####


def parse_uniref_xml(local_gz, batch_size, existing_created):
    """
    Parse UniRef XML (gzipped) and extract cluster, entity, cluster member, UniProtKB member, and cross-reference info.

    Args:
        local_gz (str): Local gzipped UniRef XML path.
        batch_size (int): Maximum number of entries to parse.
        existing_created (dict): Mapping from UniRef cluster ID to 'created' timestamp for idempotent imports.

    Returns:
        dict: Dictionary with lists for each CDM table
    """
    ns = {"ns": "http://uniprot.org/uniref"}  # Namespace for XML parsing
    entry_count = 0

    # Initialize lists to collect parsed rows for different tables
    cluster_data = []
    entity_data = []
    cluster_member_data = []
    cross_reference_data = []

    with gzip.open(local_gz, "rb") as f:
        # Stream parse the XML to avoid memory issues with big files
        context = ET.iterparse(f, events=("end",))
        for _, elem in context:
            if elem.tag.endswith("entry"):
                # Cluster basic info
                cluster_id, name = extract_cluster(elem, ns)

                # Get UniRef cluster id and timestamps
                uniref_id = elem.attrib.get("id")
                updated_time, created_time = get_timestamps(uniref_id, existing_created)

                # Populate Cluster and Entity table data
                cluster_data.append(
                    (
                        cluster_id,  # cluster_id
                        name,  # cluster name
                        "protein",  # entity_type (fixed value)
                        None,  # description (not present)
                        "UniRef 100",  # protocol_id
                    )
                )

                entity_data.append(
                    (
                        cluster_id,  # entity_id (matches cluster_id)
                        uniref_id,  # data_source_entity_id (UniRef100_xxx)
                        "Cluster",  # entity_type
                        "UniRef 100",  # data_source
                        updated_time,  # updated
                        created_time,  # created
                    )
                )

                # Extract UniProtKB member attributes and cross-references
                repr_db = elem.find("ns:representativeMember/ns:dbReference", ns)
                extract_cross_refs(repr_db, cross_reference_data, ns)

                for mem in elem.findall("ns:member/ns:dbReference", ns):
                    extract_cross_refs(mem, cross_reference_data, ns)

                # ClusterMember table (representative + members)
                add_cluster_members(cluster_id, repr_db, elem, cluster_member_data, ns)

                # Batch size limit
                entry_count += 1
                if entry_count >= batch_size:
                    break

                # Release element to save memory
                elem.clear()

    print(f"Parsed {entry_count} clusters")
    return {
        "cluster_data": cluster_data,
        "entity_data": entity_data,
        "cluster_member_data": cluster_member_data,
        "cross_reference_data": cross_reference_data,
    }


##### -------------- Save dalta table and print the preview --------------- #####


def save_delta_tables(spark, output_dir, data_dict) -> None:
    # Cluster
    cluster_schema = StructType(
        [
            StructField("cluster_id", StringType(), False),
            StructField("name", StringType(), False),
            StructField("entity_type", StringType(), False),
            StructField("description", StringType(), True),
            StructField("protocol_id", StringType(), False),
        ]
    )

    cluster_df = spark.createDataFrame(data_dict["cluster_data"], cluster_schema)
    cluster_df.write.format("delta").mode("overwrite").save(os.path.join(output_dir, "Cluster"))
    print(f"Cluster Delta table written to: {os.path.join(output_dir, 'Cluster')}")

    # Entity
    entity_schema = StructType(
        [
            StructField("entity_id", StringType(), False),
            StructField("data_source_entity_id", StringType(), False),
            StructField("entity_type", StringType(), False),
            StructField("data_source", StringType(), False),
            StructField("updated", StringType(), False),
            StructField("created", StringType(), False),
        ]
    )

    entity_df = spark.createDataFrame(data_dict["entity_data"], entity_schema)
    entity_table_path = os.path.join(output_dir, "Entity")
    entity_df.write.format("delta").mode("overwrite").save(entity_table_path)
    print(f"Entity Delta table written to: {entity_table_path}")

    # ClusterMember
    cluster_member_schema = StructType(
        [
            StructField("cluster_id", StringType(), False),
            StructField("entity_id", StringType(), False),
            StructField("is_representative", StringType(), False),
            StructField("is_seed", StringType(), False),
            StructField("score", StringType(), False),
        ]
    )

    cluster_member_df = spark.createDataFrame(data_dict["cluster_member_data"], cluster_member_schema)
    cluster_member_path = os.path.join(output_dir, "ClusterMember")
    cluster_member_df.write.format("delta").mode("overwrite").save(cluster_member_path)
    print(f"ClusterMember Delta table written to: {cluster_member_path}")

    # CrossReference
    cross_reference_schema = StructType(
        [
            StructField("entity_id", StringType(), False),
            StructField("xref_type", StringType(), False),
            StructField("xref_value", StringType(), False),
        ]
    )

    cross_reference_df = spark.createDataFrame(data_dict["cross_reference_data"], cross_reference_schema)
    cross_reference_path = os.path.join(output_dir, "CrossReference")
    cross_reference_df.write.format("delta").mode("overwrite").save(cross_reference_path)
    print(f"CrossReference Delta table written to: {cross_reference_path}")

    # Previews
    print("Sample Clusters:")
    cluster_df.createOrReplaceTempView("Cluster")
    spark.sql("SELECT * FROM Cluster LIMIT 20").show(truncate=False)

    print("Sample Entities:")
    entity_df.createOrReplaceTempView("Entity")
    spark.sql("SELECT * FROM Entity LIMIT 20").show(truncate=False)

    print("Sample ClusterMembers:")
    cluster_member_df.createOrReplaceTempView("ClusterMember")
    spark.sql("SELECT * FROM ClusterMember LIMIT 20").show(truncate=False)

    print("Sample CrossReferences:")
    cross_reference_df.createOrReplaceTempView("CrossReference")
    spark.sql("SELECT * FROM CrossReference LIMIT 20").show(truncate=False)


def build_spark_session():
    # Build a Spark session configured for Delta Lake support
    builder = (
        SparkSession.builder.appName("UniRef Cluster Extractor")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )
    return configure_spark_with_delta_pip(builder).getOrCreate()


# Click command-line interface for parameter parsing
@click.command()
@click.option("--ftp-url", required=True, help="FTP URL to UniRef100 XML file")
@click.option("--output-dir", required=True, help="Output directory for Delta table")
@click.option("--batch-size", default=1000, help="Number of UniRef entries to parse (limit)")
def main(ftp_url, output_dir, batch_size) -> None:
    # Set local path for downloaded gzipped XML file
    local_gz = os.path.join("/tmp", os.path.basename(ftp_url))

    # Download UniRef XML if not already present
    try:
        download_file(ftp_url, local_gz)
    except URLError as e:
        print(f"Error! Cannot download file: {e.reason}")
        return

    # Start Spark session with Delta Lake support
    spark = build_spark_session()

    # Load existing entity creation timestamps
    entity_table_path = os.path.join(output_dir, "Entity")
    existing_created = load_existing_created(spark, entity_table_path)

    # Parse the UniRef XML and extract all CDM table data
    data_dict = parse_uniref_xml(local_gz, batch_size, existing_created)

    # Write parsed data to Delta tables in output directory
    save_delta_tables(spark, output_dir, data_dict)

    spark.stop()


if __name__ == "__main__":
    main()
