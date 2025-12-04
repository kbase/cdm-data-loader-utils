"""
UniProt XML Delta Lake Ingestion Pipeline
=========================================

This script parses UniProt XML (.xml.gz) file and ingests the data into structured Delta Lake tables.

Typical usage:
--------------
python3 src/parsers/uniprot.py \
    --xml-url "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz" \
    --output-dir "./output" \
    --namespace "uniprot_db" \
    --batch-size 5000

Arguments:
----------
--xml-url:      URL to the UniProt XML .gz file
--output-dir:   Output directory for Delta tables and logs (default: './output')
--namespace:    Delta Lake database name (default: 'uniprot_db')
--target-date:  Process entries modified/updated since specific date
--batch-size:   Number of UniProt entries to process per write batch (default: 5000)

Functionality:
--------------
- Downloads the XML file if not present locally
- Parses UniProt entries in a memory-efficient streaming fashion
- Maps parsed data into standardized CDM tables
- Writes all tables as Delta Lake tables, supporting incremental import
- Supports overwrite of previous imports and incremental updates by unique entity_id

Typical scenario:
-----------------
- Large-scale UniProt batch parsing and warehouse ingestion
- Efficient data lake ingestion

"""

import os
import click
import datetime
import json
import requests
import gzip
import uuid
import xml.etree.ElementTree as ET
from pyspark.sql import SparkSession
from delta import configure_spark_with_delta_pip
from pyspark.sql.types import ArrayType, StringType, StructType, StructField


## XML namespace mapping for UniProt entries (used for all XPath queries)
NS = {"u": "https://uniprot.org/uniprot"}


def load_existing_identifiers(spark, output_dir, namespace):
    """
    Load the existing 'identifiers' Delta table and build a mapping from UniProt accession to CDM entity ID.
    This function enables consistent mapping of accessions to CDM IDs across multiple imports, supporting upsert and idempotent workflows

    Returns:
        dict: {accession: entity_id}
    """

    access_to_cdm_id = {}
    id_path = os.path.abspath(os.path.join(output_dir, f"{namespace}_identifiers_delta"))
    if os.path.exists(id_path):
        try:
            # Read identifier and entity_id columns from the Delta table
            df = spark.read.format("delta").load(id_path).select("identifier", "entity_id")
            for row in df.collect():
                # Identifier field: UniProt:Pxxxxx, extract the actual accession part after the colon
                accession = row["identifier"].split(":", 1)[1]
                access_to_cdm_id[accession] = row["entity_id"]
        except Exception as e:
            print(f"Couldn't load identifiers table: {e}")
    else:
        print(f"No previous identifiers delta at {id_path}.")
    return access_to_cdm_id


def generate_cdm_id():
    """
    Generate a CDM entity_id directly from UniProt accession, using 'CDM:' prefix
    Ensures that each accession is mapped to stable and unique CDM entity ID, making it easy to join across different tables by accession.
    """
    return f"CDM:{uuid.uuid4()}"


def build_datasource_record(xml_url):
    """
    Build a provenance record for the UniProt datasource without version extraction.
    """
    return {
        "name": "UniProt import",
        "source": "UniProt",
        "url": xml_url,
        "accessed": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "version": 115,
    }


def parse_identifiers(entry, cdm_id):
    """
    Extract all accession numbers in the UniProt entry and format them into a CDM identifier structure
    """
    return [
        {
            "entity_id": cdm_id,
            "identifier": f"UniProt:{acc.text}",
            "source": "UniProt",
            "description": "UniProt accession",
        }
        for acc in entry.findall("u:accession", NS)
    ]


def parse_names(entry, cdm_id):
    """
    Extract all protein names from a UniProt <entry> element, including
    - Top-level <name> elements (generic names)
    - <recommendedName> and <alternativeName> blocks within <protein> (full and short names)
    """
    names = []

    # Extract all top-level <name> tags
    for name_element in entry.findall("u:name", NS):
        if name_element.text:
            names.append(
                {
                    "entity_id": cdm_id,
                    "name": name_element.text,
                    "description": "UniProt protein name",
                    "source": "UniProt",
                }
            )

    # Extract recommended and alternative names from <protein> block
    protein = entry.find("u:protein", NS)
    if protein is not None:
        for name_type in ["recommended", "alternative"]:
            # Directly use findall for simplicity (recommendedName returns single-element list)
            name_blocks = protein.findall(f"u:{name_type}Name", NS)
            for name in name_blocks:
                for name_length in ["full", "short"]:
                    name_string = name.find(f"u:{name_length}Name", NS)
                    if name_string is None or not name_string.text:
                        continue

                    names.append(
                        {
                            "entity_id": cdm_id,
                            "name": name_string.text,
                            "description": f"UniProt {name_type} {name_length} name",
                            "source": "UniProt",
                        }
                    )
    return names


def parse_protein_info(entry, cdm_id):
    """
    Extract protein-level metadata from a UniProt XML <entry> element
    """
    protein_info = {}
    ec_numbers = []

    # Extract EC numbers from <recommendedName> and <alternativeName> in <protein>
    protein = entry.find("u:protein", NS)
    if protein is not None:
        # Find EC numbers in recommendedName
        rec = protein.find("u:recommendedName", NS)
        if rec is not None:
            for ec in rec.findall("u:ecNumber", NS):
                if ec.text:
                    ec_numbers.append(ec.text)

        # Find EC numbers in all alternativeNames
        for alt in protein.findall("u:alternativeName", NS):
            for ec in alt.findall("u:ecNumber", NS):
                if ec.text:
                    ec_numbers.append(ec.text)
        if ec_numbers:
            protein_info["ec_numbers"] = ec_numbers

    # Extract protein existence evidence type
    protein_existence = entry.find("u:proteinExistence", NS)
    if protein_existence is not None:
        protein_info["protein_id"] = cdm_id
        protein_info["evidence_for_existence"] = protein_existence.get("type")

    # Extract sequence and sequence-related attributes
    seq_elem = entry.find("u:sequence", NS)
    if seq_elem is not None and seq_elem.text:
        protein_info["length"] = seq_elem.get("length")
        protein_info["mass"] = seq_elem.get("mass")
        protein_info["checksum"] = seq_elem.get("checksum")
        protein_info["modified"] = seq_elem.get("modified")
        protein_info["sequence_version"] = seq_elem.get("version")
        protein_info["sequence"] = seq_elem.text.strip()

    # Capture the entry's modified/updated date for tracking
    entry_modified = entry.attrib.get("modified") or entry.attrib.get("updated")
    if entry_modified:
        protein_info["entry_modified"] = entry_modified

    # Return the dictionary if any protein info was extracted
    return protein_info if protein_info else None


def parse_evidence_map(entry):
    """
    Parse all <evidence> elements from a UniProt XML entry and build a mapping
    from evidence key to metadata (type, supporting objects, publications)
    """
    evidence_map = {}

    # Loop through every <evidence> element in the entry
    for evidence in entry.findall("u:evidence", NS):
        key = evidence.get("key")  # Unique evidence key (string)
        evidence_type = evidence.get("type")  # Evidence code/type (e.g., ECO:0000255)

        supporting_objects = []
        publications = []

        # Check if this evidence has a <source> element with <dbReference> children
        source = evidence.find("u:source", NS)
        if source is not None:
            for dbref in source.findall("u:dbReference", NS):
                db_type = dbref.get("type")
                db_id = dbref.get("id")
                # Add publication references as PubMed or DOI; others as supporting objects
                if db_type == "PubMed":
                    publications.append(f"PMID:{db_id}")
                elif db_type == "DOI":
                    publications.append(f"DOI:{db_id}")
                else:
                    supporting_objects.append(f"{db_type}:{db_id}")

        # Store evidence metadata, omitting empty lists for cleanliness
        evidence_map[key] = {
            "evidence_type": evidence_type,
            "supporting_objects": supporting_objects if supporting_objects else None,
            "publications": publications if publications else None,
        }

    return evidence_map


def parse_reaction_association(reaction, cdm_id, evidence_map):
    associations = []
    for dbref in reaction.findall("u:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        assoc = {
            "subject": cdm_id,
            "predicate": "catalyzes",
            "object": f"{db_type}:{db_id}",
            "evidence_type": None,
            "supporting_objects": None,
            "publications": None,
        }
        evidence_key = reaction.get("evidence")
        if evidence_key and evidence_key in evidence_map:
            assoc.update(evidence_map[evidence_key])
        associations.append(assoc)
    return associations


def parse_cofactor_association(cofactor, cdm_id):
    associations = []
    for dbref in cofactor.findall("u:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        assoc = {
            "subject": cdm_id,
            "predicate": "requires_cofactor",
            "object": f"{db_type}:{db_id}",
            "evidence_type": None,
            "supporting_objects": None,
            "publications": None,
        }
        associations.append(assoc)
    return associations


def parse_associations(entry, cdm_id, evidence_map):
    """
    Parse all relevant associations from a UniProt XML entry for the CDM model.
    Only include fields that are not None for each association.
    """
    associations = []

    def clean(d):
        """Remove None-value keys from a dict."""
        return {k: v for k, v in d.items() if v is not None}

    # Taxonomy association
    organism = entry.find("u:organism", NS)
    if organism is not None:
        taxon_ref = organism.find('u:dbReference[@type="NCBI Taxonomy"]', NS)
        if taxon_ref is not None:
            associations.append(
                clean(
                    {
                        "subject": cdm_id,
                        "object": f"NCBITaxon:{taxon_ref.get('id')}",
                        "predicate": None,
                        "evidence_type": None,
                        "supporting_objects": None,
                        "publications": None,
                    }
                )
            )

    # Database cross-references with evidence
    for dbref in entry.findall("u:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        association = {
            "subject": cdm_id,
            "object": f"{db_type}:{db_id}",
            "predicate": None,
            "evidence_type": None,
            "supporting_objects": None,
            "publications": None,
        }
        evidence_key = dbref.get("evidence")
        if evidence_key and evidence_key in evidence_map:
            association.update(evidence_map[evidence_key])
        associations.append(clean(association))

    # Catalytic/cofactor
    for comment in entry.findall("u:comment", NS):
        comment_type = comment.get("type")
        if comment_type == "catalytic activity":
            # extract catalytic associations
            for reaction in comment.findall("u:reaction", NS):
                for assoc in parse_reaction_association(reaction, cdm_id, evidence_map):
                    associations.append(clean(assoc))
        elif comment_type == "cofactor":
            # extract cofactor associations
            for cofactor in comment.findall("u:cofactor", NS):
                for assoc in parse_cofactor_association(cofactor, cdm_id):
                    associations.append(clean(assoc))
    return associations


def parse_publications(entry):
    """
    Extract all publication references from a UniProt XML <entry>
    Returns a list of standardized publication IDs (PMID and DOI)
    """
    publications = []

    # Iterate through all <reference> blocks in the entry
    for reference in entry.findall("u:reference", NS):
        citation = reference.find("u:citation", NS)
        if citation is not None:
            # Each <citation> may have multiple <dbReference> elements (e.g., PubMed, DOI)
            for dbref in citation.findall("u:dbReference", NS):
                db_type = dbref.get("type")
                db_id = dbref.get("id")
                # Standardize format for known publication types
                if db_type == "PubMed":
                    publications.append(f"PMID:{db_id}")
                elif db_type == "DOI":
                    publications.append(f"DOI:{db_id}")

    return publications


def parse_uniprot_entry(entry, cdm_id, current_timestamp, datasource_name="UniProt import", prev_created=None):
    if prev_created:
        entity_created = prev_created
        entity_updated = current_timestamp
    else:
        entity_created = current_timestamp
        entity_updated = current_timestamp

    uniprot_created = entry.attrib.get("created")
    uniprot_modified = entry.attrib.get("modified") or entry.attrib.get("updated")
    uniprot_version = entry.attrib.get("version")
    entity = {
        "entity_id": cdm_id,
        "entity_type": "protein",
        "data_source": datasource_name,
        "created": entity_created,
        "updated": entity_updated,
        "version": uniprot_version,
        "uniprot_created": uniprot_created,
        "uniprot_modified": uniprot_modified,
    }
    evidence_map = parse_evidence_map(entry)
    return {
        "entity": entity,
        "identifiers": parse_identifiers(entry, cdm_id),
        "names": parse_names(entry, cdm_id),
        "protein": parse_protein_info(entry, cdm_id),
        "associations": parse_associations(entry, cdm_id, evidence_map),
        "publications": parse_publications(entry),
    }


def download_file(url, output_path, chunk_size=8192, overwrite=False):
    """
    Download a file from a given URL to a local output path.
    """
    # Skip download if file already exists and not overwriting
    if os.path.exists(output_path) and not overwrite:
        print(f"File '{output_path}' already exists.")
        return

    # Stream download to avoid high memory usage
    try:
        with requests.get(url, stream=True, timeout=60) as response:
            response.raise_for_status()
            with open(output_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
        print(f"Downloaded '{url}' to '{output_path}'")
    except Exception as e:
        print(f"Failed to download '{url}': {e}")

        if os.path.exists(output_path):
            os.remove(output_path)  # Remove incomplete file
        raise


def stream_uniprot_xml(filepath):
    """
    Stream and parse UniProt XML entries from a local gzipped file.
    Yields each <entry> element as soon as it is parsed to avoid loading the entire XML into memory.
    """
    # Open the gzipped XML file for reading in binary mode
    with gzip.open(filepath, "rb") as f:
        # Use iterparse to process XML incrementally, triggering on element end events
        context = ET.iterparse(f, events=("end",))
        for event, element in context:
            # Check tag name, ignoring namespace
            if element.tag.endswith("entry"):
                yield element
                element.clear()


## ================================ SCHEMA =================================
"""
Defines the Spark schema for all major CDM tables derived from UniProt XML.
Each schema is tailored for protein entities, identifiers, protein details, names, associations, and linked publications.
"""

schema_entities = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("entity_type", StringType(), False),
        StructField("data_source", StringType(), False),
        StructField("created", StringType(), True),
        StructField("updated", StringType(), True),
        StructField("version", StringType(), True),
        StructField("uniprot_created", StringType(), True),
        StructField("uniprot_modified", StringType(), True),
    ]
)

schema_identifiers = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("identifier", StringType(), False),
        StructField("source", StringType(), True),
        StructField("description", StringType(), True),
    ]
)

schema_proteins = StructType(
    [
        StructField("protein_id", StringType(), False),
        StructField("ec_numbers", StringType(), True),
        StructField("evidence_for_existence", StringType(), True),
        StructField("length", StringType(), True),
        StructField("mass", StringType(), True),
        StructField("checksum", StringType(), True),
        StructField("modified", StringType(), True),
        StructField("sequence_version", StringType(), True),
        StructField("sequence", StringType(), True),
        StructField("entry_modified", StringType(), True),
    ]
)

schema_names = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("name", StringType(), False),
        StructField("description", StringType(), True),
        StructField("source", StringType(), True),
    ]
)

schema_associations = StructType(
    [
        StructField("subject", StringType(), True),
        StructField("object", StringType(), True),
        StructField("predicate", StringType(), True),
        StructField("evidence_type", StringType(), True),
        StructField("supporting_objects", ArrayType(StringType()), True),
        StructField("publications", ArrayType(StringType()), True),
    ]
)

schema_publications = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("publication", StringType(), True),
    ]
)


def save_batches_to_delta(spark, tables, output_dir, namespace):
    """
    Persist batches of parsed records for each CDM table into Delta Lake format.

    - Each table is saved into a Delta directory named '{namespace}_{table}_delta' in the output folder.
    - If the Delta directory exists, append new records. Otherwise, overwrite it.
    - Registers the table in the Spark SQL for downstream query.
    """
    for table, (records, schema) in tables.items():
        if not records:
            continue  # Skip all empty tables

        delta_dir = os.path.abspath(os.path.join(output_dir, f"{namespace}_{table}_delta"))
        # Use "append" mode if the Delta directory already exists, otherwise "overwrite"
        mode = "append" if os.path.exists(delta_dir) else "overwrite"

        print(
            f"[DEBUG] Registering table: {namespace}.{table} at {delta_dir} with mode={mode}, record count: {len(records)}"
        )

        try:
            df = spark.createDataFrame(records, schema)
            df.write.format("delta").mode(mode).option("overwriteSchema", "true").save(delta_dir)
            spark.sql(f"""
                CREATE TABLE IF NOT EXISTS {namespace}.{table}
                USING DELTA
                LOCATION '{delta_dir}'
            """)
        except Exception as e:
            print(f"Failed to save {table} to Delta: {e}")


def prepare_local_xml(xml_url, output_dir):
    """
    Download the remote UniProt XML (.xml.gz) file to the specified local output directory,
    unless the file already exists locally. Returns the full local file path.
    """
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    local_xml_path = os.path.join(output_dir, os.path.basename(xml_url))
    # Download only if file does not exist
    download_file(xml_url, local_xml_path)
    return local_xml_path


def save_datasource_record(xml_url, output_dir):
    """
    Generate and save the datasource provenance record as a JSON file in the output directory.
    """
    datasource = build_datasource_record(xml_url)
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
    output_path = os.path.join(output_dir, "datasource.json")
    with open(output_path, "w") as f:
        json.dump(datasource, f, indent=4)
    return datasource


def get_spark_session(namespace):
    """
    Initialize SparkSession with Delta Lake support, and ensure the target database exists.
    """
    # Build SparkSession with Delta extensions enabled
    builder = (
        SparkSession.builder.appName("DeltaIngestion")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config(
            "spark.sql.catalog.spark_catalog",
            "org.apache.spark.sql.delta.catalog.DeltaCatalog",
        )
    )
    spark = configure_spark_with_delta_pip(builder).getOrCreate()
    # Ensure the target namespace (database) exists
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")
    return spark


def load_existing_entity(spark, output_dir, namespace):
    """
    Load the existing entities_delta Delta table and build a mapping of entity_id to created timestamp.
    This mapping is used to support upserts and idempotent writes.
    """
    old_created_dict = {}
    entities_table_path = os.path.abspath(os.path.join(output_dir, f"{namespace}_entities_delta"))
    if os.path.exists(entities_table_path):
        try:
            # Read only the required columns for efficiency
            old_df = spark.read.format("delta").load(entities_table_path).select("entity_id", "created")
            for row in old_df.collect():
                old_created_dict[row["entity_id"]] = row["created"]
            print(f"Loaded {len(old_created_dict)} existing entity_id records for upsert.")
        except Exception as e:
            print(f"Couldn't load previous entities delta table: {e}")
    else:
        print(f"No previous entities delta at {entities_table_path}.")
    return old_created_dict


def parse_entries(local_xml_path, target_date, batch_size, spark, tables, output_dir, namespace, current_timestamp):
    """
    Parse UniProt XML entries, write to Delta Lake in batches
    Return (processed_entry_count, skipped_entry_count)

    """

    target_date_dt = None

    # Convert target_date string to datetime for comparison if provided
    if target_date:
        try:
            target_date_dt = datetime.datetime.strptime(target_date, "%Y-%m-%d")
        except Exception:
            print(f"Invalid target date is {target_date}")

    entry_count, skipped = 0, 0

    # Iterate over each <entry> element in the XML file
    for entry_elem in stream_uniprot_xml(local_xml_path):
        try:
            # Get the modification date of the entry
            mod_date = entry_elem.attrib.get("modified") or entry_elem.attrib.get("updated")
            # If target_date is set, skip entries older than target_date
            if target_date_dt and mod_date:
                try:
                    entry_date_dt = datetime.datetime.strptime(mod_date[:10], "%Y-%m-%d")
                    if entry_date_dt < target_date_dt:
                        skipped += 1
                        continue
                except Exception:
                    skipped += 1
                    continue

            # Extract main accession (skip entry if not present)
            main_accession_elem = entry_elem.find("u:accession", NS)
            if main_accession_elem is None or main_accession_elem.text is None:
                skipped += 1
                continue

            # Generate a unique CDM ID (UUID) for this entry
            cdm_id = generate_cdm_id()

            # Parse all sub-objects: entity, identifiers, names, protein, associations, publications
            record = parse_uniprot_entry(entry_elem, cdm_id, current_timestamp)
            tables["entities"][0].append(record["entity"])
            tables["identifiers"][0].extend(record["identifiers"])
            tables["names"][0].extend(record["names"])

            if record["protein"]:
                tables["proteins"][0].append(record["protein"])
            tables["associations"][0].extend(record["associations"])
            tables["publications"][0].extend(
                {"entity_id": record["entity"]["entity_id"], "publication": pub} for pub in record["publications"]
            )

            entry_count += 1
            # Write batch to Delta and clear lists every batch_size entries
            if entry_count % batch_size == 0:
                save_batches_to_delta(spark, tables, output_dir, namespace)
                for v in tables.values():
                    v[0].clear()
                print(f"{entry_count} entries processed and saved")
        except Exception as e:
            # If any error occurs in parsing this entry, skip it and count
            print(f"Error parsing entry: {e}")
            skipped += 1
            continue

    # write remaining records
    save_batches_to_delta(spark, tables, output_dir, namespace)
    return entry_count, skipped


def ingest_uniprot(xml_url, output_dir, namespace, target_date=None, batch_size=5000):
    # Generate the timestamp for the current run
    current_timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()

    # Prepare local XML
    local_xml_path = prepare_local_xml(xml_url, output_dir)

    # Save data source meta information
    save_datasource_record(xml_url, output_dir)

    # Get Spark and the existing CDM entity_id
    spark = get_spark_session(namespace)

    # Define the table structure (batch storage)
    entities, identifiers, names, proteins, associations, publications = (
        [],
        [],
        [],
        [],
        [],
        [],
    )
    tables = {
        "entities": (entities, schema_entities),
        "identifiers": (identifiers, schema_identifiers),
        "names": (names, schema_names),
        "proteins": (proteins, schema_proteins),
        "associations": (associations, schema_associations),
        "publications": (publications, schema_publications),
    }

    # Main cycle processing, transfer to current timestamp
    entry_count, skipped = parse_entries(
        local_xml_path, target_date, batch_size, spark, tables, output_dir, namespace, current_timestamp
    )
    print(f"All entries processed ({entry_count}), skipped {skipped}, writing complete tables.")
    spark.sql(f"SHOW TABLES IN {namespace}").show()
    spark.sql(f"SELECT COUNT(*) FROM {namespace}.entities").show()

    # make sql test in entity table
    spark.sql(f"SELECT * FROM {namespace}.entities LIMIT 10").show(truncate=False)

    spark.stop()

    print(f"All Delta tables are created and registered in Spark SQL under `{namespace}`.")


@click.command()
@click.option("--xml-url", required=True, help="URL to UniProt XML (.xml.gz)")
@click.option("--output-dir", default="output", help="Output directory for Delta tables")
@click.option("--namespace", default="uniprot_db", help="Delta Lake database name")
@click.option(
    "--target-date",
    default=None,
    help="Only process entries modified/updated since this date (YYYY-MM-DD)",
)
@click.option("--batch-size", default=5000, help="Batch size for writing Delta tables")
def main(xml_url, output_dir, namespace, target_date, batch_size):
    ingest_uniprot(
        xml_url=xml_url,
        output_dir=output_dir,
        namespace=namespace,
        target_date=target_date,
        batch_size=int(batch_size),
    )


if __name__ == "__main__":
    main()
