"""
UniProt XML Delta Lake Ingestion Pipeline
=========================================

This script parses UniProt XML (.xml.gz) file and ingests the data into structured Delta Lake tables.

Typical usage:
--------------
python archaea_parsers.py \
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

Requirements:
-------------
- Python 3.7+
- pyspark, delta-spark, requests, click

Typical scenario:
-----------------
- Large-scale UniProt batch parsing and warehouse ingestion
- Efficient data lake ingestion

"""

import os
import click
import re
import hashlib
import datetime
import json
import requests
import gzip
import xml.etree.ElementTree as ET
from pyspark.sql import SparkSession
from delta import configure_spark_with_delta_pip
from pyspark.sql.types import ArrayType, StringType, StructType, StructField

## XML namespace mapping for UniProt entries (used for all XPath queries)
NS = {"u": "https://uniprot.org/uniprot"}


def generate_cdm_id(accession:str) -> str:
    """ 
    Generate a deterministic CDM entity_id based on a UniProt accession.
    Uses MD5 hash to ensure fixed length and avoid exposing raw accession.
    """
    if not accession or not isinstance(accession, str):
        raise ValueError("accession must be a non-empty string")
    normalized = accession.strip()
    md5_hash = hashlib.md5(normalized.encode("utf-8")).hexdigest()
    return f"CDM:{md5_hash}"


def extract_version(url):
    """
    Extracts UniProt version string from a download URL
    """
    patterns = [(r'release[-_]?(\d{4}_\d{2})', 1)]
    for pat, group in patterns:
        match_obj = re.search(pat, url, re.IGNORECASE)
        if match_obj:
            return match_obj.group(group)
    return None
   

def build_datasource_record(xml_url):
    """
    Build a provenance record for the UniProt datasource
    """
    name = "UniProt import"
    version = extract_version(xml_url) if xml_url else None
    accessed = int(datetime.datetime.now(datetime.timezone.utc).timestamp())
    return {
        "name": name,
        "source": "UniProt",
        "url": xml_url,
        "accessed": accessed,
        "version": version,
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
            "description": "UniProt accession"
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
            names.append({
                "entity_id": cdm_id,
                "name": name_element.text,
                "description": "UniProt protein name",
                "source": "UniProt"
            })

    # Extract recommended and alternative names from <protein> block
    protein = entry.find("u:protein", NS)
    if protein is not None:
        for name_type in ["recommended", "alternative"]:
            # recommendedName is singular, alternativeName can have multiple
            if name_type == "recommended":
                name_blocks = [protein.find("u:recommendedName", NS)]
            else:
                name_blocks = protein.findall("u:alternativeName", NS)
            for name in name_blocks:
                if name is None:
                    continue
                # both <fullName> and <shortName> under each name block
                for name_length in ["full", "short"]:
                    name_string = name.find(f"u:{name_length}Name", NS)
                    if name_string is None or not name_string.text:
                        continue
                    names.append({
                        "entity_id": cdm_id,
                        "name": name_string.text,
                        "description": f"UniProt {name_type} {name_length} name",
                        "source": "UniProt"
                    })
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
            # Join multiple EC numbers with commas for storage
            protein_info["ec_numbers"] = ",".join(ec_numbers)

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


def parse_associations(entry, cdm_id, evidence_map):
    """
    Parse all relevant associations from a UniProt XML entry for the CDM model.
    Includes taxonomy, database cross-references, catalytic activity, and cofactors.
    """
    associations = []

    # Taxonomy association: protein -> NCBI Taxonomy
    organism = entry.find("u:organism", NS)
    if organism is not None:
        taxon_ref = organism.find('u:dbReference[@type="NCBI Taxonomy"]', NS)
        if taxon_ref is not None:
            associations.append({
                "subject": cdm_id,
                "object": f"NCBITaxon:{taxon_ref.get('id')}",
                "predicate": None,
                "evidence_type": None,
                "supporting_objects": None,
                "publications": None
            })

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
            "publications": None
        }
        evidence_key = dbref.get("evidence")
        if evidence_key and evidence_key in evidence_map:
            association.update(evidence_map[evidence_key])
        associations.append(association)

    # Special comments include catalytic activity, cofactor
    for comment in entry.findall("u:comment", NS):
        comment_type = comment.get("type")
        if comment_type == "catalytic activity":
            # Each reaction may have dbReferences (Rhea or ChEBI IDs)
            for reaction in comment.findall("u:reaction", NS):
                for dbref in reaction.findall("u:dbReference", NS):
                    db_type = dbref.get("type")
                    db_id = dbref.get("id")
                    catalyze_assoc = {
                        "subject": cdm_id,
                        "predicate": "catalyzes",
                        "object": f"{db_type}:{db_id}",
                        "evidence_type": None,
                        "supporting_objects": None,
                        "publications": None
                    }
                    evidence_key = reaction.get("evidence")
                    if evidence_key and evidence_key in evidence_map:
                        catalyze_assoc.update(evidence_map[evidence_key])
                    associations.append(catalyze_assoc)

        elif comment_type == "cofactor":
            # Each cofactor may have dbReferences (ChEBI)
            for cofactor in comment.findall("u:cofactor", NS):
                for dbref in cofactor.findall("u:dbReference", NS):
                    db_type = dbref.get("type")
                    db_id = dbref.get("id")
                    cofactor_assoc = {
                        "subject": cdm_id,
                        "predicate": "requires_cofactor",
                        "object": f"{db_type}:{db_id}",
                        "evidence_type": None,
                        "supporting_objects": None,
                        "publications": None
                    }
                    associations.append(cofactor_assoc)
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


def parse_uniprot_entry(entry, datasource_name="UniProt import", prev_created=None):
    """
    Parse a single UniProt <entry> XML element into CDM-compatible records
    """

    # Extract primary accession as CDM unique id 
    main_accession = entry.find("u:accession", NS)
    if main_accession is None or main_accession.text is None:
        raise ValueError("Cannot generate entity_id")
    cdm_id = generate_cdm_id(main_accession.text)

    # Handle the CDM created/updated timestamps
    desire_date = datetime.datetime.now(datetime.timezone.utc).isoformat()
    if prev_created:
        entity_created = prev_created          # Use previous created date if updating existing record
        entity_updated = desire_date           # Update the 'updated' timestamp
    else:
        entity_created = desire_date          
        entity_updated = desire_date

    # Extract UniProt native metadata for provenance 
    uniprot_created = entry.attrib.get("created")
    uniprot_modified = entry.attrib.get("modified") or entry.attrib.get("updated")
    uniprot_version = entry.attrib.get("version")

    # Build CDM entity record, includes metadata extension 
    entity = {
        "entity_id": cdm_id,
        "entity_type": "protein",
        "data_source": datasource_name,
        "created": entity_created,                # CDM record created timestamp
        "updated": entity_updated,                # CDM record last updated timestamp
        "version": uniprot_version,               # UniProt version string
        "uniprot_created": uniprot_created,       # UniProt record creation date
        "uniprot_modified": uniprot_modified      # UniProt last modified date
    }

    # Parse and collect related records for each CDM table 
    evidence_map = parse_evidence_map(entry)  
    return {
        "entity": entity,
        "identifiers": parse_identifiers(entry, cdm_id),                  # All UniProt accessions
        "names": parse_names(entry, cdm_id),                              # Protein names
        "protein": parse_protein_info(entry, cdm_id),                     # Sequence, EC, existence
        "associations": parse_associations(entry, cdm_id, evidence_map),  # Cross-refs, taxonomy, evidence
        "publications": parse_publications(entry)                         # Literature links
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
            with open(output_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:  
                        f.write(chunk)
        print(f"[download_file] Downloaded '{url}' to '{output_path}'")
    except Exception as e:
        print(f"[download_file] Failed to download '{url}': {e}")

        if os.path.exists(output_path):
            os.remove(output_path)  # Remove incomplete file
        raise


def stream_uniprot_xml(filepath):
    """
    Stream and parse UniProt XML entries from a local gzipped file.
    Yields each <entry> element as soon as it is parsed to avoid loading the entire XML into memory.
    """
    # Open the gzipped XML file for reading in binary mode
    with gzip.open(filepath, 'rb') as f:
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

schema_entities = StructType([
    StructField("entity_id", StringType(), False),
    StructField("entity_type", StringType(), False),
    StructField("data_source", StringType(), False),
    StructField("created", StringType(), True),
    StructField("updated", StringType(), True),
    StructField("version", StringType(), True),
    StructField("uniprot_created", StringType(), True),     
    StructField("uniprot_modified", StringType(), True),   
])
schema_identifiers = StructType([
    StructField("entity_id", StringType(), False),
    StructField("identifier", StringType(), False),
    StructField("source", StringType(), True),
    StructField("description", StringType(), True)
])
schema_proteins = StructType([
    StructField("protein_id", StringType(), False),
    StructField("ec_numbers", StringType(), True),
    StructField("evidence_for_existence", StringType(), True),
    StructField("length", StringType(), True),
    StructField("mass", StringType(), True),
    StructField("checksum", StringType(), True),
    StructField("modified", StringType(), True),
    StructField("sequence_version", StringType(), True),
    StructField("sequence", StringType(), True),
    StructField("entry_modified", StringType(), True)
])
schema_names = StructType([
    StructField("entity_id", StringType(), False),
    StructField("name", StringType(), False),
    StructField("description", StringType(), True),
    StructField("source", StringType(), True)
])
schema_associations = StructType([
    StructField("subject", StringType(), True),
    StructField("object", StringType(), True),
    StructField("predicate", StringType(), True),
    StructField("evidence_type", StringType(), True),
    StructField("supporting_objects", ArrayType(StringType()), True),
    StructField("publications", ArrayType(StringType()), True)
])
schema_publications = StructType([
    StructField("entity_id", StringType(), False),
    StructField("publication", StringType(), True)
])


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
        try:
            df = spark.createDataFrame(records, schema)
            df.write.format("delta").mode(mode).option("overwriteSchema", "true").save(delta_dir)
            spark.sql(f"""
                CREATE TABLE IF NOT EXISTS {namespace}.{table}
                USING DELTA
                LOCATION '{delta_dir}'
            """)
        except Exception as e:
            print(f"[WARN] Failed to save {table} to Delta: {e}")


def prepare_local_xml(xml_url, output_dir):
    """
    Download the remote UniProt XML (.xml.gz) file to the specified local output directory,
    unless the file already exists locally. Returns the full local file path.
    """
    os.makedirs(output_dir, exist_ok=True)  # Ensure output directory exists
    local_xml_path = os.path.join(output_dir, os.path.basename(xml_url))
    download_file(xml_url, local_xml_path)   # Download only if file does not exist
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
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
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


def parse_entries(local_xml_path, old_created_dict, target_date, batch_size, spark, tables, output_dir, namespace):
    """
    Parses XML entry, writes to Delta in batches, returns entry_count, skipped_count
    """
    target_date_dt = None
    if target_date:
        try:
            target_date_dt = datetime.datetime.strptime(target_date, "%Y-%m-%d")
        except Exception:
            print(f"Invalid target date: {target_date}")

    entry_count, skipped = 0, 0

    for entry_elem in stream_uniprot_xml(local_xml_path):
        try:
            mod_date = entry_elem.attrib.get("modified") or entry_elem.attrib.get("updated")
            if target_date_dt and mod_date:
                try:
                    entry_date_dt = datetime.datetime.strptime(mod_date[:10], "%Y-%m-%d")
                    if entry_date_dt < target_date_dt:
                        skipped += 1
                        continue
                except Exception:
                    skipped += 1
                    continue

            main_accession_elem = entry_elem.find("u:accession", NS)
            if main_accession_elem is None or main_accession_elem.text is None:
                skipped += 1
                continue
            main_accession = main_accession_elem.text
            cdm_id = generate_cdm_id(main_accession)
            prev_created = old_created_dict.get(cdm_id)
            record = parse_uniprot_entry(entry_elem, prev_created=prev_created)
            tables["entities"][0].append(record["entity"])
            tables["identifiers"][0].extend(record["identifiers"])
            tables["names"][0].extend(record["names"])
            if record["protein"]:
                tables["proteins"][0].append(record["protein"])
            tables["associations"][0].extend(record["associations"])
            tables["publications"][0].extend(
                {"entity_id": record["entity"]["entity_id"], "publication": pub}
                for pub in record["publications"]
            )
            entry_count += 1
            if entry_count % batch_size == 0:
                save_batches_to_delta(spark, tables, output_dir, namespace)
                for v in tables.values():
                    v[0].clear()
                print(f"{entry_count} entries processed and saved...")

        except Exception as e:
            print(f"Error parsing entry: {e}")
            skipped += 1
            continue

    save_batches_to_delta(spark, tables, output_dir, namespace)
    return entry_count, skipped


def ingest_uniprot(xml_url,output_dir,namespace,target_date=None,batch_size=5000):
    # Prepare local XML
    local_xml_path = prepare_local_xml(xml_url, output_dir)

    # Save data source meta information
    save_datasource_record(xml_url, output_dir)

    # Get Spark and the existing CDM entity_id
    spark = get_spark_session(namespace)
    old_created_dict = load_existing_entity(spark, output_dir, namespace)

    # Define the table structure (batch storage)
    entities, identifiers, names, proteins, associations, publications = [], [], [], [], [], []
    tables = {
        "entities": (entities, schema_entities),
        "identifiers": (identifiers, schema_identifiers),
        "names": (names, schema_names),
        "proteins": (proteins, schema_proteins),
        "associations": (associations, schema_associations),
        "publications": (publications, schema_publications)
    }

    # Main cycle processing
    entry_count, skipped = parse_entries(
        local_xml_path, old_created_dict, target_date, batch_size, spark, tables, output_dir, namespace
    )
    print(f"All entries processed ({entry_count}), skipped {skipped}, writing complete tables.")
    spark.sql(f"SHOW TABLES IN {namespace}").show()
    spark.sql(f"SELECT COUNT(*) FROM {namespace}.entities").show()
    spark.stop()
    print(f"All Delta tables are created and registered in Spark SQL under `{namespace}`.")


@click.command()
@click.option('--xml-url', required=True, help='URL to UniProt XML (.xml.gz)')
@click.option('--output-dir', default='output', help='Output directory for Delta tables')
@click.option('--namespace', default='uniprot_db', help='Delta Lake database name')
@click.option('--target-date', default=None, help='Only process entries modified/updated since this date (YYYY-MM-DD)')
@click.option('--batch-size', default=5000, help='Batch size for writing Delta tables')

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
