"""UniProt XML Delta Lake Ingestion Pipeline.

This script parses UniProt XML (.xml.gz) file and ingests the data into structured Delta Lake tables.

Typical usage:
--------------
python3 src/parsers/uniprot.py \
    --xml-url "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz" \
    --output-dir "./output" \
    --namespace "uniprot_db" \
    --batch-size 5000

Use it in the local computer as:
python3 uniprot.py \
  --xml-url "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz" \
  --output-dir "./output_archaea" \
  --namespace "uniprot_archaea_db" \
  --batch-size 5000

Use it in the local computer as:
python3 uniprot.py \
  --xml-url "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz" \
  --output-dir "./output_archaea" \
  --namespace "uniprot_archaea_db" \
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

import datetime
import gzip
import json
import os
import uuid
import xml.etree.ElementTree as ET
from typing import Optional

import click
import requests
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.types import ArrayType, StringType, StructField, StructType

from cdm_data_loader_utils.parsers.shared_identifiers import parse_identifiers_generic
from cdm_data_loader_utils.parsers.xml_utils import clean_dict, find_all_text, get_attr, get_text, parse_db_references

## XML namespace mapping for UniProt entries (used for all XPath queries)
NS = {"ns": "https://uniprot.org/uniprot"}


def load_existing_identifiers(spark, output_dir, namespace):
    """
    Load the existing 'identifiers' Delta table and build a mapping from UniProt accession to CDM entity ID.
    This function enables consistent mapping of accessions to CDM IDs across multiple imports, supporting upsert and idempotent workflows.

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


def generate_cdm_id(prefix="cdm_prot_"):
    """
    Generate a CDM internal ID with a configurable prefix.
    Example: cdm_prot_<uuid>
    """
    return f"{prefix}{uuid.uuid4()}"


def build_datasource_record(xml_url):
    """
    Build a provenance record for the UniProt datasource without version extraction.
    """
    return {
        "name": "UniProt import",
        "source": "UniProt",
        "url": xml_url,
        "accessed": datetime.datetime.now(datetime.UTC).isoformat(),
        "version": 115,
    }


def parse_identifiers(entry, cdm_id):
    out = parse_identifiers_generic(entry=entry, xpath="ns:accession", prefix="UniProt", ns=NS)
    for row in out:
        row["entity_id"] = cdm_id
    return out


def _make_name_record(cdm_id: str, name_text: str, description: str) -> dict:
    """
    Small helper to create a standardized name record.
    """
    return {
        "entity_id": cdm_id,
        "name": name_text,
        "description": description,
        "source": "UniProt",
    }


def parse_names(entry, cdm_id):
    """
    Extract all protein names from a UniProt <entry> element, including
    - Top-level <name> elements (generic names)
    - <recommendedName> and <alternativeName> blocks within <protein> (full and short names)

    Returns:
        List[Dict[str, str]] with keys:
          - entity_id
          - name
          - description
          - source
    """

    names = []

    # Use `find_all_text` to automatically perform strip and deduplication
    top_level_names = find_all_text(entry, "ns:name", NS)
    for txt in top_level_names:
        names.append(
            _make_name_record(
                cdm_id=cdm_id,
                name_text=txt,
                description="UniProt protein name",
            )
        )

    # Recommended/alternative names in <protein> block
    protein = entry.find("ns:protein", NS)
    if protein is not None:
        # (tag_name, logical_type)
        name_blocks_spec = [
            ("recommendedName", "recommended"),
            ("alternativeName", "alternative"),
        ]

        # (xml tag, label used in description)
        length_spec = [
            ("fullName", "full"),
            ("shortName", "short"),
        ]

        for tag_name, logical_type in name_blocks_spec:
            for name_block in protein.findall(f"ns:{tag_name}", NS):
                for xml_tag, length_label in length_spec:
                    elem = name_block.find(f"ns:{xml_tag}", NS)
                    text = get_text(elem)
                    if not text:
                        continue
                    desc = f"UniProt {logical_type} {length_label} name"
                    names.append(
                        _make_name_record(
                            cdm_id=cdm_id,
                            name_text=text,
                            description=desc,
                        )
                    )

    return names


def parse_protein_info(entry, cdm_id):
    """
    Extract protein-level metadata from a UniProt XML <entry> element
    using shared XML utilities.
    """

    protein_info = {}

    # Extract EC numbers (recommended + alternative names)
    protein = entry.find("ns:protein", NS)
    if protein is not None:
        # full list of EC number search paths
        ec_paths = [
            "ns:recommendedName/ns:ecNumber",
            "ns:alternativeName/ns:ecNumber",
        ]

        ec_numbers = []
        for path in ec_paths:
            ec_numbers.extend(find_all_text(protein, path, NS))

        if ec_numbers:
            protein_info["ec_numbers"] = ec_numbers

    # Protein existence evidence
    protein_existence = entry.find("ns:proteinExistence", NS)
    if protein_existence is not None:
        protein_info["protein_id"] = cdm_id
        protein_info["evidence_for_existence"] = get_attr(protein_existence, "type")

    # Sequence block
    seq_elem = entry.find("ns:sequence", NS)
    if seq_elem is not None:
        protein_info.update(
            clean_dict(
                {
                    "length": get_attr(seq_elem, "length"),
                    "mass": get_attr(seq_elem, "mass"),
                    "checksum": get_attr(seq_elem, "checksum"),
                    "modified": get_attr(seq_elem, "modified"),
                    "sequence_version": get_attr(seq_elem, "version"),
                    "sequence": get_text(seq_elem),
                }
            )
        )

    # XML entry metadata (modified timestamp)
    entry_modified = get_attr(entry, "modified") or get_attr(entry, "updated")
    if entry_modified:
        protein_info["entry_modified"] = entry_modified

    # Return clean dict
    return protein_info if protein_info else None


def parse_evidence_map(entry):
    """
    Parse all <evidence> elements and return:
        { evidence_key : { evidence_type, supporting_objects, publications } }
    Using shared XML helpers for dbReference parsing and dict cleaning.
    """

    evidence_map = {}

    # Iterate over <evidence> elements
    for ev in entry.findall("ns:evidence", NS):
        key = get_attr(ev, "key")
        if not key:
            continue

        evidence_type = get_attr(ev, "type")

        # Parse <source><dbReference>
        pubs: list[str] = []
        others: list[str] = []

        source = ev.find("ns:source", NS)
        if source is not None:
            raw_pubs, raw_others = parse_db_references(source, NS)

            # Normalize publications:
            normalized_pubs: list[str] = []
            for p in raw_pubs:
                up = p.upper()
                if up.startswith("PUBMED:"):
                    _, acc = p.split(":", 1)
                    normalized_pubs.append(f"PMID:{acc}")
                else:
                    normalized_pubs.append(p)

            pubs = normalized_pubs
            others = raw_others

        # Build clean evidence metadata
        evidence_map[key] = clean_dict(
            {
                "evidence_type": evidence_type,
                "publications": pubs or None,
                "supporting_objects": others or None,
            }
        )

    return evidence_map


def parse_reaction_association(reaction, cdm_id, evidence_map):
    associations = []
    for dbref in reaction.findall("ns:dbReference", NS):
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
    for dbref in cofactor.findall("ns:dbReference", NS):
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


def _make_association(
    cdm_id: str,
    obj: str,
    predicate: Optional[str] = None,
    evidence_key: Optional[str] = None,
    evidence_map: Optional[dict] = None,
):
    """
    Helper to construct a CDM association dict and merge evidence fields from evidence_map.
    """
    assoc = {
        "subject": cdm_id,
        "object": obj,
        "predicate": predicate,
        "evidence_type": None,
        "supporting_objects": None,
        "publications": None,
    }

    if evidence_key and evidence_map and evidence_key in evidence_map:
        assoc.update(evidence_map[evidence_key])

    return clean_dict(assoc)


def parse_associations(entry, cdm_id, evidence_map):
    """
    Parse all relevant associations from a UniProt XML entry for the CDM model.
    Only include fields that are not None for each association.
    """
    associations = []

    # ---------------- Taxonomy association ----------------
    organism = entry.find("ns:organism", NS)
    if organism is not None:
        taxon_ref = organism.find('ns:dbReference[@type="NCBI Taxonomy"]', NS)
        if taxon_ref is not None:
            tax_id = taxon_ref.get("id")
            if tax_id:
                associations.append(
                    _make_association(
                        cdm_id=cdm_id,
                        obj=f"NCBITaxon:{tax_id}",
                        predicate=None,
                        evidence_key=None,
                        evidence_map=None,
                    )
                )

    # ------------- Database cross-references --------------
    for dbref in entry.findall("ns:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        if not db_type or not db_id:
            continue

        evidence_key = dbref.get("evidence")
        associations.append(
            _make_association(
                cdm_id=cdm_id,
                obj=f"{db_type}:{db_id}",
                predicate=None,
                evidence_key=evidence_key,
                evidence_map=evidence_map,
            )
        )

    # ------------- Catalytic activity / cofactor ----------
    for comment in entry.findall("ns:comment", NS):
        comment_type = comment.get("type")
        if comment_type == "catalytic activity":
            for reaction in comment.findall("ns:reaction", NS):
                for assoc in parse_reaction_association(reaction, cdm_id, evidence_map):
                    associations.append(clean_dict(assoc))
        elif comment_type == "cofactor":
            for cofactor in comment.findall("ns:cofactor", NS):
                for assoc in parse_cofactor_association(cofactor, cdm_id):
                    associations.append(clean_dict(assoc))

    return associations


def parse_publications(entry):
    """
    Extract all publication references from a UniProt XML <entry>.
    Uses shared parse_db_references() to gather PubMed/DOI IDs from <citation> blocks.

    Returns:
        List[str] of standardized publication IDs, e.g. ["PMID:12345", "DOI:10.1000/xyz"]
    """
    publications: list[str] = []

    for reference in entry.findall("ns:reference", NS):
        citation = reference.find("ns:citation", NS)
        if citation is None:
            continue

        # Reuse the shared dbReference parser
        raw_pubs, _ = parse_db_references(citation, NS)

        for p in raw_pubs:
            up = p.upper()
            # Expect formats like "PUBMED:" or "DOI:"
            if up.startswith("PUBMED:"):
                _, acc = p.split(":", 1)
                publications.append(f"PMID:{acc}")
            elif up.startswith("DOI:"):
                # Preserve DOI value, but normalize prefix to upper
                _, acc = p.split(":", 1)
                publications.append(f"DOI:{acc}")

    return list(dict.fromkeys(publications))


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


def download_file(url, output_path, chunk_size=8192, overwrite=False) -> None:
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
        for _event, element in context:
            # Check tag name, ignoring namespace
            if element.tag.endswith("entry"):
                yield element
                element.clear()


## ================================ SCHEMA =================================
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


def save_batches_to_delta(spark, tables, output_dir, namespace) -> None:
    """
    Persist batches of parsed records for each CDM table into Delta Lake format,
    and register them into Spark SQL as managed Delta tables under {namespace}.

    """

    # Ensure the namespace(database) exists
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")

    for table_name, (records, schema) in tables.items():
        if not records:
            continue  # Skip all empty tables

        # Directory path: <output_dir>/<namespace>/<table_name>
        delta_dir = os.path.abspath(os.path.join(output_dir, namespace, table_name))

        # Determine write mode based on directory presence
        mode = "append" if os.path.exists(delta_dir) else "overwrite"

        print(
            f"[DEBUG] Saving table '{namespace}.{table_name}' to {delta_dir} "
            f"with mode={mode}, record count={len(records)}"
        )

        try:
            # Convert to Spark DataFrame
            df = spark.createDataFrame(records, schema)

            # Write to Delta directory
            (df.write.format("delta").mode(mode).option("overwriteSchema", "true").save(delta_dir))

            # Register using Spark SQL
            spark.sql(f"""
                CREATE TABLE IF NOT EXISTS {namespace}.{table_name}
                USING DELTA
                LOCATION '{delta_dir}'
            """)

        except Exception as e:
            print(f"[ERROR] Failed to save '{table_name}' to Delta: {e}")


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
    Return (processed_entry_count, skipped_entry_count).

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
            main_accession_elem = entry_elem.find("ns:accession", NS)
            if main_accession_elem is None or main_accession_elem.text is None:
                skipped += 1
                continue

            # Generate a unique CDM ID (UUID) for this entry
            cdm_id = generate_cdm_id(prefix="cdm_prot_")

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


def ingest_uniprot(xml_url, output_dir, namespace, target_date=None, batch_size=5000) -> None:
    # Generate the timestamp for the current run
    current_timestamp = datetime.datetime.now(datetime.UTC).isoformat()

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
def main(xml_url, output_dir, namespace, target_date, batch_size) -> None:
    ingest_uniprot(
        xml_url=xml_url,
        output_dir=output_dir,
        namespace=namespace,
        target_date=target_date,
        batch_size=int(batch_size),
    )


if __name__ == "__main__":
    main()
