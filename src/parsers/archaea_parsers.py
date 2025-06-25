import uuid
import xml.etree.ElementTree as ET
from datetime import date
import json
import requests
import gzip
from pyspark.sql import SparkSession
from delta import configure_spark_with_delta_pip
from pyspark.sql.types import StringType
import os
import shutil
from pyspark.sql.types import (ArrayType, StructType, StructField)


NS = {"u": "https://uniprot.org/uniprot"}

def generate_unique_id():
    """
    Generate a unique identifier for the CDM.
    """
    return f"CDM:{uuid.uuid4()}"


def build_datasource_record():
    """
    Construct provenance information for the data source.
    Add in DataSource Table with provenance information for Uniprot download.
    """
    return {
        "name": "UniProt archaea",
        "source": "UniProt",
        "url": "https://uniprot.org/downloads/path/to/file",
        "accessed": date.today().strftime("%Y-%m-%d"),
        "version": 115,
    }


def parse_identifiers(entry, cdm_id):
    """
    Extract UniProt accession numbers from the XML entry and format them as CDM identifier records.

    Args:
        entry (xml.etree.ElementTree.Element): The UniProt XML <entry> element.
        cdm_id (str): The unique CDM entity ID assigned to this protein.

    Returns:
        list of dict: Each dictionary represents a CDM identifier record, including the entity_id, identifier (with 'UniProt:' prefix), source, and description.
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
    This function extracts all protein names from the UniProt XML <entry>:

    - Adds all top-level <name> as generic names
    - Iterates over <recommendedName> and all <alternativeName> blocks under <protein>
    - For each, collects both <fullName> and <shortName> if present
    - Results are labeled with name type (recommended/alternative) and name length (full/short)
    - Output is a flat list of dictionaries for CDM ingestion
    """

    names = []

    # Extract all top-level <name> tags as generic protein names
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
    
    # Parses the <protein> block with recommended and alternative names
    protein = entry.find("u:protein", NS)
    if protein is not None:
        for name_type in ["recommended", "alternative"]:
            # Only one recommended name
            # have several alternative names
            if name_type == "recommended":
                name_blocks = [protein.find("u:recommendedName", NS)]
            else: 
                name_blocks = protein.findall("u:alternativeName", NS)
            
            # Iterate over all name_blocks (recommended names/alternative names)
            for name in name_blocks:
                if name is None:
                    continue
                
                # extract fullName and shortName respectively
                for name_length in ["full", "short"]:
                    name_string = name.find(f"u:{name_length}Name", NS)
                    # If the field is not present or the content is empty, then skip
                    if name_string is None or not name_string.text:
                        continue

                    # Add to names list, description field reflects type and length
                    names.append({
                        "entity_id": cdm_id,
                        "name": name_string.text,
                        "description": f"UniProt {name_type} {name_length} name",
                        "source": "UniProt"
                    })
    return names


def parse_protein_info(entry, cdm_id):
    """
    Extract protein-level metadata from a UniProt XML entry, including EC numbers,
    existence evidence, and sequence information.

    Args:
        entry (xml.etree.ElementTree.Element): The UniProt XML <entry> element.
        cdm_id (str): The unique CDM entity ID assigned to this protein.

    Returns:
        dict or None: Protein info record for the CDM 'protein' table, or None if empty.
    """

    protein_info = {}
    ec_numbers = []

    # Extract EC numbers from <protein>
    protein = entry.find("u:protein", NS)
    if protein is not None:
        # From <recommendedName> block: ecNumbers are usually here
        rec = protein.find("u:recommendedName", NS)
        if rec is not None:
            for ec in rec.findall("u:ecNumber", NS):
                if ec.text:
                    ec_numbers.append(ec.text)

        # From <alternativeName> block: ecNumbers can also be here
        for alt in protein.findall("u:alternativeName", NS):
            for ec in alt.findall("u:ecNumber", NS):
                if ec.text:
                    ec_numbers.append(ec.text)
        if ec_numbers:
            protein_info["ec_numbers"] = ",".join(ec_numbers)

    # Extract protein existence evidence
    protein_existence = entry.find("u:proteinExistence", NS)
    if protein_existence is not None:
        protein_info["protein_id"] = cdm_id
        protein_info["evidence_for_existence"] = protein_existence.get("type")

    # Extract sequence information
    seq_elem = entry.find("u:sequence", NS)
    if seq_elem is not None and seq_elem.text:
        # Add all relevant sequence attributes according to examples
        # Note: 'length' and 'mass' are optional, use get() to avoid KeyError
        protein_info["length"] = seq_elem.get("length")
        protein_info["mass"] = seq_elem.get("mass")
        protein_info["checksum"] = seq_elem.get("checksum")
        protein_info["modified"] = seq_elem.get("modified")
        protein_info["sequence_version"] = seq_elem.get("version")
        protein_info["sequence"] = seq_elem.text.strip()

    ## Add modified/update in entry-level
    entry_modified = entry.attrib.get("modified") or entry.attrib.get("updated")
    if entry_modified:
        protein_info["entry_modified"] = entry_modified

    # Return dictionary if info was collected
    return protein_info if protein_info else None


def parse_evidence_map(entry):
    """
    Parse the <evidence> elements from a UniProt XML entry and create a mapping from evidence keys to evidence metadata.
    This mapping can be used to enrich associations and other data objects, with evidence types, supporting objects, and publication references.

    Args:
        entry (xml.etree.ElementTree.Element): The <entry> element from UniProt XML.

    Returns:
        dict: A mapping from evidence keys to dictionaries with evidence_type, supporting_objects, and publications.
    """

    evidence_map = {}

    # Iterate over all <evidence> elements in the entry
    for evidence in entry.findall("u:evidence", NS):
        key = evidence.get("key")  # Unique key for this evidence
        types = evidence.get("type")  # Evidence type (ECO)

        supporting_objects = []
        publications = []

        # If <source> sub-element exists, extract all <dbReference> children
        source = evidence.find("u:source", NS)
        if source is not None:
            for dbref in source.findall("u:dbReference", NS):
                db_type = dbref.get("type")
                db_id = dbref.get("id")

                # Differentiate publication (PubMed) from other supporting objects
                if db_type == "PubMed":
                    publications.append(f"PMID:{db_id}")
                else:
                    supporting_objects.append(f"{db_type}:{db_id}")

        # Store evidence information, omit empty lists
        evidence_map[key] = {
            "evidence_type": types,
            "supporting_objects": supporting_objects if supporting_objects else None,
            "publications": publications if publications else None,
        }

    return evidence_map


def parse_associations(entry, cdm_id, evidence_map):
    """
    Parse associations for a UniProt entry, including:
    - Organism (taxonomy ID) associations
    - Cross-database references with evidence
    - Special comment-derived associations (catalytic activity, cofactors)

    """

    associations = []

    # Organism association (link protein to NCBI taxonomy)
    organism = entry.find("u:organism", NS)
    if organism is not None:
        taxon_ref = organism.find('u:dbReference[@type="NCBI Taxonomy"]', NS)
        if taxon_ref is not None:
            associations.append(
                {"subject": cdm_id, "object": f"NCBITaxon:{taxon_ref.get('id')}"}
            )

    # dbReference associations
    for dbref in entry.findall("u:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        assoc = {"subject": cdm_id, "object": f"{db_type}:{db_id}"}

        # If there is evidence, merge in evidence attributes
        evidence_key = dbref.get("evidence")
        if evidence_key is not None and evidence_key in evidence_map:
            assoc.update(evidence_map[evidence_key])
        associations.append(assoc)

    # Associations from comment types, include catalytic activity and cofactor
    for comment in entry.findall("u:comment", NS):
        comment_type = comment.get("type")

        if comment_type == "catalytic activity":
            for reaction in comment.findall("u:reaction", NS):
                # For each dbReference under the reaction, create a catalyzes association
                for dbref in reaction.findall("u:dbReference", NS):
                    db_type = dbref.get("type")
                    db_id = dbref.get("id")
                    catalyze_assoc = {
                        "subject": cdm_id,
                        "predicate": "catalyzes",
                        "object": f"{db_type}:{db_id}",
                    }

                    # Evidence can be linked from the reaction element
                    evidence_key = reaction.get("evidence")
                    if evidence_key is not None:
                        if evidence_key in evidence_map:
                            catalyze_assoc.update(evidence_map[evidence_key])
                    associations.append(catalyze_assoc)

        elif comment_type == "cofactor":
            for cofactor in comment.findall("u:cofactor", NS):
                # Each cofactor can have multiple dbReferences
                for dbref in cofactor.findall("u:dbReference", NS):
                    db_type = dbref.get("type")
                    db_id = dbref.get("id")
                    cofactor_assoc = {
                        "subject": cdm_id,
                        "predicate": "requires_cofactor",
                        "object": f"{db_type}:{db_id}",
                    }
                    # Add more evidence logic if needed
                    associations.append(cofactor_assoc)

    return associations


def parse_publications(entry):
    """
    Extract publication from a UniProt XML entry.

    This function parses <reference> elements and their <citation> sub-elements,
    retrieving all related dbReference IDs.

    Returns a list of formatted publication identifiers.

    """

    publications = []

    # Loop through each <reference> element in the UniProt XML <entry>
    for refer in entry.findall("u:reference", NS):
        citation = refer.find("u:citation", NS)
        # Get the <citation> element inside the reference
        if citation is not None:
            # Iterate over all <dbReference> elements within the <citation> block
            for dbref in citation.findall("u:dbReference", NS):
                db_type = dbref.get("type") ## PubMed, DOI, EMBL
                db_id = dbref.get("id")

                # Format and append based on database type
                if db_type == "PubMed":
                    publications.append(f"PMID:{db_id}")
                elif db_type == "DOI":
                    publications.append(f"DOI:{db_id}")
                elif db_type in ["EMBL", "GenBank", "DDBJ"]:
                    publications.append(f"{db_type}:{db_id}")

    return publications


def parse_uniprot_entry(entry):
    """
    Parse a single <entry> element from a UniProt XML file into structured CDM-compatible records.

    This function orchestrates parsing across different UniProt XML substructures, extracting
    information such as entity metadata, identifiers, protein names, sequence-level, associations, 
    and linked publications.

    Args:
        entry (xml.etree.ElementTree.Element): The UniProt <entry> XML element.

    Returns:
        dict: A dictionary containing extracted records across CDM tables:
              - 'entity': core metadata
              - 'identifiers': UniProt accessions
              - 'names': protein naming data
              - 'protein': sequence and EC info
              - 'associations': taxonomy, cross-references, evidence-backed links
              - 'publications': linked literature
    """
    ## Generate a globally unique CDM entity ID (UUID-based)
    cdm_id = generate_unique_id()

    # Build the core entity record for the protein, capturing creation/modification metadata
    entity = {
        "entity_id": cdm_id, 
        "entity_type": "protein", # CDM classification
        "data_source": "UniProt archaea", # Provenance
        "created": entry.attrib.get("created"), # Date entry was first added
        "updated": entry.attrib.get("modified") or entry.attrib.get("updated"), # Last modified
        "version": entry.attrib.get("version"), # UniProt version number
    }
    
    # Parse all <evidence> tags into a mapping from evidence keys to metadata
    evidence_map = parse_evidence_map(entry)
    
    # Return all extracted components organized under CDM table keys
    return {
        "entity": entity,
        "identifiers": parse_identifiers(entry, cdm_id),
        "names": parse_names(entry, cdm_id),
        "protein": parse_protein_info(entry, cdm_id),
        "associations": parse_associations(entry, cdm_id, evidence_map),
        "publications": parse_publications(entry),
    }


def stream_uniprot_xml_gz(url):
    """
    Stream and parse a UniProt XML .gz file from a remote URL.

    This function uses streaming HTTP and on-the-fly GZIP decompression to efficiently
    iterate over large UniProt XML files without loading the entire content into memory.

    Args:
        url (str): URL to the gzipped UniProt XML file.

    Yields:
        xml.etree.ElementTree.Element: Each complete <entry> element in the XML file.

    """

    # Make a streaming HTTP GET request to avoid loading the whole file
    with requests.get(url, stream=True) as r:
        r.raise_for_status()

        with gzip.GzipFile(fileobj=r.raw) as f:
            # Use iterparse to stream through the XML
            context = ET.iterparse(f, events=("end",))
            for event, element in context:
                # Yield only when a full <entry> element is parsed
                if element.tag.endswith("entry"): 
                    yield element
                    element.clear()
    

## Define the Spark schema for each CDM table
## These schemas ensure correct data types and structure when reading JSONL files into Spark DataFrames

# Schema for the 'entities' table
# basic metadata about each protein entity
schema_entities = StructType([
    StructField("entity_id", StringType(), False),
    StructField("entity_type", StringType(), False),
    StructField("data_source", StringType(), False),
    StructField("created", StringType(), True),
    StructField("updated", StringType(), True),
    StructField("version", StringType(), True)
    ])

# Schema for the 'identifiers' table: UniProt accession and related external IDs
schema_identifiers = StructType([
    StructField("entity_id", StringType(), False),
    StructField("identifier", StringType(), False),
    StructField("source", StringType(), True),
    StructField("description", StringType(), True)
    ])

# Schema for the 'proteins' table: detailed protein information including sequence and EC numbers
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

# Schema for the 'names' table: recommended and alternative names
schema_names = StructType([
    StructField("entity_id", StringType(), False),
    StructField("name", StringType(), False),
    StructField("description", StringType(), True),
    StructField("source", StringType(), True)
    ])

# Schema for the 'associations' table: relationships to taxonomy, references, cofactors
schema_associations = StructType([
    StructField("subject", StringType(), True),
    StructField("object", StringType(), True),
    StructField("predicate", StringType(), True),
    StructField("evidence_type", StringType(), True),
    StructField("supporting_objects", ArrayType(StringType()), True),
    StructField("publications", ArrayType(StringType()), True)
    ])

# Schema for the 'publications' table: standalone list of publication IDs
schema_publications = StructType([
    StructField("entity_id", StringType(), False),
    StructField("publication", StringType(), True)
    ])


if __name__ == "__main__":
    xml_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz"
    
    datasource = build_datasource_record()
    target_date = None  ## change to a specific date if needed

    ### Step 1. Print and check all jsonl files for each CDM table
    with open("datasource.json", "w") as f:
        json.dump(datasource, f, indent=4)
    
    # Open output JSONL files for each CDM entity table
    with (
        open("entities.jsonl", "w") as entity_out,
        open("identifiers.jsonl", "w") as identifier_out,
        open("names.jsonl", "w") as names_out,
        open("proteins.jsonl", "w") as proteins_out,
        open("associations.jsonl", "w") as associations_out,
        open("publications.jsonl", "w") as publications_out
        ):

        for entry_elem in stream_uniprot_xml_gz(xml_url):
            ## Skip the entries that don't match a specific modified date
            mod_date = entry_elem.attrib.get("modified") or entry_elem.attrib.get("updated")
            if target_date and mod_date != target_date:
                continue
            
            ## Parse the UniProt entry into CDM-format record dictionaries
            record = parse_uniprot_entry(entry_elem)

            ## List each parsed object to its corresponding JSONL file
            json.dump(record["entity"], entity_out)
            entity_out.write("\n")

            for iden in record["identifiers"]:
                json.dump(iden, identifier_out)
                identifier_out.write("\n")
            for nm in record["names"]:
                json.dump(nm, names_out)
                names_out.write("\n")
            if record["protein"]:
                json.dump(record["protein"], proteins_out)
                proteins_out.write("\n")
            for assoc in record["associations"]:
                json.dump(assoc, associations_out)
                associations_out.write("\n")
            for pub in record["publications"]:
                json.dump({"entity_id": record["entity"]["entity_id"], "publication": pub}, publications_out)
                publications_out.write("\n")
    print("Exported JSONL files.")
    
    ## Initialize Spark session with Delta Lake 
    builder = (
        SparkSession.builder.appName("DeltaIngestion")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
        )
    spark = configure_spark_with_delta_pip(builder).getOrCreate()
    
    ## Dictionary mapping table names to the Spark schemas
    schemas = {
        "entities": schema_entities,
        "identifiers": schema_identifiers,
        "names": schema_names,
        "proteins": schema_proteins,
        "associations": schema_associations,
        "publications": schema_publications
        }
    
    ## Create a Delta Lake database archaebase_db 
    namespace = "archaebase_db"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")
    
    ## For each table, read its JSONL data with Delta format
    for table, schema in schemas.items():
        df = spark.read.schema(schema).json(f"{table}.jsonl")
        delta_dir = os.path.abspath(f"./{table}_delta")
        
        ## Clean up any existing Delta directory
        if os.path.exists(delta_dir):
            shutil.rmtree(delta_dir)
        
        ## Write the DataFrame as a Delta table
        df.write.format("delta").mode("overwrite").save(delta_dir)
        
        ## Register the Delta table in Spark SQL 
        spark.sql(f"""
            CREATE TABLE IF NOT EXISTS {namespace}.{table}
            USING DELTA
            LOCATION '{delta_dir}'
        """)

    ##  Show all tables created in the namespace and count of entities
    spark.sql(f"SHOW TABLES IN {namespace}").show()
    df_count = spark.sql(f"SELECT COUNT(*) FROM {namespace}.entities")
    df_count.show()

    spark.stop()
    print("All Delta tables are created and registered in Spark SQL.")
