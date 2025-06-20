import uuid
import xml.etree.ElementTree as ET
from datetime import date
import json
import logging
import requests
import gzip
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, udf
from pyspark.sql.types import StringType
import os


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
        # From <recommendedName> block
        rec = protein.find("u:recommendedName", NS)
        if rec is not None:
            for ec in rec.findall("u:ecNumber", NS):
                if ec.text:
                    ec_numbers.append(ec.text)
        # From <alternativeName> block
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
    This mapping can be used to enrich associations and other data objects,
    with evidence types, supporting objects, and publication references.

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

    # Loop each <reference> element in the entry
    for refer in entry.findall("u:reference", NS):
        citation = refer.find("u:citation", NS)
        if citation is not None:
            # Loop all <dbReference> elements within <citation>
            for dbref in citation.findall("u:dbReference", NS):
                db_type = dbref.get("type")
                db_id = dbref.get("id")

                # Format according to type
                if db_type == "PubMed":
                    publications.append(f"PMID:{db_id}")
                elif db_type == "DOI":
                    publications.append(f"DOI:{db_id}")
                elif db_type in ["EMBL", "GenBank", "DDBJ"]:
                    publications.append(f"{db_type}:{db_id}")

    return publications


def parse_uniprot_entry(entry):
    cdm_id = generate_unique_id()
    entity = {
        "entity_id": cdm_id,
        "entity_type": "protein",
        "data_source": "UniProt archaea",
        "created": entry.attrib.get("created"),
        "updated": entry.attrib.get("modified") or entry.attrib.get("updated"),
        "version": entry.attrib.get("version"),
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


def stream_uniprot_xml_gz(url):
    ## Download and decompress data from URL
    with requests.get(url, stream=True) as r:
        r.raise_for_status()

        with gzip.GzipFile(fileobj=r.raw) as f:
            context = ET.iterparse(f, events=("end",))
            for event, element in context:
                # Only process <entry> element
                if element.tag.endswith("entry"):
                    yield element
                    element.clear()


def read_jsonl(file_path):
    with open(file_path, "r") as f:
        return [json.loads(line) for line in f]


if __name__ == "__main__":
    xml_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz"
    datasource = build_datasource_record()
    target_date = None  

    ### Step 1. Print all jsonl files to be created
    with open("datasource.json", "w") as f:
        json.dump(datasource, f, indent=4)

    ## Parsing XML to jsonl
    with (
        open("entities.jsonl", "w") as entity_out,
        open("identifiers.jsonl", "w") as identifier_out,
        open("names.jsonl", "w") as names_out,
        open("proteins.jsonl", "w") as proteins_out,
        open("associations.jsonl", "w") as associations_out,
        open("publications.jsonl", "w") as publications_out,
    ):
        try:
            for entry_elem in stream_uniprot_xml_gz(xml_url):
                mod_date = entry_elem.attrib.get("modified") or entry_elem.attrib.get("updated")

                if target_date and mod_date != target_date:
                    continue

                try:
                    record = parse_uniprot_entry(entry_elem)
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
                        json.dump(
                            {
                                "entity_id": record["entity"]["entity_id"],
                                "publication": pub
                            },
                            publications_out,
                        )
                        publications_out.write("\n")
                except Exception:
                    logging.exception("Failed to process one entry skipped.")
        except Exception:
            logging.exception("Fatal error in streaming or file writing.")

    print("All CDM tables exported as .jsonl files in the folder")

    # Step 2. Read all jsonl files and write to a single cdm_data.json
    cdm_tables = {
        "datasource": json.load(open("datasource.json")),
        "entities": read_jsonl("entities.jsonl"),
        "identifiers": read_jsonl("identifiers.jsonl"),
        "names": read_jsonl("names.jsonl"),
        "proteins": read_jsonl("proteins.jsonl"),
        "associations": read_jsonl("associations.jsonl"),
        "publications": read_jsonl("publications.jsonl"),
    }
    with open("cdm_data.json", "w") as f:
        json.dump(cdm_tables, f, indent=2)
    print("cdm_data.json completed with all tables.")

    def array_to_string(arr):
        if arr is None:
            return None
        return ";".join(arr)
    
    array_to_string_udf = udf(array_to_string, StringType())

    # Step 3. jsonl -> CSV for delta 
    csv_tables = ["entities","identifiers","names","proteins","associations","publications"]

    ### builder = SparkSession.builder.appName("DeltaIngestion")

    builder = (SparkSession.builder.appName("DeltaIngestion")
               .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
               .config("spark.sql.catalog.spark_catalog", "org.apache.spark.sql.delta.catalog.DeltaCatalog")
    )

    spark = configure_spark_with_delta_pip(builder).getOrCreate()
    csv_output_dir = "./"

    for key in csv_tables:
        # Read all jsonl files
        df = spark.read.json(f"{key}.jsonl")

        for field in df.schema.fields:
            if str(field.dataType).startswith("ArrayType"):
                print(f"Converting array column '{field.name}' to string in {key} table")
                df = df.withColumn(field.name, array_to_string_udf(col(field.name)))

        # Write to a csv directory with a folder for each table to import directly into delta
        csv_dir = f"{csv_output_dir}{key}.csv"
        df.write.mode("overwrite").option("header", True).csv(csv_dir)
        print(f"{key}.csv written via Spark.")

    # Step 4. CSV -> Delta table 
    namespace = "archaebase_db"
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")
    spark.sql(f"USE {namespace}")

    for key in csv_tables:
        #csv_dir = f"{csv_output_dir}{key}.csv"
        #delta_dir = f"{csv_output_dir}{key}_delta"
        csv_dir = os.path.abspath(f"{csv_output_dir}/{key}.csv")
        delta_dir = os.path.abspath(f"{csv_output_dir}/{key}_delta")
        df = spark.read.option("header", True).csv(csv_dir)
        df.write.format("delta").mode("overwrite").save(delta_dir)

        spark.sql(f"""
            CREATE TABLE IF NOT EXISTS {namespace}.{key}
            USING DELTA
            LOCATION '{delta_dir}'
        """)

        spark.sql(f"REFRESH TABLE {namespace}.{key}")
        print(f"{key} Delta Table created at {delta_dir}")

    spark.sql(f"SHOW TABLES IN {namespace}").show(truncate=False)
    df = spark.sql(f"SELECT COUNT(*) FROM {namespace}.entities")
    df.show()


    print("All Delta tables created")
    spark.stop()

