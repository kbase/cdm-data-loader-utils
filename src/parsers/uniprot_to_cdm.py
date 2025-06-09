import uuid
import xml.etree.ElementTree as ET
from datetime import date
import json
import logging
from pyspark.sql import SparkSession
import requests
import gzip


NS = {'u': 'https://uniprot.org/uniprot'}

def generate_unique_id():
    """
    Generate a unique identifier for the CDM.
    """
    return f'CDM:{uuid.uuid4()}'

def build_datasource_record(): 
    """
    Construct provenance information for the data source.
    Add in DataSource Table with provenance information for Uniprot download.
    """
    return {    
        'name': 'UniProt archaea',
        'source': 'UniProt',
        'url': 'https://uniprot.org/downloads/path/to/file',
        'accessed': date.today().strftime('%Y-%m-%d'),
        'version': 115
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
            'entity_id': cdm_id,
            'identifier': f"UniProt:{acc.text}",
            'source': 'UniProt',
            'description': 'UniProt accession'
        }
        for acc in entry.findall('u:accession', NS)
    ]

def parse_names(entry, cdm_id):
    """
    Extract all relevant protein names from a UniProt XML entry.
    Include <name> elements and all recommended and alternative names from the <protein> element.
    
    Args:
        entry (xml.etree.ElementTree.Element): The UniProt XML <entry> element.
        cdm_id (str): The unique CDM entity ID assigned to this protein.

    Returns:
        list of dict: Each dictionary represents a CDM 'name' record with entity_id, name, description, and source.
    """
    names = []

    # Extract top-level (simple protein name) elements 
    for name_element in entry.findall('u:name', NS):
        if name_element.text:
            names.append({
                'entity_id': cdm_id,
                'name': name_element.text,
                'description': 'UniProt protein name',
                'source': 'UniProt'
            })

    # Parse the <protein> block for recommended and alternative names
    protein = entry.find('u:protein', NS)
    if protein is not None:
        # List the recommendedName
        recommend = protein.find('u:recommendedName', NS)
        if recommend is not None:
            # Recommended fullName 
            full = recommend.find('u:fullName', NS)
            if full is not None and full.text:
                names.append({
                    'entity_id': cdm_id, 
                    'name': full.text, 
                    'description': 'UniProt recommended full name', 
                    'source': 'UniProt'})
                
            # Recommended shortName
            short = recommend.find('u:shortName', NS)
            if short is not None and short.text:
                names.append({
                    'entity_id': cdm_id, 
                    'name': short.text, 
                    'description': 'UniProt recommended short name', 
                    'source': 'UniProt'})
                
        # List the alternativeName
        for alt in protein.findall('u:alternativeName', NS):
            # Alternative fullName
            alternative_full = alt.find('u:fullName', NS)
            if alternative_full is not None and alternative_full.text:
                names.append({
                    'entity_id': cdm_id, 
                    'name': alternative_full.text, 
                    'description': 'UniProt alternative full name', 
                    'source': 'UniProt'})
                
            # Alternative shortName
            alternative_short = alt.find('u:shortName', NS)
            if alternative_short is not None and alternative_short.text:
                names.append({
                    'entity_id': cdm_id, 
                    'name': alternative_short.text, 
                    'description': 'UniProt alternative short name', 
                    'source': 'UniProt'})
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
    protein = entry.find('u:protein', NS)
    if protein is not None:
        # From <recommendedName> block 
        rec = protein.find('u:recommendedName', NS)
        if rec is not None:
            for ec in rec.findall('u:ecNumber', NS):
                if ec.text:
                    ec_numbers.append(ec.text)
        # From <alternativeName> block
        for alt in protein.findall('u:alternativeName', NS):
            for ec in alt.findall('u:ecNumber', NS):
                if ec.text:
                    ec_numbers.append(ec.text)
        if ec_numbers:
            protein_info['ec_numbers'] = ec_numbers

    # Extract protein existence evidence 
    protein_existence = entry.find('u:proteinExistence', NS)
    if protein_existence is not None:
        protein_info['protein_id'] = cdm_id
        protein_info['evidence_for_existence'] = protein_existence.get('type')

    # Extract sequence information 
    seq_elem = entry.find('u:sequence', NS)
    if seq_elem is not None and seq_elem.text:
        # Add all relevant sequence attributes according to examples
        protein_info['length'] = seq_elem.get('length')
        protein_info['mass'] = seq_elem.get('mass')
        protein_info['checksum'] = seq_elem.get('checksum')
        protein_info['modified'] = seq_elem.get('modified')
        protein_info['sequence_version'] = seq_elem.get('version')
        protein_info['sequence'] = seq_elem.text.strip()

    ## Add modified/update in entry-level 
    entry_modified = entry.attrib.get("modified") or entry.attrib.get("updated")
    if entry_modified:
        protein_info['entry_modified'] = entry_modified
    
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
    for evidence in entry.findall('u:evidence', NS):
        key = evidence.get('key') # Unique key for this evidence
        types = evidence.get('type') # Evidence type (ECO)

        supporting_objects = []
        publications = []

        # If <source> sub-element exists, extract all <dbReference> children
        source = evidence.find('u:source', NS)
        if source is not None:
            for dbref in source.findall('u:dbReference', NS):
                db_type = dbref.get('type')
                db_id = dbref.get('id')

                # Differentiate publication (PubMed) from other supporting objects
                if db_type == 'PubMed':
                    publications.append(f"PMID:{db_id}")
                else:
                    supporting_objects.append(f"{db_type}:{db_id}")

        # Store evidence information, omit empty lists
        evidence_map[key] = {
            "evidence_type": types,
            "supporting_objects": supporting_objects if supporting_objects else None,
            "publications": publications if publications else None
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
    organism = entry.find('u:organism', NS)
    if organism is not None:
        taxon_ref = organism.find('u:dbReference[@type="NCBI Taxonomy"]', NS)
        if taxon_ref is not None:
            associations.append({
                "subject": cdm_id, 
                "object": f"NCBITaxon:{taxon_ref.get('id')}"
                })
            
    # dbReference associations 
    for dbref in entry.findall('u:dbReference', NS):
        db_type = dbref.get('type')
        db_id = dbref.get('id')
        assoc = {
            "subject": cdm_id,
            "object": f"{db_type}:{db_id}"
            }

        # If there is evidence, merge in evidence attributes
        evidence_key = dbref.get('evidence')
        if evidence_key is not None and evidence_key in evidence_map:
            assoc.update(evidence_map[evidence_key])
        associations.append(assoc)

    # Associations from comment types, include catalytic activity and cofactor 
    for comment in entry.findall('u:comment', NS):
        comment_type = comment.get('type')

        if comment_type == "catalytic activity":
            for reaction in comment.findall('u:reaction', NS):
                # For each dbReference under the reaction, create a catalyzes association
                for dbref in reaction.findall('u:dbReference', NS):
                    db_type = dbref.get('type')
                    db_id = dbref.get('id')
                    catalyze_assoc = {
                        "subject": cdm_id,
                        "predicate": "catalyzes",
                        "object": f"{db_type}:{db_id}"
                    }

                    # Evidence can be linked from the reaction element
                    evidence_key = reaction.get('evidence')
                    if evidence_key is not None:
                        if evidence_key in evidence_map:
                            catalyze_assoc.update(evidence_map[evidence_key])
                    associations.append(catalyze_assoc)

        elif comment_type == "cofactor":
            for cofactor in comment.findall('u:cofactor', NS):
                # Each cofactor can have multiple dbReferences
                for dbref in cofactor.findall('u:dbReference', NS):
                    db_type = dbref.get('type')
                    db_id = dbref.get('id')
                    cofactor_assoc = {
                        "subject": cdm_id,
                        "predicate": "requires_cofactor",
                        "object": f"{db_type}:{db_id}"
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
    for refer in entry.findall('u:reference', NS):
        citation = refer.find('u:citation', NS)
        if citation is not None:
            # Loop all <dbReference> elements within <citation>
            for dbref in citation.findall('u:dbReference', NS):
                db_type = dbref.get('type')
                db_id = dbref.get('id')

                # Format according to type
                if db_type == 'PubMed':
                    publications.append(f'PMID:{db_id}')
                elif db_type == 'DOI':
                    publications.append(f'DOI:{db_id}')
                elif db_type in ['EMBL', 'GenBank', 'DDBJ']:
                    publications.append(f'{db_type}:{db_id}')

    return publications


def parse_uniprot_entry(entry):
    cdm_id = generate_unique_id()
    entity = {
        'entity_id': cdm_id,
        'entity_type': 'protein',
        'data_source': 'UniProt archaea',
        'created': entry.attrib.get('created'),
        'updated': entry.attrib.get('modified') or entry.attrib.get('updated'),
        'version': entry.attrib.get('version')
    }

    evidence_map = parse_evidence_map(entry)

    return {
        'entity': entity,
        'identifiers': parse_identifiers(entry, cdm_id),
        'names': parse_names(entry, cdm_id),
        'protein': parse_protein_info(entry, cdm_id),
        'associations': parse_associations(entry, cdm_id, evidence_map),
        'publications': parse_publications(entry)
    }

def stream_uniprot_xml_gz(url): 
    ## Download and decompress data from URL 
    with requests.get(url, stream = True) as r: 
        r.raise_for_status()

        with gzip.GzipFile(fileobj=r.raw) as f:
            context = ET.iterparse(f, events=('end',))
            for event, element in context:
                # Only process <entry> element
                if element.tag.endswith('entry'):
                    yield element
                    element.clear() 


def read_jsonl(file_path):
    with open(file_path, 'r') as f:
        return [json.loads(line) for line in f]


if __name__ == "__main__":
    xml_url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz"
    datasource = build_datasource_record()
    target_date = None ## In SwissProt, "2025-04-09" for filtering 

    ### Step 1. Print all jsonl files to be created
    with open("datasource.json", "w") as f:
        json.dump(datasource, f, indent=4)

    with open("entities.jsonl", "w") as entity_out, \
         open("identifiers.jsonl", "w") as identifier_out, \
         open("names.jsonl", "w") as names_out, \
         open("proteins.jsonl", "w") as proteins_out, \
         open("associations.jsonl", "w") as associations_out, \
         open("publications.jsonl", "w") as publications_out:

        try:
            for entry_elem in stream_uniprot_xml_gz(xml_url):
                mod_date = entry_elem.attrib.get("modified") or entry_elem.attrib.get("updated")

                if target_date and mod_date != target_date:
                    continue

                try:
                    record = parse_uniprot_entry(entry_elem)
                    json.dump(record['entity'], entity_out); entity_out.write('\n')
                    for iden in record['identifiers']:
                        json.dump(iden, identifier_out); identifier_out.write('\n')
                    for nm in record['names']:
                        json.dump(nm, names_out); names_out.write('\n')
                    if record['protein']:
                        json.dump(record['protein'], proteins_out); proteins_out.write('\n')
                    for assoc in record['associations']:
                        json.dump(assoc, associations_out); associations_out.write('\n')
                    for pub in record['publications']:
                        json.dump({'entity_id': record['entity']['entity_id'], 'publication': pub}, publications_out); publications_out.write('\n')
                except Exception as e:
                    logging.exception("Failed to process one entry (skipped).")
        except Exception as e:
            logging.exception("Fatal error in streaming or file writing.")

    print("Done! All CDM tables exported as .jsonl files in the folder")


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


    # Step 3. Convert to Parquet format
    parquet_table = ['entities', 'identifiers', 'names', 'proteins', 'associations', 'publications']
    spark = SparkSession.builder.appName("CDM Ingestion").getOrCreate()

    for key in parquet_table:
        # Read all jsonl files
        df = spark.read.json(f"{key}.jsonl")
        # Write as parquet
        df.write.mode("overwrite").parquet(f"{key}.parquet")
        print(f"{key}.parquet written via Spark.")

    spark.stop()
    print("Table written as parquet")

