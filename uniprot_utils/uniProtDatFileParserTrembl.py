import pandas as pd
import hashlib
import ast
import uuid
import os

# This code was specifically written to handle the trempl .dat file and able to parse into a tab delimited file
# Given the large size of trembl file, this code would able to process the trembl dataset as chunks and produce consolidated output file
# Janaka E

# Function to calculate SHA256 hash for sequences
def calculate_hash(sequence):
    return hashlib.sha256(sequence.encode('utf-8')).hexdigest() if pd.notnull(sequence) else None

# Function to generate UUID
def generate_uuid():
    return str(uuid.uuid4())

# Function to process each chunk
def process_chunk(chunk, hash_to_uuid_map, output_dir, chunk_number):
    # Treat missing values as NULL
    chunk = chunk.fillna("NULL")

    # Calculate SHA256 hash for sequences
    chunk['hash'] = chunk['Sequence'].apply(calculate_hash)

    # Remove duplicate sequences based on the hash
    unique_chunk = chunk.drop_duplicates(subset=['hash'])

    # Map UUIDs to each unique hash
    new_uuids = {hash_val: generate_uuid() for hash_val in unique_chunk['hash'] if hash_val not in hash_to_uuid_map}
    hash_to_uuid_map.update(new_uuids)

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Generate and save the 'protein' table for the chunk
    protein_data = pd.DataFrame({
        'protein_id': unique_chunk['hash'].map(hash_to_uuid_map),
        'name': unique_chunk['Entry'],
        'length': pd.to_numeric(unique_chunk['Length'], errors='coerce').fillna(0).astype(int),
        'sequence': unique_chunk['Sequence'],
        'hash': unique_chunk['hash'],
        'description': unique_chunk['Protein names']
    })
    protein_data.to_parquet(os.path.join(output_dir, f"protein_table_chunk_{chunk_number}.parquet"), index=False)

    # Generate and save the 'name' table for the chunk
    def extract_source(gene_names):
        if gene_names == "NULL" or not gene_names.strip():
            return "NULL"
        try:
            parsed_data = ast.literal_eval(gene_names)
            if isinstance(parsed_data, dict) and 'ORFNames' in parsed_data:
                orf_names = parsed_data.get('ORFNames', [])
                if isinstance(orf_names, list) and len(orf_names) > 0:
                    return orf_names[0]
            return "NULL"
        except (ValueError, SyntaxError):
            return "NULL"

    name_data = pd.DataFrame({
        'protein_id': unique_chunk['hash'].map(hash_to_uuid_map),
        'name': unique_chunk['Entry'],
        'entry': unique_chunk['Entry Name'],
        'source': unique_chunk['Gene Names'].apply(extract_source),
        'description': unique_chunk['Protein names']
    })
    name_data.to_parquet(os.path.join(output_dir, f"name_table_chunk_{chunk_number}.parquet"), index=False)

    # Generate and save the 'identifier' table for the chunk
    identifier_data = pd.DataFrame({
        'protein_id': unique_chunk['hash'].map(hash_to_uuid_map),
        'identifier': unique_chunk['Entry Name'],
        'source': unique_chunk['Entry'].apply(lambda x: f"https://www.uniprot.org/uniprotkb/{x}/entry"),
        'description': unique_chunk['Protein names']
    })
    identifier_data.to_parquet(os.path.join(output_dir, f"identifier_table_chunk_{chunk_number}.parquet"), index=False)

    # Generate and save the 'association' table for the chunk
    def parse_ontologies(row):
        ontologies = []
        if row['KEGG'] != "NULL":
            ontologies.append(f"KEGG: {row['KEGG']}")
        if row['GO'] != "NULL":
            go_terms = [term.split(' ')[0] for term in row['GO'].split('; ') if term]
            ontologies.extend(go_terms)
        return ontologies

    association_data = []
    for _, row in unique_chunk.iterrows():
        protein_id = hash_to_uuid_map[row['hash']]
        ontologies = parse_ontologies(row)
        for ontology in ontologies:
            association_data.append({
                'subject': protein_id,
                'ontology_id': ontology,
                'publications': row['Publications'] if row['Publications'] != "NULL" else "NULL",
                'evidence_type': row['Evidence Codes'] if row['Evidence Codes'] != "NULL" else "NULL"
            })
    pd.DataFrame(association_data).to_parquet(os.path.join(output_dir, f"association_table_chunk_{chunk_number}.parquet"), index=False)

    # Generate and save the 'feature_x_protein' table for the chunk
    feature_x_protein_data = []
    for _, row in unique_chunk.iterrows():
        protein_id = hash_to_uuid_map[row['hash']]
        gene_ids = row['GeneID'] if row['GeneID'] != "NULL" else None

        if gene_ids:
            for gene_id in gene_ids.split("; "):  # Handle multiple GeneIDs
                feature_x_protein_data.append({
                    'protein_id': protein_id,
                    'feature_id': gene_id.strip(),
                    'protocol_id': "SwissProt/NCBI"
                })
        else:
            feature_x_protein_data.append({
                'protein_id': protein_id,
                'feature_id': "NULL",
                'protocol_id': "SwissProt/NCBI"
            })
    pd.DataFrame(feature_x_protein_data).to_parquet(os.path.join(output_dir, f"feature_x_protein_chunk_{chunk_number}.parquet"), index=False)

    return len(new_uuids)

# Input and output file paths
input_file = "Full_parsed_trembl_data.tsv"
chunk_size = 1000000  # Process 500,000 lines at a time
output_dir = "output_parquet_files"

# Initialize hash-to-UUID mapping
hash_to_uuid_map = {}

# Process the input file in chunks
for chunk_number, chunk in enumerate(pd.read_csv(input_file, sep='\t', chunksize=chunk_size, dtype=str), start=1):
    print(f"Processing chunk {chunk_number}, lines read so far: {chunk_number * chunk_size}")
    process_chunk(chunk, hash_to_uuid_map, output_dir, chunk_number)

print(f"Tables generated and saved in {output_dir}.")
