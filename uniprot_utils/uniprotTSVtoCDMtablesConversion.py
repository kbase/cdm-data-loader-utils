import pandas as pd
import hashlib
import ast  # Safe evaluation of strings to Python literals


## This code utilize the Full_parsed_swissprot_data.tsv that was generated from uniprotDatFileParser.py
## and generate the CDM schema compliant parquet files and sample tsv files foreach table type

# Function to calculate SHA256 hash for sequences
def calculate_hash(sequence):
    return hashlib.sha256(sequence.encode('utf-8')).hexdigest()

# Load the parsed_swissprot_data.tsv file - a sample file that has the "Full_parsed_swissprot_data.tsv" is located in uniprotTest folder
input_file = "Full_parsed_swissprot_data.tsv"
data = pd.read_csv(input_file, sep='\t')

# Calculate SHA256 hash for sequences
data['hash'] = data['Sequence'].apply(calculate_hash)

# Remove duplicate sequences based on the hash
unique_data = data.drop_duplicates(subset=['hash'])

# 1. Generate the 'protein' table
protein_table = pd.DataFrame({
    'protein_hash': unique_data['hash'],  # Use hash as the primary key
    'name': unique_data['Entry'],
    'length': unique_data['Length'],
    'sequence': unique_data['Sequence'],
    'description': unique_data['Protein names']
})

# Save the 'protein' table as Parquet and TSV for QA/QC
protein_table.to_parquet("protein_table.parquet", index=False)
protein_table.head(20).to_csv("protein_table_sample.tsv", sep='\t', index=False)

# 2. Generate the 'name' table
def extract_source(gene_names):
    if pd.isna(gene_names) or not gene_names.strip():
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

name_table = pd.DataFrame({
    'protein_hash': unique_data['hash'],  # Use hash as the foreign key
    'name': unique_data['Entry'],
    'entry': unique_data['Entry Name'],
    'source': unique_data['Gene Names'].apply(extract_source),
    'description': unique_data['Protein names']
})

# Save the 'name' table as Parquet and TSV for QA/QC
name_table.to_parquet("name_table.parquet", index=False)
name_table.head(20).to_csv("name_table_sample.tsv", sep='\t', index=False)

# 3. Generate the 'identifier' table
identifier_table = pd.DataFrame({
    'protein_hash': unique_data['hash'],  # Use hash as the foreign key
    'identifier': unique_data['Entry Name'],
    'source': unique_data['Entry'].apply(lambda x: f"https://www.uniprot.org/uniprotkb/{x}/entry"),
    'description': unique_data['Protein names']
})

# Save the 'identifier' table as Parquet and TSV for QA/QC
identifier_table.to_parquet("identifier_table.parquet", index=False)
identifier_table.head(20).to_csv("identifier_table_sample.tsv", sep='\t', index=False)

# 4. Generate the 'association' table
def parse_ontologies(row):
    ontologies = []
    if not pd.isna(row['KEGG']):
        ontologies.append(f"KEGG: {row['KEGG']}")
    if not pd.isna(row['GO']):
        go_terms = [term.split(' ')[0] for term in row['GO'].split('; ') if term]
        ontologies.extend(go_terms)
    return ontologies

association_data = []
for _, row in unique_data.iterrows():
    protein_hash = row['hash']  # Use hash as the foreign key
    ontologies = parse_ontologies(row)
    for ontology in ontologies:
        association_data.append({
            'subject': protein_hash,
            'ontology_id': ontology,
            'predicate': row['Entry'],  # Protein name as the predicate
            'publications': row['Publications'] if not pd.isna(row['Publications']) else "NULL",
            'evidence_type': row['Evidence Codes'] if not pd.isna(row['Evidence Codes']) else "NULL"
        })

association_table = pd.DataFrame(association_data)

# Save the 'association' table as Parquet and TSV for QA/QC
association_table.to_parquet("association_table.parquet", index=False)
association_table.head(20).to_csv("association_table_sample.tsv", sep='\t', index=False)

# 5. Generate the 'feature_x_protein' table
feature_x_protein_data = []

for _, row in unique_data.iterrows():
    protein_hash = row['hash']  # Use hash as the foreign key
    gene_ids = row['GeneID'] if not pd.isna(row['GeneID']) else "NULL"

    if gene_ids != "NULL":
        for gene_id in gene_ids.split("; "):  # Handle multiple GeneIDs
            feature_x_protein_data.append({
                'protein_hash': protein_hash,
                'feature_id': gene_id.strip(),
                'protocol_id': "SwissProt/NCBI"
            })
    else:
        feature_x_protein_data.append({
            'protein_hash': protein_hash,
            'feature_id': "NULL",
            'protocol_id': "SwissProt/NCBI"
        })

feature_x_protein_table = pd.DataFrame(feature_x_protein_data)

# Save the 'feature_x_protein' table as Parquet and TSV for QA/QC
feature_x_protein_table.to_parquet("feature_x_protein.parquet", index=False)
feature_x_protein_table.head(20).to_csv("feature_x_protein_sample.tsv", sep='\t', index=False)

print("Tables generated and saved as Parquet and QA/QC TSV files.")
