import pandas as pd
from Bio import SwissProt
from io import StringIO

# This code parses the Uniprot .dat dumps into a tab delimited file
# Janaka E

def parse_swissprot_file_in_chunks(file_path, output_file, chunk_size=5000000):
    def write_to_file(records, output_file, write_mode):
        """Write processed records to the output file."""
        df = pd.DataFrame(records)
        df.to_csv(output_file, sep="\t", index=False, mode=write_mode, header=(write_mode == 'w'))

    def process_buffer(buffer, output_file, chunk_number):
        """Process lines in the buffer and write to the output file."""
        records = []
        try:
            buffer_as_string = "".join(buffer)
            buffer_stream = StringIO(buffer_as_string)
            for record in SwissProt.parse(buffer_stream):
                try:
                    # Extract evidence codes from features
                    evidence_codes = []
                    for feature in record.features:
                        if isinstance(feature.qualifiers, dict) and "evidence" in feature.qualifiers:
                            evidence_values = feature.qualifiers.get("evidence", [])
                            evidence_codes.extend(evidence_values if isinstance(evidence_values, list) else [evidence_values])
                    evidence_codes = "; ".join(map(str, evidence_codes)) if evidence_codes else "NULL"

                    # Extract publication details
                    publications = []
                    for ref in record.references:
                        authors = ref.authors.strip() if ref.authors else "Unknown Authors"
                        title = ref.title.strip() if ref.title else "No Title"
                        location = ref.location.strip() if ref.location else "No Journal Info"
                        pubmed = ""
                        for db, id_ in ref.references:
                            if db.lower() == "pubmed":
                                pubmed = f"PubMed:{id_}"
                        publication = f"{authors}. {title}. {location}. {pubmed}".strip()
                        publications.append(publication)
                    publications = "; ".join(publications) if publications else "NULL"

                    # Extract GO terms
                    go_terms = [
                        f"{xref[1]} ({xref[2]})"
                        for xref in record.cross_references
                        if len(xref) >= 3 and xref[0] == "GO"
                    ]
                    go_terms = "; ".join(map(str, go_terms)) if go_terms else "NULL"

                    # Consolidate all key data
                    entry_data = {
                        "Entry": record.accessions[0] if record.accessions else "NULL",
                        "Entry Name": record.entry_name if record.entry_name else "NULL",
                        "Reviewed": "Reviewed" if record.data_class == "Reviewed" else "Unreviewed",
                        "Protein names": record.description if record.description else "NULL",
                        "Gene Names": "; ".join(map(str, record.gene_name)) if record.gene_name else "NULL",
                        "Organism": record.organism if record.organism else "NULL",
                        "Taxonomy": "; ".join(map(str, record.organism_classification)) if record.organism_classification else "NULL",
                        "Length": len(record.sequence) if record.sequence else "NULL",
                        "Sequence": record.sequence if record.sequence else "NULL",
                        "PE": record.protein_existence if hasattr(record, "protein_existence") and record.protein_existence else "NULL",
                        "EMBL": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "EMBL") or "NULL",
                        "RefSeq": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "RefSeq") or "NULL",
                        "GeneID": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "GeneID") or "NULL",
                        "PDB": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "PDB") or "NULL",
                        "KEGG": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "KEGG") or "NULL",
                        "Reactome": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "Reactome") or "NULL",
                        "HGNC": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "HGNC") or "NULL",
                        "STRING": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "STRING") or "NULL",
                        "BioCyc": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "BioCyc") or "NULL",
                        "Pfam": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "Pfam") or "NULL",
                        "InterPro": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "InterPro") or "NULL",
                        "GO": go_terms,
                        "Proteomes": "; ".join(id_ for db, id_, *_ in record.cross_references if db == "Proteomes") or "NULL",
                        "Keywords": "; ".join(map(str, record.keywords)) if record.keywords else "NULL",
                        "Evidence Codes": evidence_codes,
                        "Publications": publications,
                    }
                    records.append(entry_data)
                except Exception as e:
                    print(f"Error parsing record: {e}")
                    continue
        except Exception as e:
            print(f"Error processing buffer: {e}")

        # Write to file
        if records:
            write_mode = 'a' if chunk_number > 1 else 'w'
            write_to_file(records, output_file, write_mode)

    # Main logic to read in chunks
    chunk_number = 1
    line_buffer = []
    with open(file_path, "r") as handle:
        for line in handle:
            line_buffer.append(line)
            if line.startswith("//") and len(line_buffer) >= chunk_size:
                print(f"Processing chunk {chunk_number}, lines read so far: {len(line_buffer)}")
                process_buffer(line_buffer, output_file, chunk_number)
                line_buffer = []  # Clear buffer
                chunk_number += 1

        # Process any remaining lines
        if line_buffer:
            print(f"Processing final chunk, total lines read: {len(line_buffer)}")
            process_buffer(line_buffer, output_file, chunk_number)

    print(f"Parsed data saved to {output_file}")

# Input and output file paths
# File paths are written to match the file paths in Sequoia -
#input_file = "/home/janakae/scratch/Uniprot/Trembl/uniprot_trembl.dat"  # Path to your UniProt text file
#output_file = "/home/janakae/scratch/Uniprot/Trembl/Full_parsed_trembl_data.tsv"

#Test files in the repo
input_file = "/uniprotTest/uniprotTest.dat"  # Path to your UniProt text file
output_file = "/uniprotTest/Full_parsed_swissprot_data_test.tsv"

# Parse the UniProt/SwissProt file in chunks
parse_swissprot_file_in_chunks(input_file, output_file)
