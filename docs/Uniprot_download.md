# Uniprot API & Downloads Overview <br> 
Date: 2025-05-08 <br> 

## Background 
To understand how the UniProt API can be used, it’s important to first understand what UniProtKB is. <br>

UniProtKB (Universal Protein Resource Knowledgebase) is a comprehensive, high-quality database of protein functional information. It consists of two main sections:
- Swiss-Prot: Manually reviewed entries with high accuracy and stable field structure
- TrEMBL: Automatically annotated entries, including a large number of unreviewed records <br>

UniProt data is widely used in fields such as data analysis, drug discovery, and protein structure prediction.


## 🔗 Official Documentation
- [UniProt API Documentation](https://www.uniprot.org/api-documentation/uniprotkb#operations-UniProtKB-getByAccession)
- [Text Search|Query syntax|Query builder](https://www.uniprot.org/help/text-search)<br>

This document outlines how to query, download, and parse Uniprot protein data using the UniProt REST API. <br> 

## 📌 What the API is for <br>
The UniProt API enables programmatic access to:
- Full protein entry details (function, structure, taxonomy, sequence, etc.)
- Complex searches with filters
- Large-scale batch downloads
- Data export in multiple formats: json, tsv, fasta, txt

## 📥 Key API Endpoints

| Endpoint | Descriptions | 
|------|------|
| /uniprotkb/{accession} | Get a single protein entry by accession ID   | 
| /uniprotkb/search | Search UniProt entries using flexible queries  | 
| /uniprotkb/stream | Stream and download large datasets efficiently   | 

## URL Parameters 
- query (str): UniProt research string.
- format (str): JSON, FASTA, etc.
- size (int): Number of results to return (e.g. size = 500)

## Code Usage 
### 1. Get Protein Information by Accession ID <br>

```python
import requests

accession = "P68871"  # Hemoglobin subunit beta
url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"

response = requests.get(url)
data = response.json()

print(data)
```
</pre>

### 2. Search Reviewed Human Proteins <br>

```python
import requests

query = "organism_id:9606 AND reviewed:true"
url = "https://rest.uniprot.org/uniprotkb/search"
params = {"query": query, "format": "json", "size": 100}

response = requests.get(url, params=params)
results = response.json().get("results", [])

for entry in results:
    print(entry["primaryAccession"])
```
</pre>

### 3. Batch Download Reviewed Human Proteins <br> 

```python
import requests

url = "https://rest.uniprot.org/uniprotkb/stream"
params = { "query": "organism_id:9606 AND reviewed:true", "format": "json", "size": 500 }

with requests.get(url, params=params, stream=True) as r:
    r.raise_for_status()
    with open("uniprot_bulk.json", "wb") as f:
        for chunk in r.iter_content(chunk_size = 500):
            f.write(chunk)

print("Download Complete")
```
</pre>

## Supported Query Fields
- organism_id: NCBI taxonomy ID (e.g. 9606 for Homo sapiens)
- reviewed:true: Swiss-Prot manually reviewed entries
- gene: gene name search (e.g. TP53)

## Advanced Queries 
- All reviewed proteins related to p53 gene: gene:TP53 AND reviewed:true
- Proteins with length > 300: length:[300 TO *]
- Newly added human proteins: organism_id:9606 AND first_public:[2024-04-01 TO 2025-05-01]

## Notes 
- Use .jsonl or flat .json for easy comparison and merging.
- Combine custom pipeline for bulk + incremental updates.




