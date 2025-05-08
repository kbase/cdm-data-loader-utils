# Uniprot API & Downloads Overview
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


