# UniProt Incremental Import Pipeline
This document outlines the step-by-step workflow of the UniProt incremental data pipeline based on current API capabilities and your local processing logic. <br> 

## Objective 
To simulate incremental (changeset) imports from UniProt by locally filtering and merging newly downloaded entries with existing bulk data. <br> 

## Pipeline Steps 
### Step 1: Download from UniProt REST API 
- Query: organism_id:9606 AND reviewed:true
- API Endpoint: https://rest.uniprot.org/uniprotkb/search
- Save raw JSON to: data/uniprot_raw.json

### Step 2: Flatten JSON Records
- Parse nested JSON into flattened structure
- Each row: Single protein entry
- Save as:
  - JSON: data/uniprot_flat.json
  - JSON Lines: data/uniprot_flat.jsonl

### Step 3: Local Incremental Filtering
- Use field entryAudit.lastAnnotationUpdateDate
- Keep only records with update date ≥ from_date
- Save filtered results to: data/uniprot_merged_incremental_incremental.jsonl

### Step 4: Merge with Existing Dataset
- Load existing file (if exists): data/uniprot_merged.jsonl
- Append new entries
- Sort by entryAudit.lastAnnotationUpdateDate
- Deduplicate by primaryAccession
- Save final merged output to: data/uniprot_merged.jsonl

## Run Example 
```python
Running INCREMENTAL import from 2024-04-01 to 2024-05-01 
[Requesting] https://rest.uniprot.org/uniprotkb/search
[Params] {'query': 'organism_id:9606 AND reviewed:true', 'format': 'json', 'size': 100}
Raw JSON saved to: data/uniprot_raw.json
Flattened JSONL saved to: data/uniprot_merged_incremental_incremental.jsonl
[Filtered] Keeping 100 updated rows since 2024-04-01
No existing data found. Starting fresh.
Combined total: 100 rows
Merged file saved to: data/uniprot_merged.jsonl
```
</pre>

## Note 
1. The current UniProt API does not support direct filtering by date.
```python
[ERROR] Failed to retrieve data: 400 Client Error
...
'publication_date' is not a valid search field
```
</pre>

The official documentation for UniProt's REST API describes the supported fields as:<br>
- organism_id
- reviewed
- gene
- accession
- length
- cc_*, ft_* (for annotation queries)

but not supported: <br>
- publication_date
- first_public
- last_modified
- structured fields like entryAudit. <br>

2. The current pipeline simulates incremental behavior using local filters.
3. Adjust from_date in CLI to change the cutoff for updates.

## Description <br> 
first_public is used to query entries published in uniprotkb. However, it only contains the first appearance of proteins and does not contain content updates for existing entries. 
Display: first_public: [yyyy-mm-dd TO yyyy-mm-dd] <br> 
For the official site maintenance changelog in the Uniprot Release Notes. When the version is updated, it shows how many data entries have been added, modified, deleted or merged.

The JSONL format is more suitable than traditional JSON for batch import and data update (especially incremental merge) scenarios. <br>
Reason: 
1. Line-by-line processing, no need to load as a whole
	- JSONL can be processed line by line, which is suitable for comparing the difference between old and new data.
	- JSON is a large array, which must be read into memory and then processed, which is inefficient and error-prone. <br>
2. Convenient de-duplication + merging
Use pandas.read_json(... , lines=True) to load JSONL, you can easily complete the incremental merging and de-weighting.
```python
pd.concat([old_df, new_df]).drop_duplicates(subset=“primaryAccession”)
```
</pre> 
3. JSON needs to be parsed into a list first, and then transferred to DataFrame, the process is more complicated.<br>
4. Each record of JSONL is naturally a row of the database/table, which is very suitable for data import, cleaning and version control.

## Next steps 
- Schedule regular incremental runs (e.g. weekly)
- Monitor merged file size and schema consistency
- Extend to export to databases (e.g. PostgreSQL, SQLite)



