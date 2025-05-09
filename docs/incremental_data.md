UniProt's REST API does not support publication_date as a query field. 
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

