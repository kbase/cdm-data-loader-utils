SOP: RefSeq Annotation Parser CDM Pipeline

Prefix Handling Simplified:
Refactored apply_prefix() to use a centralized PREFIX_TRANSLATION mapping with clearer conditional logic, reducing redundancy.

ContigCollection Loader Improved:
load_contig_collections() now uses schema-aware construction with row = [None] * len(schema), ensuring future-proof field alignment.

Delta Table Writer Unified:
write_to_table() now detects special schema (CONTIG_COLLECTION_MIN_SCHEMA) conditionally and uses consistent DataFrame creation and writing logic.

Cleaner SQL Preview Logic:
run_sql_query() uses clearer formatting and error handling when querying each CDM table.

Safer Type Conversion:
Added to_int() utility function with graceful fallback to None on failure.



