"""
python3 src/cdm_data_loader_utils/parsers/prefix_governance_report.py
"""

from pyspark.sql import SparkSession
from pyspark.sql.functions import lower, col
from prefix_normalizer import normalize_prefix
from collections import defaultdict
import json


spark = SparkSession.builder.appName("PrefixGovernanceReport").getOrCreate()


# Load registry
with open("registry.json") as f:
    registry = json.load(f)

registry_set = {k.lower() for k in registry.keys()}

print("panther in registry_set", "panther" in registry_set)
print("Registry size:", len(registry_set))
print([x for x in registry_set if "panther" in x])


# Load parquet
df = spark.read.parquet("dataset/part-00000-0a0d0261-1fee-477d-90d8-1df048058fbf-c000.snappy.parquet")
prefixes = df.select(lower(col("db")).alias("db")).distinct().collect()


# Run normalization
results = []
category_buckets = defaultdict(list)

for row in prefixes:
    db = row["db"]
    result = normalize_prefix(db, registry_set)

    category = result["category"]
    normalized = result["normalized"]

    results.append((db, category, normalized))
    category_buckets[category].append(db)


# Print detailed table
print("\n=== Detailed Prefix Classification ===\n")

for db, category, normalized in sorted(results):
    print(f"{db:20} | {category:15} | {normalized}")


# Print summary statistics
print("\n=== Category Summary ===\n")

priority_order = [
    "unknown",
    "registry_gap",
    "internal",
    "map",
    "synonym",
    "exact",
    "null",
]

for category in priority_order:
    if category in category_buckets:
        print(f"{category:15} : {len(category_buckets[category])}")


# unknown prefixes
print("\n=== Registry Fuzzy Matches for Unknown Prefixes ===\n")

for prefix in sorted(category_buckets.get("unknown", [])):
    matches = [r for r in registry_set if prefix in r]
    if matches:
        print(f"{prefix:20} → {matches}")


print("\n=== Sample rows for ensemblbacteria ===\n")

df.filter(lower(col("db")) == "ensemblbacteria").show(20, truncate=False)
df.filter(lower(col("db")) == "panther").show(10, truncate=False)
df.filter(lower(col("db")) == "proteomes").show(10, truncate=False)
df.filter(lower(col("db")) == "phosphositeplus").show(10, truncate=False)


print("\n=== Unknown Prefix Deep Classification ===\n")

needs_subtype = {"ctd", "gramene"}

for prefix in sorted(category_buckets.get("unknown", [])):
    matches = [r for r in registry_set if prefix in r]

    if prefix in needs_subtype:
        print(f"{prefix:20} → NEEDS SUBTYPE MAPPING")

    elif matches:
        print(f"{prefix:20} → HAS REGISTRY CANDIDATE: {matches}")

    else:
        print(f"{prefix:20} → NO REGISTRY MATCH (CHECK UNIPROT DBLIST)")


print("\n=== UNKNOWN PREFIX SAMPLE INSPECTION ===\n")

unknown_prefixes = sorted(category_buckets.get("unknown", []))

for prefix in unknown_prefixes:
    print(f"\n--- {prefix.upper()} ---")

    sample_df = df.filter(lower(col("db")) == prefix).select("entity_id", "db", "xref").limit(10)

    rows = sample_df.collect()

    if not rows:
        print("No rows found")
        continue

    xrefs = set(r["xref"] for r in rows if r["xref"])

    print("Sample xrefs:", list(xrefs)[:5])

    # heuristic: check if xref equals UniProt accession
    membership_like = True

    for r in rows:
        entity = r["entity_id"].split(":")[-1]
        xref = r["xref"]
        if entity != xref:
            membership_like = False
            break

    if membership_like:
        print("→ Likely ANNOTATION SOURCE (xref equals UniProt accession)")
    else:
        print("→ Likely EXTERNAL DATABASE (independent identifier)")


print("\n=== REGISTRY EXACT CHECK FOR REMAINING EXTERNAL DB ===\n")

external_candidates = [
    "biogrid-orcs",
    "dnasu",
    "esther",
    "funfam",
    "gene3d",
    "ncbifam",
    "patric",
    "sfld",
    "veupathdb",
    "wbparasite",
]

for prefix in external_candidates:
    exact_match = prefix in registry_set
    startswith_matches = [r for r in registry_set if r.startswith(prefix)]
    contains_matches = [r for r in registry_set if prefix in r]

    print(f"\n{prefix}")
    print("exact in registry:", exact_match)
    print("startswith matches:", startswith_matches)
    print("contains matches:", contains_matches[:5])


# Highlight risk categories
print("\n=== Unknown Prefixes (Require Investigation) ===\n")
for p in sorted(category_buckets.get("unknown", [])):
    print(p)

print("\n=== Registry Gaps (Not in BioRegistry) ===\n")
for p in sorted(category_buckets.get("registry_gap", [])):
    print(p)

spark.stop()
