from pyspark.sql import SparkSession
from pyspark.sql.functions import lower, col
from prefix_normalizer import normalize_prefix
import json


spark = SparkSession.builder.appName("PrefixTest").getOrCreate()

with open("registry.json") as f:
    registry = json.load(f)

registry_set = {k.lower() for k in registry.keys()}


# Load parquet
df = spark.read.parquet("dataset/part-00000-0a0d0261-1fee-477d-90d8-1df048058fbf-c000.snappy.parquet")
prefixes = df.select(lower(col("db")).alias("db")).distinct().collect()


# Run normalization
results = []

for row in prefixes:
    db = row["db"]
    result = normalize_prefix(db, registry_set)
    results.append((db, result["category"], result["normalized"]))


for db, category, normalized in sorted(results):
    print(f"{db:20} | {category:15} | {normalized}")

spark.stop()
