from pyspark.sql import SparkSession
from delta.tables import DeltaTable
from minio import Minio
from minio.error import S3Error

#Janaka E

# Initialize Spark session and MinIO client
spark = get_spark_session()
minio_client = get_minio_client()

# File-to-table mapping
file_table_mapping = {
    "janaka_db-source/sp_proteins_Jan8_parquet/feature_x_protein.parquet": "feature_x_protein",
    "janaka_db-source/sp_proteins_Jan8_parquet/protein_table.parquet": "protein",
    "janaka_db-source/sp_proteins_Jan8_parquet/name_table.parquet": "name",
    "janaka_db-source/sp_proteins_Jan8_parquet/identifier_table.parquet": "identifier",
    "janaka_db-source/sp_proteins_Jan8_parquet/association_table.parquet": "association",
}

# Bucket and namespace information
bucket_name = "cdm-lake"
namespace = "janaka_db"

# Existing Delta table locations (aligned with the metastore)
existing_delta_locations = {
    "protein": "janaka_db-deltalake/protein_table_delta",
    "feature_x_protein": "janaka_db-deltalake/feature_x_protein_delta",
    "name": "janaka_db-deltalake/name_table_delta",
    "identifier": "janaka_db-deltalake/identifier_table_delta",
    "association": "janaka_db-deltalake/association_table_delta",
}

# Process each file and upload/update to its respective Delta table
for file_name, table_name in file_table_mapping.items():
    try:
        # Derive paths
        parquet_file_path = f"s3a://{bucket_name}/{file_name}"
        delta_table_path = existing_delta_locations[table_name]
        delta_table_s3_path = f"s3a://{bucket_name}/{delta_table_path}"
        spark_table = f"{namespace}.{table_name}"

        print(f"Processing file: {file_name} -> Table: {table_name}")

        # Load the Parquet file into Spark
        df_spark = spark.read.parquet(parquet_file_path)
        print(f"Loaded Parquet file from {parquet_file_path} into Spark DataFrame.")

        # Check if Delta table exists
        if DeltaTable.isDeltaTable(spark, delta_table_s3_path):
            # Perform upsert (merge) into the Delta table
            delta_table = DeltaTable.forPath(spark, delta_table_s3_path)
            # Define merge condition based on primary keys (adjust based on schema)
            if table_name == "feature_x_protein":
                merge_condition = "existing.protein_id = updates.protein_id AND existing.feature_id = updates.feature_id"
            else:
                merge_condition = "existing.protein_id = updates.protein_id"

            # Perform the merge
            delta_table.alias("existing").merge(
                df_spark.alias("updates"),
                merge_condition
            ).whenMatchedUpdateAll(
            ).whenNotMatchedInsertAll(
            ).execute()
            print(f"Table '{table_name}' updated successfully with new data.")
        else:
            # Create a new Delta table
            (df_spark.write.mode("overwrite")
                .format("delta")
                .option("path", delta_table_s3_path)
                .saveAsTable(spark_table))
            print(f"New table '{table_name}' created successfully in namespace '{namespace}'.")

        # Verify the table contents
        result = spark.sql(f"SELECT * FROM {spark_table} LIMIT 5")
        result.show(truncate=False)
        print(f"Verification successful for table: {table_name}.")

    except Exception as e:
        print(f"Error processing file {file_name} for table {table_name}: {e}")
