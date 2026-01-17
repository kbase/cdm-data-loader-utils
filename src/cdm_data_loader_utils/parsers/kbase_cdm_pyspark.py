"""Automated conversion of cdm_schema to PySpark."""

from pyspark.sql.types import BooleanType, DateType, FloatType, IntegerType, StringType, StructField, StructType

schema = {
    "Association": StructType(
        [
            StructField("association_id", StringType(), nullable=False),
            StructField("subject", StringType(), nullable=False),
            StructField("object", StringType(), nullable=False),
            StructField("predicate", StringType(), nullable=False),
            StructField("negated", BooleanType(), nullable=True),
            StructField("evidence_type", StringType(), nullable=True),
            StructField("primary_knowledge_source", StringType(), nullable=True),
            StructField("aggregator_knowledge_source", StringType(), nullable=True),
            StructField("annotation_date", DateType(), nullable=True),
            StructField("comments", StringType(), nullable=True),
        ]
    ),
    "Association_x_SupportingObject": StructType(
        [
            StructField("association_id", StringType(), nullable=False),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "Cluster": StructType(
        [
            StructField("cluster_id", StringType(), nullable=False),
            StructField("description", StringType(), nullable=True),
            StructField("name", StringType(), nullable=True),
            StructField("entity_type", StringType(), nullable=False),
            StructField("protocol_id", StringType(), nullable=True),
        ]
    ),
    "ClusterMember": StructType(
        [
            StructField("cluster_id", StringType(), nullable=False),
            StructField("entity_id", StringType(), nullable=False),
            StructField("is_representative", BooleanType(), nullable=True),
            StructField("is_seed", BooleanType(), nullable=True),
            StructField("score", FloatType(), nullable=True),
        ]
    ),
    "Contig": StructType(
        [
            StructField("contig_id", StringType(), nullable=False),
            StructField("hash", StringType(), nullable=True),
            StructField("gc_content", FloatType(), nullable=True),
            StructField("length", IntegerType(), nullable=True),
        ]
    ),
    "ContigCollection": StructType(
        [
            StructField("contig_collection_id", StringType(), nullable=False),
            StructField("hash", StringType(), nullable=True),
            StructField("asm_score", FloatType(), nullable=True),
            StructField("checkm_completeness", FloatType(), nullable=True),
            StructField("checkm_contamination", FloatType(), nullable=True),
            StructField("checkm_version", StringType(), nullable=True),
            StructField("contig_bp", IntegerType(), nullable=True),
            StructField("contig_collection_type", StringType(), nullable=True),
            StructField("contig_l50", IntegerType(), nullable=True),
            StructField("contig_l90", IntegerType(), nullable=True),
            StructField("contig_n50", IntegerType(), nullable=True),
            StructField("contig_n90", IntegerType(), nullable=True),
            StructField("contig_logsum", FloatType(), nullable=True),
            StructField("contig_max", IntegerType(), nullable=True),
            StructField("contig_powersum", FloatType(), nullable=True),
            StructField("gap_percent", FloatType(), nullable=True),
            StructField("gc_average", FloatType(), nullable=True),
            StructField("gc_std", FloatType(), nullable=True),
            StructField("gtdb_taxon_id", StringType(), nullable=True),
            StructField("n_chromosomes", IntegerType(), nullable=True),
            StructField("n_contigs", IntegerType(), nullable=True),
            StructField("n_scaffolds", IntegerType(), nullable=True),
            StructField("ncbi_taxon_id", StringType(), nullable=True),
            StructField("scaffold_l50", IntegerType(), nullable=True),
            StructField("scaffold_l90", IntegerType(), nullable=True),
            StructField("scaffold_n50", IntegerType(), nullable=True),
            StructField("scaffold_n90", IntegerType(), nullable=True),
            StructField("scaffold_bp", IntegerType(), nullable=True),
            StructField("scaffold_logsum", FloatType(), nullable=True),
            StructField("scaffold_maximum_length", IntegerType(), nullable=True),
            StructField("scaffold_powersum", FloatType(), nullable=True),
            StructField("scaffolds_n_over_50K", IntegerType(), nullable=True),
            StructField("scaffolds_percent_over_50K", FloatType(), nullable=True),
            StructField("scaffolds_total_length_over_50k", IntegerType(), nullable=True),
        ]
    ),
    "ContigCollection_x_EncodedFeature": StructType(
        [
            StructField("contig_collection_id", StringType(), nullable=False),
            StructField("encoded_feature_id", StringType(), nullable=False),
        ]
    ),
    "ContigCollection_x_Feature": StructType(
        [
            StructField("contig_collection_id", StringType(), nullable=False),
            StructField("feature_id", StringType(), nullable=False),
        ]
    ),
    "ContigCollection_x_Protein": StructType(
        [
            StructField("contig_collection_id", StringType(), nullable=False),
            StructField("protein_id", StringType(), nullable=False),
        ]
    ),
    "Contig_x_ContigCollection": StructType(
        [
            StructField("contig_id", StringType(), nullable=False),
            StructField("contig_collection_id", StringType(), nullable=False),
        ]
    ),
    "Contig_x_EncodedFeature": StructType(
        [
            StructField("contig_id", StringType(), nullable=False),
            StructField("encoded_feature_id", StringType(), nullable=False),
        ]
    ),
    "Contig_x_Feature": StructType(
        [
            StructField("contig_id", StringType(), nullable=False),
            StructField("feature_id", StringType(), nullable=False),
        ]
    ),
    "Contig_x_Protein": StructType(
        [
            StructField("contig_id", StringType(), nullable=False),
            StructField("protein_id", StringType(), nullable=False),
        ]
    ),
    "Contributor": StructType(
        [
            StructField("contributor_id", StringType(), nullable=False),
            StructField("contributor_type", StringType(), nullable=True),
            StructField("name", StringType(), nullable=True),
            StructField("given_name", StringType(), nullable=True),
            StructField("family_name", StringType(), nullable=True),
        ]
    ),
    "ContributorAffiliation": StructType(
        [
            StructField("contributor_id", StringType(), nullable=False),
            StructField("affiliation_id", StringType(), nullable=True),
        ]
    ),
    "Contributor_x_DataSource": StructType(
        [
            StructField("contributor_id", StringType(), nullable=False),
            StructField("data_source_id", StringType(), nullable=False),
            StructField("contributor_role", StringType(), nullable=True),
        ]
    ),
    "Contributor_x_Role_x_Project": StructType(
        [
            StructField("contributor_id", StringType(), nullable=False),
            StructField("project_id", StringType(), nullable=False),
            StructField("contributor_role", StringType(), nullable=True),
        ]
    ),
    "ControlledTermValue": StructType(
        [
            StructField("value_cv_label", StringType(), nullable=False),
            StructField("raw_value", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
            StructField("attribute_cv_id", StringType(), nullable=True),
            StructField("attribute_cv_label", StringType(), nullable=True),
            StructField("attribute_string", StringType(), nullable=True),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "ControlledVocabularyTermValue": StructType(
        [
            StructField("value_cv_label", StringType(), nullable=True),
            StructField("value_cv_id", StringType(), nullable=True),
            StructField("raw_value", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
            StructField("attribute_cv_id", StringType(), nullable=True),
            StructField("attribute_cv_label", StringType(), nullable=True),
            StructField("attribute_string", StringType(), nullable=True),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "DataSource": StructType(
        [
            StructField("data_source_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=True),
        ]
    ),
    "DataSourceNew": StructType(
        [
            StructField("data_source_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=True),
            StructField("comments", StringType(), nullable=True),
            StructField("date_accessed", DateType(), nullable=False),
            StructField("date_published", DateType(), nullable=True),
            StructField("date_updated", DateType(), nullable=True),
            StructField("license", StringType(), nullable=True),
            StructField("publisher", StringType(), nullable=True),
            StructField("resource_type", StringType(), nullable=False),
            StructField("url", StringType(), nullable=True),
            StructField("version", StringType(), nullable=True),
        ]
    ),
    "DataSource_x_Description": StructType(
        [
            StructField("data_source_id", StringType(), nullable=False),
            StructField("resource_description_id", StringType(), nullable=False),
        ]
    ),
    "DataSource_x_FundingReference": StructType(
        [
            StructField("data_source_id", StringType(), nullable=False),
            StructField("funding_reference_id", StringType(), nullable=False),
        ]
    ),
    "DataSource_x_License": StructType(
        [
            StructField("data_source_id", StringType(), nullable=False),
            StructField("license_id", StringType(), nullable=False),
        ]
    ),
    "DataSource_x_Title": StructType(
        [
            StructField("data_source_id", StringType(), nullable=False),
            StructField("resource_title_id", StringType(), nullable=False),
        ]
    ),
    "DateTimeValue": StructType(
        [
            StructField("date_time", DateType(), nullable=False),
            StructField("raw_value", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
            StructField("attribute_cv_id", StringType(), nullable=True),
            StructField("attribute_cv_label", StringType(), nullable=True),
            StructField("attribute_string", StringType(), nullable=True),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "EncodedFeature": StructType(
        [
            StructField("encoded_feature_id", StringType(), nullable=False),
            StructField("hash", StringType(), nullable=True),
            StructField("has_stop_codon", BooleanType(), nullable=True),
            StructField("type", StringType(), nullable=True),
        ]
    ),
    "EncodedFeature_x_Feature": StructType(
        [
            StructField("encoded_feature_id", StringType(), nullable=False),
            StructField("feature_id", StringType(), nullable=False),
        ]
    ),
    "EncodedFeature_x_Protein": StructType(
        [
            StructField("encoded_feature_id", StringType(), nullable=False),
            StructField("protein_id", StringType(), nullable=False),
        ]
    ),
    "EntailedEdge": StructType(
        [
            StructField("subject", StringType(), nullable=True),
            StructField("predicate", StringType(), nullable=True),
            StructField("object", StringType(), nullable=True),
        ]
    ),
    "Entity": StructType(
        [
            StructField("entity_id", StringType(), nullable=False),
            StructField("entity_type", StringType(), nullable=False),
            StructField("data_source_id", StringType(), nullable=True),
            StructField("data_source_entity_id", StringType(), nullable=True),
            StructField("data_source_created", DateType(), nullable=False),
            StructField("data_source_updated", DateType(), nullable=True),
            StructField("created", DateType(), nullable=False),
            StructField("updated", DateType(), nullable=False),
        ]
    ),
    "Event": StructType(
        [
            StructField("event_id", StringType(), nullable=False),
            StructField("created_at", DateType(), nullable=True),
            StructField("description", StringType(), nullable=True),
            StructField("name", StringType(), nullable=True),
            StructField("location", StringType(), nullable=True),
        ]
    ),
    "Experiment": StructType(
        [
            StructField("experiment_id", StringType(), nullable=False),
            StructField("protocol_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=True),
            StructField("description", StringType(), nullable=True),
            StructField("created_at", DateType(), nullable=True),
        ]
    ),
    "ExperimentCondition": StructType(
        [
            StructField("experiment_condition_id", StringType(), nullable=False),
            StructField("experiment_id", StringType(), nullable=False),
            StructField("variable_id", StringType(), nullable=False),
            StructField("value", StringType(), nullable=True),
        ]
    ),
    "ExperimentConditionSet": StructType(
        [
            StructField("experiment_condition_set_id", StringType(), nullable=False),
            StructField("experiment_condition_id", StringType(), nullable=False),
        ]
    ),
    "Feature": StructType(
        [
            StructField("feature_id", StringType(), nullable=False),
            StructField("hash", StringType(), nullable=True),
            StructField("cds_phase", StringType(), nullable=True),
            StructField("e_value", FloatType(), nullable=True),
            StructField("end", IntegerType(), nullable=True),
            StructField("p_value", FloatType(), nullable=True),
            StructField("start", IntegerType(), nullable=True),
            StructField("strand", StringType(), nullable=True),
            StructField("source_database", StringType(), nullable=True),
            StructField("protocol_id", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
        ]
    ),
    "Feature_x_Protein": StructType(
        [
            StructField("feature_id", StringType(), nullable=False),
            StructField("protein_id", StringType(), nullable=False),
        ]
    ),
    "FundingReference": StructType(
        [
            StructField("funding_reference_id", StringType(), nullable=False),
            StructField("funder", StringType(), nullable=True),
            StructField("grant_id", StringType(), nullable=True),
            StructField("grant_title", StringType(), nullable=True),
            StructField("grant_url", StringType(), nullable=True),
        ]
    ),
    "Geolocation": StructType(
        [
            StructField("latitude", FloatType(), nullable=False),
            StructField("longitude", FloatType(), nullable=False),
            StructField("raw_value", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
            StructField("attribute_cv_id", StringType(), nullable=True),
            StructField("attribute_cv_label", StringType(), nullable=True),
            StructField("attribute_string", StringType(), nullable=True),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "GoldEnvironmentalContext": StructType(
        [
            StructField("gold_environmental_context_id", StringType(), nullable=False),
            StructField("ecosystem", StringType(), nullable=True),
            StructField("ecosystem_category", StringType(), nullable=True),
            StructField("ecosystem_subtype", StringType(), nullable=True),
            StructField("ecosystem_type", StringType(), nullable=True),
            StructField("specific_ecosystem", StringType(), nullable=True),
        ]
    ),
    "Identifier": StructType(
        [
            StructField("entity_id", StringType(), nullable=False),
            StructField("identifier", StringType(), nullable=False),
            StructField("description", StringType(), nullable=True),
            StructField("source", StringType(), nullable=True),
            StructField("relationship", StringType(), nullable=True),
        ]
    ),
    "License": StructType(
        [
            StructField("license_id", StringType(), nullable=False),
            StructField("id", StringType(), nullable=True),
            StructField("name", StringType(), nullable=True),
            StructField("url", StringType(), nullable=True),
        ]
    ),
    "Measurement": StructType(
        [
            StructField("measurement_id", StringType(), nullable=False),
            StructField("measurement_set_id", StringType(), nullable=False),
            StructField("experiment_condition_set_id", StringType(), nullable=False),
            StructField("value", StringType(), nullable=True),
        ]
    ),
    "MeasurementSet": StructType(
        [
            StructField("measurement_set_id", StringType(), nullable=False),
            StructField("variable_id", StringType(), nullable=False),
            StructField("quality", StringType(), nullable=True),
            StructField("created_at", DateType(), nullable=True),
        ]
    ),
    "MixsEnvironmentalContext": StructType(
        [
            StructField("mixs_environmental_context_id", StringType(), nullable=False),
            StructField("env_broad_scale", StringType(), nullable=True),
            StructField("env_local_scale", StringType(), nullable=True),
            StructField("env_medium", StringType(), nullable=True),
        ]
    ),
    "Name": StructType(
        [
            StructField("entity_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=False),
            StructField("description", StringType(), nullable=True),
            StructField("source", StringType(), nullable=True),
        ]
    ),
    "OrderedProtocolStep": StructType(
        [
            StructField("protocol_id", StringType(), nullable=False),
            StructField("protocol_step_id", StringType(), nullable=False),
            StructField("step_index", IntegerType(), nullable=False),
        ]
    ),
    "Parameter": StructType(
        [
            StructField("parameter_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=True),
            StructField("description", StringType(), nullable=True),
            StructField("value_type", StringType(), nullable=True),
            StructField("required", BooleanType(), nullable=True),
            StructField("cardinality", StringType(), nullable=True),
            StructField("default", StringType(), nullable=True),
            StructField("parameter_type", StringType(), nullable=True),
        ]
    ),
    "Prefix": StructType(
        [
            StructField("prefix", StringType(), nullable=True),
            StructField("base", StringType(), nullable=True),
        ]
    ),
    "Project": StructType(
        [
            StructField("project_id", StringType(), nullable=False),
            StructField("description", StringType(), nullable=True),
        ]
    ),
    "Protein": StructType(
        [
            StructField("protein_id", StringType(), nullable=False),
            StructField("hash", StringType(), nullable=True),
            StructField("description", StringType(), nullable=True),
            StructField("evidence_for_existence", StringType(), nullable=True),
            StructField("length", IntegerType(), nullable=True),
            StructField("sequence", StringType(), nullable=True),
        ]
    ),
    "Protocol": StructType(
        [
            StructField("protocol_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=True),
            StructField("description", StringType(), nullable=True),
            StructField("doi", StringType(), nullable=True),
            StructField("url", StringType(), nullable=True),
            StructField("version", StringType(), nullable=True),
        ]
    ),
    "ProtocolExecution": StructType(
        [
            StructField("protocol_execution_id", StringType(), nullable=False),
            StructField("protocol_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=True),
            StructField("description", StringType(), nullable=True),
            StructField("created_at", DateType(), nullable=True),
        ]
    ),
    "ProtocolInput": StructType(
        [
            StructField("parameter_id", StringType(), nullable=False),
            StructField("protocol_input_id", StringType(), nullable=False),
            StructField("protocol_execution_id", StringType(), nullable=False),
            StructField("value", StringType(), nullable=False),
        ]
    ),
    "ProtocolInputSet": StructType(
        [
            StructField("protocol_input_id", StringType(), nullable=False),
            StructField("protocol_input_set_id", StringType(), nullable=False),
        ]
    ),
    "ProtocolOutput": StructType(
        [
            StructField("protocol_output_id", StringType(), nullable=False),
            StructField("protocol_input_set_id", StringType(), nullable=False),
            StructField("value", StringType(), nullable=False),
        ]
    ),
    "ProtocolStep": StructType(
        [
            StructField("protocol_step_id", StringType(), nullable=False),
            StructField("step", StringType(), nullable=True),
        ]
    ),
    "ProtocolVariable": StructType(
        [
            StructField("protocol_id", StringType(), nullable=False),
            StructField("variable_id", StringType(), nullable=False),
        ]
    ),
    "Publication": StructType(
        [
            StructField("publication_id", StringType(), nullable=False),
        ]
    ),
    "QuantityRangeValue": StructType(
        [
            StructField("maximum_numeric_value", FloatType(), nullable=False),
            StructField("minimum_numeric_value", FloatType(), nullable=False),
            StructField("unit_cv_id", StringType(), nullable=True),
            StructField("unit_cv_label", StringType(), nullable=True),
            StructField("unit_string", StringType(), nullable=True),
            StructField("raw_value", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
            StructField("attribute_cv_id", StringType(), nullable=True),
            StructField("attribute_cv_label", StringType(), nullable=True),
            StructField("attribute_string", StringType(), nullable=True),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "QuantityValue": StructType(
        [
            StructField("numeric_value", FloatType(), nullable=False),
            StructField("unit_cv_id", StringType(), nullable=True),
            StructField("unit_cv_label", StringType(), nullable=True),
            StructField("unit_string", StringType(), nullable=True),
            StructField("raw_value", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
            StructField("attribute_cv_id", StringType(), nullable=True),
            StructField("attribute_cv_label", StringType(), nullable=True),
            StructField("attribute_string", StringType(), nullable=True),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "ResourceDescription": StructType(
        [
            StructField("resource_description_id", StringType(), nullable=False),
            StructField("description_text", StringType(), nullable=False),
            StructField("description_type", StringType(), nullable=True),
            StructField("language", StringType(), nullable=True),
        ]
    ),
    "ResourceTitle": StructType(
        [
            StructField("resource_title_id", StringType(), nullable=False),
            StructField("language", StringType(), nullable=True),
            StructField("title", StringType(), nullable=False),
            StructField("title_type", StringType(), nullable=True),
        ]
    ),
    "Sample": StructType(
        [
            StructField("sample_id", StringType(), nullable=False),
            StructField("description", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
        ]
    ),
    "Sequence": StructType(
        [
            StructField("sequence_id", StringType(), nullable=False),
            StructField("entity_id", StringType(), nullable=False),
            StructField("type", StringType(), nullable=True),
            StructField("length", IntegerType(), nullable=True),
            StructField("checksum", StringType(), nullable=True),
        ]
    ),
    "Statement": StructType(
        [
            StructField("subject", StringType(), nullable=True),
            StructField("predicate", StringType(), nullable=True),
            StructField("object", StringType(), nullable=True),
            StructField("value", StringType(), nullable=True),
            StructField("datatype", StringType(), nullable=True),
            StructField("language", StringType(), nullable=True),
        ]
    ),
    "TextValue": StructType(
        [
            StructField("text_value", StringType(), nullable=False),
            StructField("language", StringType(), nullable=True),
            StructField("raw_value", StringType(), nullable=True),
            StructField("type", StringType(), nullable=True),
            StructField("attribute_cv_id", StringType(), nullable=True),
            StructField("attribute_cv_label", StringType(), nullable=True),
            StructField("attribute_string", StringType(), nullable=True),
            StructField("entity_id", StringType(), nullable=False),
        ]
    ),
    "Variable": StructType(
        [
            StructField("variable_id", StringType(), nullable=False),
            StructField("name", StringType(), nullable=True),
            StructField("description", StringType(), nullable=True),
            StructField("name_cv_id", StringType(), nullable=True),
            StructField("unit", StringType(), nullable=True),
            StructField("value_type", StringType(), nullable=False),
        ]
    ),
    "VariableValue": StructType(
        [
            StructField("variable_value_id", StringType(), nullable=False),
            StructField("variable_id", StringType(), nullable=False),
            StructField("value_type", StringType(), nullable=True),
        ]
    ),
}
