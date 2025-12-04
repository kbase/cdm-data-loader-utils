import uuid
from pyspark.sql.types import StructType, StructField, StringType, LongType, DoubleType

# CDM_NAMESPACE = uuid.uuid4()
CDM_NAMESPACE = uuid.UUID("ff56e4aa-b9e8-4bcf-a43b-5b1d31394b9b")
NCBI_BASE_V2 = "https://api.ncbi.nlm.nih.gov/datasets/v2"

REFSEQ_ASSEMBLY_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt"
REFSEQ_ASM_REPORTS = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"

EXPECTED_COLS = [
    "cdm_id",
    "n_contigs",
    "contig_n50",
    "contig_l50",
    "n_scaffolds",
    "scaffold_n50",
    "scaffold_l50",
    "n_component_sequences",
    "gc_percent",
    "n_chromosomes",
    "contig_bp",
    "checkm_completeness",
    "checkm_contamination",
    "checkm_version",
]

DEFAULT_HASH_TABLE = "assembly_hashes"

CDM_SCHEMA = StructType(
    [
        StructField("cdm_id", StringType(), True),
        StructField("n_contigs", LongType(), True),
        StructField("contig_n50", LongType(), True),
        StructField("contig_l50", LongType(), True),
        StructField("n_scaffolds", LongType(), True),
        StructField("scaffold_n50", LongType(), True),
        StructField("scaffold_l50", LongType(), True),
        StructField("n_component_sequences", LongType(), True),
        StructField("gc_percent", DoubleType(), True),
        StructField("n_chromosomes", DoubleType(), True),
        StructField("contig_bp", LongType(), True),
        StructField("checkm_completeness", DoubleType(), True),
        StructField("checkm_contamination", DoubleType(), True),
        StructField("checkm_version", StringType(), True),
    ]
)
