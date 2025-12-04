import uuid
from datetime import datetime, date
from typing import Any, Dict, Literal, Tuple, Optional, List
from pyspark.sql import SparkSession, DataFrame
from pyspark.sql.types import StructType, StructField, StringType


from .extractors import (
    extract_organism_name,
    extract_assembly_name,
    extract_biosample_ids,
    extract_bioproject_ids,
    extract_taxid,
    extract_assembly_accessions,
)


# ============================================================
#                      CDM DATASOURCE
# ============================================================


def build_cdm_datasource(spark: SparkSession) -> DataFrame:
    """
    Build datasource table in Spark.
    """
    schema = StructType(
        [
            StructField("name", StringType(), False),
            StructField("source", StringType(), False),
            StructField("url", StringType(), False),
            StructField("accessed", StringType(), False),
            StructField("version", StringType(), False),
        ]
    )

    row = {
        "name": "RefSeq",
        "source": "NCBI RefSeq",
        "url": "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/",
        "accessed": date.today().isoformat(),
        "version": "231",
    }

    return spark.createDataFrame([row], schema=schema)


# ============================================================
#                         ENTITY
# ============================================================

CDM_NAMESPACE = uuid.UUID("11111111-2222-3333-4444-555555555555")


def build_entity_id(key: str) -> str:
    """Generate deterministic CDM ID."""
    return f"CDM:{uuid.uuid5(CDM_NAMESPACE, (key or '').strip())}"


def build_cdm_entity(
    spark: SparkSession,
    key_for_uuid: str,
    created_date: Optional[str],
    *,
    entity_type: Literal["contig_collection", "genome", "protein", "gene"] = "contig_collection",
    data_source: str = "RefSeq",
) -> Tuple[DataFrame, str]:
    entity_id = build_entity_id(key_for_uuid)

    schema = StructType(
        [
            StructField("entity_id", StringType(), False),
            StructField("entity_type", StringType(), True),
            StructField("data_source", StringType(), True),
            StructField("created", StringType(), True),
            StructField("updated", StringType(), True),
        ]
    )

    row = {
        "entity_id": entity_id,
        "entity_type": entity_type,
        "data_source": data_source,
        "created": created_date or date.today().isoformat(),
        "updated": datetime.now().isoformat(timespec="seconds"),
    }

    return spark.createDataFrame([row], schema=schema), entity_id


# ============================================================
#                   CONTIG_COLLECTION
# ============================================================


def build_cdm_contig_collection(
    spark: SparkSession, entity_id: str, taxid: Optional[str] = None, collection_type: str = "isolate"
) -> DataFrame:
    schema = StructType(
        [
            StructField("collection_id", StringType(), False),
            StructField("contig_collection_type", StringType(), True),
            StructField("ncbi_taxon_id", StringType(), True),
            StructField("gtdb_taxon_id", StringType(), True),
        ]
    )

    row = {
        "collection_id": entity_id,
        "contig_collection_type": collection_type,
        "ncbi_taxon_id": f"NCBITaxon:{taxid}" if taxid else None,
        "gtdb_taxon_id": None,
    }

    return spark.createDataFrame([row], schema=schema)


# ============================================================
#                        NAME TABLE
# ============================================================


def build_cdm_name_rows(spark: SparkSession, entity_id: str, rep: Dict[str, Any]) -> Optional[DataFrame]:
    schema = StructType(
        [
            StructField("entity_id", StringType(), False),
            StructField("name", StringType(), False),
            StructField("description", StringType(), True),
            StructField("source", StringType(), True),
        ]
    )

    rows = []

    org = extract_organism_name(rep)
    if org:
        rows.append(
            {"entity_id": entity_id, "name": org.strip(), "description": "RefSeq organism name", "source": "RefSeq"}
        )

    asm = extract_assembly_name(rep)
    if asm:
        rows.append(
            {"entity_id": entity_id, "name": asm.strip(), "description": "RefSeq assembly name", "source": "RefSeq"}
        )

    if not rows:
        return None

    return spark.createDataFrame(rows, schema=schema)


# ============================================================
#                   IDENTIFIER TABLE
# ============================================================

IDENTIFIER_PREFIXES = {
    "biosample": ("Biosample", "BioSample ID"),
    "bioproject": ("BioProject", "BioProject ID"),
    "taxon": ("NCBITaxon", "NCBI Taxon ID"),
    "gcf": ("ncbi.assembly", "NCBI Assembly ID"),
    "gca": ("insdc.gca", "GenBank Assembly ID"),
}


def build_cdm_identifier_rows(
    entity_id: str, rep: Dict[str, Any], request_taxid: Optional[str]
) -> List[Dict[str, Any]]:
    """
    Identifiers remain as list[dict].
    Spark conversion happens later in finalize_tables.

    """

    rows = []

    # ---- BioSample IDs ----
    for bs in extract_biosample_ids(rep):
        rows.append(
            {
                "entity_id": entity_id,
                "identifier": f"{IDENTIFIER_PREFIXES['biosample'][0]}:{bs.strip()}",
                "source": "RefSeq",
                "description": IDENTIFIER_PREFIXES["biosample"][1],
            }
        )

    # ---- BioProject IDs ----
    for bp in extract_bioproject_ids(rep):
        rows.append(
            {
                "entity_id": entity_id,
                "identifier": f"{IDENTIFIER_PREFIXES['bioproject'][0]}:{bp.strip()}",
                "source": "RefSeq",
                "description": IDENTIFIER_PREFIXES["bioproject"][1],
            }
        )

    # ---- Taxon ----
    tx = extract_taxid(rep) or (str(request_taxid).strip() if request_taxid else None)
    if tx and tx.isdigit():
        rows.append(
            {
                "entity_id": entity_id,
                "identifier": f"{IDENTIFIER_PREFIXES['taxon'][0]}:{tx}",
                "source": "RefSeq",
                "description": IDENTIFIER_PREFIXES["taxon"][1],
            }
        )

    # ---- Assembly Accessions ----
    gcf_list, gca_list = extract_assembly_accessions(rep)

    for gcf in gcf_list:
        rows.append(
            {
                "entity_id": entity_id,
                "identifier": f"{IDENTIFIER_PREFIXES['gcf'][0]}:{gcf.strip()}",
                "source": "RefSeq",
                "description": IDENTIFIER_PREFIXES["gcf"][1],
            }
        )

    for gca in gca_list:
        rows.append(
            {
                "entity_id": entity_id,
                "identifier": f"{IDENTIFIER_PREFIXES['gca'][0]}:{gca.strip()}",
                "source": "RefSeq",
                "description": IDENTIFIER_PREFIXES["gca"][1],
            }
        )

    # ---- Deduplicate ----
    uniq = []
    seen = set()
    for r in rows:
        key = (r["entity_id"], r["identifier"])
        if key not in seen:
            seen.add(key)
            uniq.append(r)

    return uniq
