from typing import Any, Dict, List, Set
from pyspark.sql import SparkSession, DataFrame
from .datasets_api import fetch_reports_by_taxon


from .extractors import (
    _coalesce,
    extract_created_date,
    extract_assembly_name,
    extract_organism_name,
    extract_assembly_accessions
)

from .cdm_builders import (
    build_cdm_entity,
    build_cdm_contig_collection,
    build_cdm_name_rows,
    build_cdm_identifier_rows,
)


# ----------------------------------------------------------
# Pure PySpark version
# ----------------------------------------------------------
def process_report(
    spark: SparkSession,
    rep: Dict[str, Any],
    tx: str,
    seen: Set[str],
    debug: bool,
    allow_genbank_date: bool
):
    from datetime import date

    entities = []
    collections = []
    names = []
    identifiers = []

    # --- choose accession ---
    gcf_list, gca_list = extract_assembly_accessions(rep)
    acc = gcf_list[0] if gcf_list else (gca_list[0] if gca_list else None)

    asm_name = extract_assembly_name(rep)
    org_name = extract_organism_name(rep)

    created = extract_created_date(rep, allow_genbank_date)
    if not created:
        if debug:
            print(f"[WARN] No RefSeq date for accession={acc}")
        created = date.today().isoformat()

    # unique key
    key = _coalesce(acc, asm_name, org_name, tx)
    if not key or key in seen:
        return entities, collections, names, identifiers
    seen.add(key)

    # --- entity ---
    df_entity, entity_id = build_cdm_entity(spark, key, created)
    entities.append(df_entity)

    # --- contig_collection ---
    df_coll = build_cdm_contig_collection(spark, entity_id, taxid=tx)
    collections.append(df_coll)

    # --- names: convert Spark DF → list[dict] ---
    df_name = build_cdm_name_rows(spark, entity_id, rep)
    if df_name.count() > 0:
        names.extend(df_name.toPandas().to_dict("records"))

    # --- identifiers: list[dict] ---
    rows_id = build_cdm_identifier_rows(entity_id, rep, tx)
    identifiers.extend(rows_id)

    return entities, collections, names, identifiers


# ----------------------------------------------------------
# Pure PySpark version
# ----------------------------------------------------------
def process_taxon(
    spark: SparkSession,
    tx: str,
    api_key: str,
    debug: bool,
    allow_genbank_date: bool,
    unique_per_taxon: bool,
    seen: Set[str]
):
    """
    Process ONE TaxID → return lists of Spark DataFrames.

    Output:
        entities:     List[Spark DF]
        collections:  List[Spark DF]
        names:        List[Spark DF]
        identifiers:  List[Spark DF]
    """

    entities: List[DataFrame] = []
    collections: List[DataFrame] = []
    names: List[DataFrame] = []
    identifiers: List[DataFrame] = []

    # ===== fetch reports from NCBI API =====
    reports = list(fetch_reports_by_taxon(taxon=tx, api_key=api_key))

    # ===== keep only latest one per taxon =====
    if unique_per_taxon and reports:
        reports.sort(
            key=lambda r: extract_created_date(r, allow_genbank_date) or "0000-00-00",
            reverse=True
        )
        reports = [reports[0]]

    # ===== parse each assembly =====
    for rep in reports:
        e, c, n, i = process_report(
            spark=spark,
            rep=rep,
            tx=tx,
            seen=seen,
            debug=debug,
            allow_genbank_date=allow_genbank_date
        )

        entities.extend(e)
        collections.extend(c)
        names.extend(n)
        identifiers.extend(i)

    return entities, collections, names, identifiers

