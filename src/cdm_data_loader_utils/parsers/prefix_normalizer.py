from __future__ import annotations
from typing import Optional, Dict


SYNONYM_MAP = {
    "geneid": "ncbigene",
    "unipathway": "upa",
    "ctd": "ctd.gene",
    "gramene": "gramene.gene",
}

MAP_NAMESPACE = {
    "merops": "merops.entry",
    "ensemblbacteria": "ensembl",
    "ensemblmetazoa": "ensembl",
    "ensemblplants": "ensembl",
    "panther": "panther.family",
    "pro": "pr",
    "oma": "oma.protein",
    "paxdb": "paxdb.protein",
    "pir": "pirsf",
    "peptideatlas": "peptideatlas.peptide",
    "proteomicsdb": "proteomicsdb.protein",
    "proteomes": "uniprot.proteome",
}

ANNOTATION_SOURCE = {
    "expressionatlas",
    "funcoup",
    "glycosmos",
    "glygen",
    "inparanoid",
    "iptmnet",
    "metosite",
    "phosphositeplus",
    "smr",
    "swisspalm",
    "topdownproteomics",
}

INTERNAL_METADATA = {
    "gene_name",
    "gene_orfname",
    "gene_orderedlocusname",
    "crc64",
    "uniprotkb-id",
    "ensemblgenome_pro",
    "ensemblgenome_trs",
    "ensemblgenome",
}

# Registry gaps identified during audit
REGISTRY_GAP = {
    "collectf",
    "alphafolddb",
    "agr",
    "antibodypedia",
    "bgee",
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
}


def normalize_prefix(
    db: Optional[str],
    registry_set: set[str],
) -> Dict[str, Optional[str]]:
    if db is None:
        return {"normalized": None, "category": "null"}

    db = db.strip().lower()

    if not db:
        return {"normalized": None, "category": "null"}

    if db in INTERNAL_METADATA:
        return {"normalized": None, "category": "internal"}

    if db in ANNOTATION_SOURCE:
        return {"normalized": db, "category": "annotation"}

    if db in SYNONYM_MAP:
        return {"normalized": SYNONYM_MAP[db], "category": "synonym"}

    if db in MAP_NAMESPACE:
        return {"normalized": MAP_NAMESPACE[db], "category": "map"}

    if db in registry_set:
        return {"normalized": db, "category": "exact"}

    if db in REGISTRY_GAP:
        return {"normalized": db, "category": "registry_gap"}

    return {"normalized": db, "category": "unknown"}
