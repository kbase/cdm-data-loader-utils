from __future__ import annotations
from typing import Optional, Dict


SYNONYM_MAP = {
    "geneid": "ncbigene",
    "unipathway": "upa",
}

MAP_NAMESPACE = {
    "merops": "merops.entry",
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

    if db in SYNONYM_MAP:
        return {"normalized": SYNONYM_MAP[db], "category": "synonym"}

    if db in MAP_NAMESPACE:
        return {"normalized": MAP_NAMESPACE[db], "category": "map"}

    if db in registry_set:
        return {"normalized": db, "category": "exact"}

    if db in REGISTRY_GAP:
        return {"normalized": db, "category": "registry_gap"}

    return {"normalized": db, "category": "unknown"}
