import re
from typing import Any

# ---------------- Robust extractors set up ----------------

# regex patterns
PAT_BIOSAMPLE = re.compile(r"\bSAMN\d+\b")
PAT_BIOPROJECT = re.compile(r"\bPRJNA\d+\b")
PAT_GCF = re.compile(r"\bGCF_\d{9}\.\d+\b")
PAT_GCA = re.compile(r"\bGCA_\d{9}\.\d+\b")


def _coalesce(*vals: Any) -> str | None:
    """
    Return the first non-empty, non-whitespace string from a list of inputs.
    """
    for v in vals:
        if isinstance(v, str):
            trimmed = v.strip()
            if trimmed:
                return trimmed
    return None


def _deep_find_str(obj: Any, target_keys: set[str]) -> str | None:
    """
    Recursively search a nested dict/list structure for the first non-empty string value under any of the target_keys.
    """
    # If the object is a dictionary
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(k, str) and k in target_keys and isinstance(v, str) and v.strip():
                return v.strip()

            # recursively search the v
            found = _deep_find_str(v, target_keys)
            if found:
                return found

    # If the object is a list, search each element
    elif isinstance(obj, list):
        for it in obj:
            found = _deep_find_str(it, target_keys)
            if found:
                return found
    return None


def _deep_collect_regex(obj: Any, pattern: re.Pattern) -> list[str]:
    """
    Recursively collect all unique regex matches from a nested structure

    Args:
        obj: The nested object to search (can be dict, list, or string).
        pattern: Compiled regex pattern to search for within string values.

    Returns:
        A sorted list of unique regex matches found anywhere inside the object.
    """
    results = set()  # Use a set to avoid duplicate matches

    def _walk(x):
        if isinstance(x, dict):
            # Recursively process each value in the dictionary
            for v in x.values():
                _walk(v)
        elif isinstance(x, list):
            # Recursively process each element in the list
            for v in x:
                _walk(v)
        elif isinstance(x, str):
            # Apply regex to the string and add matches to results
            for m in pattern.findall(x):
                results.add(m)

    # Start recursion
    _walk(obj)

    # Convert set to sorted list for consistent ordering
    return sorted(results)


def extract_created_date(rep: dict[str, Any], allow_genbank_date: bool = False, debug: bool = False) -> str | None:
    """
    Extract creation/release date for a genome assembly.

    Priority:
    - For RefSeq: release_date > assembly_date > submission_date
    - For GenBank: submission_date

    Returns None if no valid date is found.

    """
    # Normalize assembly info
    assem_data = rep.get("assembly_info") or rep.get("assemblyInfo") or {}
    src_db = rep.get("source_database") or assem_data.get("sourceDatabase")

    # Collect candidate dates
    candidates: dict[str, str] = {}
    for src in (assem_data, rep.get("assembly") or {}, rep):
        for key in [
            "releaseDate",
            "assemblyDate",
            "submissionDate",
            "release_date",
            "assembly_date",
            "submission_date",
        ]:
            v = src.get(key)  # safely fetch the value

            # Only accept non-empty string values
            if isinstance(v, str) and v.strip():
                # Normalize the key
                # "release_date" -> "releasedate"
                norm_key = key.lower().replace("_", "")

                # Store the cleaned value in candidates under the normalized key
                candidates[norm_key] = v.strip()

    if debug and candidates:
        print(f"[DEBUG] found candidates={candidates}, source={src_db}")

    # RefSeq: prioritize release > assembly > submission
    if src_db == "SOURCE_DATABASE_REFSEQ":
        for pref in ("releasedate", "assemblydate", "submissiondate"):
            if pref in candidates:
                return candidates[pref]

    # GenBank: fallback only submission date
    if allow_genbank_date and src_db == "SOURCE_DATABASE_GENBANK":
        if "submissiondate" in candidates:
            return candidates["submissiondate"]

    return None


def extract_assembly_name(rep: dict[str, Any]) -> str | None:
    """
    Extract the assembly name from a genome report record.
    Supports:
      - assemblyInfo.assemblyName
      - assembly.assemblyName
      - assembly.display_name / displayName
      - rep.assemblyName
      - deep search fallback
    """
    assembly_info = rep.get("assemblyInfo") or {}
    a = rep.get("assembly") or {}

    # Try direct fields including display_name
    v = _coalesce(
        assembly_info.get("assemblyName"),
        a.get("assemblyName"),
        a.get("display_name"),
        a.get("displayName"),
        rep.get("assemblyName"),
    )
    if v:
        return v.strip()

    # Fallback deep search for assemblyName or display_name
    return _deep_find_str(rep, {"assemblyName", "assembly_name", "display_name", "displayName"})


def extract_organism_name(rep: dict[str, Any]) -> str | None:
    assembly_info = rep.get("assemblyInfo") or {}
    a = rep.get("assembly") or {}
    org_top = rep.get("organism") or {}

    candidates = [
        org_top.get("organismName"),
        org_top.get("scientificName"),
        org_top.get("taxName"),
        org_top.get("name"),
        (assembly_info.get("organism") or {}).get("organismName")
        if isinstance(assembly_info.get("organism"), dict)
        else None,
        (a.get("organism") or {}).get("organismName") if isinstance(a.get("organism"), dict) else None,
    ]

    for v in candidates:
        if isinstance(v, str) and v.strip():
            return v.strip()

    return _deep_find_str(
        rep, {"organismName", "scientificName", "sciName", "taxName", "displayName", "organism_name", "name"}
    )


def extract_taxid(rep: dict[str, Any]) -> str | None:
    """
    Extract the NCBI Taxonomy ID (taxid) from a genome assembly report.
    """

    def is_numeric_id(value) -> bool:
        """Check if a value is a valid numeric taxid."""
        return isinstance(value, (int, float)) or (isinstance(value, str) and value.strip().isdigit())

    # --- Try top-level organism block ---
    org_top = rep.get("organism") or {}
    v = org_top.get("taxId") or org_top.get("taxid") or org_top.get("taxID")
    if is_numeric_id(v):
        return str(int(v))

    # --- Recursive search ---
    def _deep_find_taxid(x):
        if isinstance(x, dict):
            for k, v in x.items():
                # normalize key (lowercase) and check if it contains 'taxid'
                if isinstance(k, str) and k.lower().replace("_", "") in {"taxid", "taxidvalue"}:
                    if is_numeric_id(v):
                        return str(int(v))
                found = _deep_find_taxid(v)
                if found:
                    return found
        elif isinstance(x, list):
            for it in x:
                found = _deep_find_taxid(it)
                if found:
                    return found
        return None

    return _deep_find_taxid(rep)


def extract_biosample_ids(rep: dict[str, Any]) -> list[str]:
    accs = set()

    direct_bs = rep.get("biosample")
    if isinstance(direct_bs, list) and all(isinstance(x, str) for x in direct_bs):
        accs.update([x.strip() for x in direct_bs if x.strip()])
        return sorted(accs)

    for path in [
        rep.get("assemblyInfo", {}).get("biosample"),
        rep.get("assembly", {}).get("biosample"),
        rep.get("biosample"),
    ]:
        if isinstance(path, dict):
            v = path.get("accession") or path.get("biosampleAccession")
            if isinstance(v, str) and v.strip():
                accs.add(v.strip())
        elif isinstance(path, list):
            for it in path:
                if isinstance(it, dict):
                    v = it.get("accession") or it.get("biosampleAccession")
                    if isinstance(v, str) and v.strip():
                        accs.add(v.strip())

    if not accs:
        accs.update(_deep_collect_regex(rep, PAT_BIOSAMPLE))

    return sorted(accs)


def extract_bioproject_ids(rep: dict[str, Any]) -> list[str]:
    accs = set()

    direct_bp = rep.get("bioproject")
    if isinstance(direct_bp, list) and all(isinstance(x, str) for x in direct_bp):
        accs.update([x.strip() for x in direct_bp if x.strip()])
        return sorted(accs)

    for path in [
        rep.get("assemblyInfo", {}).get("bioproject"),
        rep.get("assembly", {}).get("bioproject"),
        rep.get("bioproject"),
    ]:
        if isinstance(path, dict):
            v = path.get("accession") or path.get("bioprojectAccession")
            if isinstance(v, str) and v.strip():
                accs.add(v.strip())
        elif isinstance(path, list):
            for it in path:
                if isinstance(it, dict):
                    v = it.get("accession") or it.get("bioprojectAccession")
                    if isinstance(v, str) and v.strip():
                        accs.add(v.strip())

    if not accs:
        accs.update(_deep_collect_regex(rep, PAT_BIOPROJECT))

    return sorted(accs)


def extract_assembly_accessions(rep: dict[str, Any]) -> tuple[list[str], list[str]]:
    gcf, gca = set(), set()

    def _add(acc):
        if not isinstance(acc, str):
            return
        acc = acc.strip()
        if acc.startswith("GCF_"):
            gcf.add(acc)
        elif acc.startswith("GCA_"):
            gca.add(acc)

    asm = rep.get("assembly", {})
    if isinstance(asm, dict):
        for k in ["assembly_accession", "insdc_assembly_accession"]:
            val = asm.get(k)
            if isinstance(val, list):
                for x in val:
                    _add(x)

    _add(rep.get("accession"))
    ai = rep.get("assembly_info") or rep.get("assemblyInfo") or {}
    paired = ai.get("paired_assembly", {})
    if isinstance(paired, dict):
        _add(paired.get("accession"))

    return sorted(gcf), sorted(gca)
