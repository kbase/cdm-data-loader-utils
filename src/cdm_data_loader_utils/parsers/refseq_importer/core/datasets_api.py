from typing import Optional
import requests


def fetch_reports_by_taxon(
    taxon: str,
    api_key: Optional[str] = None,
    page_size: int = 500,
    refseq_only: bool = True,
    current_only: bool = True,
    debug: bool = False,
):
    """
    Generator to iterate through genome dataset reports from NCBI Datasets v2 API by TaxID.

    Features:
    - Calls the NCBI Datasets v2 endpoint for genome reports.
    - Applies filters: RefSeq only / current assemblies only.
    - Handles pagination via `next_page_token`.
    - Yields report dicts for each assembly.
    
    """

    # ---------------- API endpoint ----------------
    # Base URL for NCBI Datasets REST API v2
    base = "https://api.ncbi.nlm.nih.gov/datasets/v2"
    url = f"{base}/genome/taxon/{taxon}/dataset_report"

    # ---------------- Request params ----------------
    # metadata + assembly report text
    params = {
        "page_size": page_size,
        "returned_content": "COMPLETE",
        "filters.report_type": "assembly_report"
    }

    if current_only:
        params["filters.assembly_version"] = "current"
    if refseq_only:
        params["filters.assembly_source"] = "refseq"

    # ---------------- Headers ----------------
    headers = {"Accept": "application/json"}
    if api_key:
        headers["api-key"] = api_key

    # ---------------- Pagination loop ----------------
    token = None
    while True:
        if token:
            params["page_token"] = token

        # ---- request ----
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=60)
            resp.raise_for_status()
            payload = resp.json()
        except (requests.RequestException, ValueError) as e:
            print(f"Request failed for taxon {taxon}: {e}")
            break

        # ---- Extract reports ----
        reports = payload.get("reports", [])
        if not reports:
            print(f"No reports returned for taxon {taxon}")
            break

        # ---------------- Filter loop ----------------
        for rep in reports:
            info = rep.get("assemblyInfo") or rep.get("assembly_info") or {}
            src_db = info.get("sourceDatabase")

            # Skip if explicitly marked as GenBank
            if src_db and src_db != "SOURCE_DATABASE_REFSEQ":
                continue

            # Print source info if debugging
            if debug:
                if src_db is None:
                    print(f"[DEBUG] accession={rep.get('accession')} has no sourceDatabase field")
                else:
                    print(f"[DEBUG] accession={rep.get('accession')} sourceDatabase={src_db}")

                # Print the first 200 chars of assemblyReport for inspection
                if info.get("assemblyReport"):
                    snippet = info["assemblyReport"][:200].replace("\n", " ")
                    print(f"[DEBUG] {snippet}")

            # Yield one assembly report to caller
            yield rep

        # ---------------- Handle pagination ----------------
        token = payload.get("next_page_token")
        if not token:
            break


