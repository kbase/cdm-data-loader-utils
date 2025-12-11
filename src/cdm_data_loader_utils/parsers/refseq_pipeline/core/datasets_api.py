import logging
from collections.abc import Iterable
from typing import Any

import requests

from cdm_data_loader_utils.parsers.refseq_pipeline.core.config import NCBI_BASE_V2

# -------------------------------
# Logging setup
# -------------------------------
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# -------------------------------
# Shared requests session
# -------------------------------
_session = None


def get_session():
    """
    Return shared requests.Session with retries for stability.
    """
    global _session
    if _session is None:
        from requests.adapters import HTTPAdapter, Retry

        s = requests.Session()
        retries = Retry(
            total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504], allowed_methods=["GET"]
        )
        adapter = HTTPAdapter(max_retries=retries, pool_connections=16, pool_maxsize=16)
        s.mount("http://", adapter)
        s.mount("https://", adapter)
        _session = s
    return _session


# -------------------------------
# NCBI Datasets API fetcher
# -------------------------------
def fetch_reports_by_taxon(
    taxon: str,
    *,
    api_key: str | None = None,
    page_size: int = 500,
    refseq_only: bool = True,
    current_only: bool = True,
    debug: bool = False,
    max_pages: int | None = None,
) -> Iterable[dict[str, Any]]:
    """
    Generator that yields genome dataset reports for a given NCBI taxon ID.

    Args:
        taxon (str): NCBI Taxonomy ID
        api_key (str): Optional NCBI API key
        page_size (int): Number of results per page
        refseq_only (bool): Filter to RefSeq assemblies
        current_only (bool): Filter to current version assemblies
        debug (bool): Print debug info during fetch
        max_pages (int): Optional limit for number of pages to fetch

    Yields:
        dict: Single genome dataset report
    """
    url = f"{NCBI_BASE_V2}/genome/taxon/{taxon}/dataset_report"

    params = {
        "page_size": page_size,
        "returned_content": "COMPLETE",
        "filters.report_type": "assembly_report",
    }
    if current_only:
        params["filters.assembly_version"] = "current"
    if refseq_only:
        params["filters.assembly_source"] = "refseq"

    headers = {"Accept": "application/json"}
    if api_key:
        headers["api-key"] = api_key

    session = get_session()
    token = None
    page_idx = 0
    error_count = 0
    MAX_ERRORS = 3

    while True:
        if token:
            params["page_token"] = token

        try:
            resp = session.get(url, params=params, headers=headers, timeout=60)
            resp.raise_for_status()
            error_count = 0  # reset on success
        except requests.RequestException as e:
            logger.warning(f"[datasets] Error fetching taxon={taxon} page={page_idx}: {e}")
            error_count += 1
            if error_count >= MAX_ERRORS:
                logger.error(f"[datasets] Exceeded max retries ({MAX_ERRORS}); aborting.")
                break
            continue

        payload = resp.json()
        reports = payload.get("reports", []) or []

        if debug:
            logger.info(f"[datasets] taxon={taxon} page={page_idx} reports={len(reports)}")

        if not reports:
            logger.warning(f"[datasets] taxon={taxon} page={page_idx} returned no reports.")
            break

        for r in reports:
            if debug:
                acc = r.get("current_accession") or r.get("accession")
                logger.info(f"[datasets]  - accession={acc}")
            yield r

        token = payload.get("next_page_token")
        page_idx += 1

        if not token:
            logger.info(f"[datasets] taxon={taxon} finished pagination at page={page_idx}")
            break
        if max_pages is not None and page_idx >= max_pages:
            logger.info(f"[datasets] Reached max_pages={max_pages}; stopping.")
            break
