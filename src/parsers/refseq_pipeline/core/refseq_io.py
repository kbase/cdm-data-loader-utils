import requests
import hashlib
import logging
from typing import Optional, TypedDict, Dict
from functools import lru_cache
import pandas as pd

from refseq_pipeline.core.config import REFSEQ_ASSEMBLY_SUMMARY_URL

"""
python -m refseq_pipeline.core.refseq_io

DONE 

"""


# ----------------------------------------
# Logger setup (configurable externally)
# ----------------------------------------
logger = logging.getLogger(__name__)


# ----------------------------------------
# Typed structure for parsed metadata
# ----------------------------------------
class AssemblyMeta(TypedDict):
    ftp_path: str
    taxid: str
    species_taxid: str


# ----------------------------------------
# Shared requests session with retry
# ----------------------------------------
_session = None
def get_session() -> requests.Session:
    """
    Return a shared session with retry logic for stable downloads.
    Avoids creating a new session for every request.
    """
    global _session
    if _session is None:
        from requests.adapters import HTTPAdapter, Retry
        s = requests.Session()
        retries = Retry(
            total=3,
            backoff_factor=0.5,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"]
        )
        adapter = HTTPAdapter(max_retries=retries, pool_connections=16, pool_maxsize=16)
        s.mount("http://", adapter)
        s.mount("https://", adapter)
        _session = s
    return _session


# ----------------------------------------
# Download + normalization helpers
# ----------------------------------------
def download_text(url: str, timeout: int = 60, session: Optional[requests.Session] = None) -> str:
    """
    Download raw text from a URL using retry-enabled session.
    """
    s = session or get_session()
    r = s.get(url, timeout=timeout)
    r.raise_for_status()
    return r.text

def normalize_multiline_text(txt: str) -> str:
    """
    Clean up multiline text by trimming trailing whitespace from each line.
    """
    return "\n".join([ln.rstrip() for ln in txt.splitlines()])


# ----------------------------------------
# RefSeq Assembly Summary Parser
# ----------------------------------------
def parse_assembly_summary(content: str) -> Dict[str, AssemblyMeta]:
    """
    Parse assembly_summary_refseq.txt contents into a structured dict.
    Returns accession → metadata dictionary.
    """
    acc2meta: Dict[str, AssemblyMeta] = {}
    for line in content.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 20:
            continue
        accession = parts[0].strip()
        taxid = parts[5].strip() if len(parts) > 6 else ""
        species_taxid = parts[6].strip() if len(parts) > 7 else ""
        ftp_path = parts[19].strip().rstrip("/")
        if accession and ftp_path:
            acc2meta[accession] = {
                "ftp_path": ftp_path,
                "taxid": taxid,
                "species_taxid": species_taxid,
            }
    return acc2meta

def load_refseq_assembly_index(url: Optional[str] = None) -> Dict[str, AssemblyMeta]:
    """
    Load and parse the RefSeq assembly summary file.
    Downloads from default URL unless a custom URL is provided.
    """
    try:
        txt = download_text(url or REFSEQ_ASSEMBLY_SUMMARY_URL)
        return parse_assembly_summary(txt)
    except Exception as e:
        logger.error(f"[refseq_io] Failed to load assembly index: {e}")
        return {}

def load_local_refseq_assembly_index(path: str) -> pd.DataFrame:
    """
    Load RefSeq assembly_summary_refseq.tsv from local file.
    Skips comment lines and reads as pandas DataFrame.
    """
    try:
        df = pd.read_csv(path, sep="\t", comment="#", dtype=str)
        return df
    except Exception as e:
        logger.error(f"[refseq_io] Failed to load local TSV: {path} → {e}")
        raise

# ----------------------------------------
# Fetch file hashes from annotation FTP
# ----------------------------------------
@lru_cache(maxsize=None)
def fetch_annotation_hash(ftp_path: str, timeout: int = 30) -> Optional[str]:
    """
    Retrieve and normalize the contents of annotation_hashes.txt under a given FTP path.
    Uses LRU cache to avoid redundant network calls.
    """
    url = f"{ftp_path}/annotation_hashes.txt"
    try:
        return normalize_multiline_text(download_text(url, timeout))
    except Exception as e:
        logger.warning(f"[fetch] Failed to fetch {url}: {e}")
        return None

@lru_cache(maxsize=None)
def fetch_md5_checksums(ftp_path: str, timeout: int = 30) -> Optional[str]:
    """
    Retrieve and normalize the contents of md5checksums.txt under a given FTP path.
    Uses LRU cache to avoid redundant network calls.
    """
    url = f"{ftp_path}/md5checksums.txt"
    try:
        return normalize_multiline_text(download_text(url, timeout))
    except Exception as e:
        logger.warning(f"[fetch] Failed to fetch {url}: {e}")
        return None


# ----------------------------------------
# SHA256 hash calculator for text content
# ----------------------------------------
def text_sha256(s: str) -> str:
    """
    Compute SHA256 fingerprint of the given string.
    Used for file change tracking.
    """
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


if __name__ == "__main__":
    print("Downloading and parsing RefSeq index..")
    acc2meta = load_refseq_assembly_index()
    
    print(f"Total accessions parsed: {len(acc2meta)}")

    # Sample print
    for acc, meta in list(acc2meta.items())[:10]:  # show first 10
        print(f"- {acc}: {meta}")



