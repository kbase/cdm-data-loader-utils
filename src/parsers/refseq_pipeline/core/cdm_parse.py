import uuid
import pandas as pd
from typing import Any, Dict, List, Optional, Union, TYPE_CHECKING
from refseq_pipeline.core.config import CDM_NAMESPACE, CDM_SCHEMA, EXPECTED_COLS

if TYPE_CHECKING:
    from pyspark.sql import DataFrame as SparkDF

# === Safe type conversions ===

def safe_int(v: Any) -> Optional[int]:
    try:
        return int(str(v).replace(",", "").strip()) if v not in (None, "") else None
    except Exception:
        return None

def safe_float(v: Any) -> Optional[float]:
    try:
        return float(str(v).replace(",", "").strip()) if v not in (None, "") else None
    except Exception:
        return None

def percent_to_fraction_strict(v: Any) -> Optional[float]:
    f = safe_float(v)
    return f / 100.0 if f is not None else None


# === Field access utilities ===

def get_first(d: dict, *keys, default=None):
    for k in keys:
        if isinstance(d, dict) and k in d:
            return d[k]
    return default

def pick_section(rep: dict, snake: str, camel: str) -> dict:
    return (rep.get(snake) or rep.get(camel) or {}) or {}

def _get_accession_for_cdm(report: dict) -> str:
    ai = pick_section(report, "assembly_info", "assemblyInfo")
    paired = pick_section(ai, "paired_assembly", "pairedAssembly")
    return (
        report.get("current_accession")
        or report.get("accession")
        or paired.get("accession")
        or ""
    )


# === CDM ID ===

def generate_cdm_id(accession: str, report: Dict[str, Any]) -> str:
    if not accession:
        return f"CDM:{uuid.uuid4()}"
    info = report.get("assembly_info") or report.get("assemblyInfo") or {}
    src  = info.get("sourceDatabase") or report.get("source_database") or ""
    name = info.get("assembly_name") or ""
    ver  = info.get("assembly_version") or ""
    root = accession.split()[0]
    key  = "|".join([root, str(src), str(name), str(ver)])
    return f"CDM:{uuid.uuid5(CDM_NAMESPACE, key)}"


# === Main parser ===

def parse_report_to_row(report: Dict[str, Any]) -> Dict[str, Any]:
    asm = pick_section(report, "assembly_stats", "assemblyStats")
    chk = pick_section(report, "checkm_info", "checkmInfo")
    return {
        "cdm_id": generate_cdm_id(_get_accession_for_cdm(report), report),
        "n_contigs": safe_int(get_first(asm, "number_of_contigs", "numberOfContigs")),
        "contig_n50": safe_int(get_first(asm, "contig_n50", "contigN50")),
        "contig_l50": safe_int(get_first(asm, "contig_l50", "contigL50")),
        "n_scaffolds": safe_int(get_first(asm, "number_of_scaffolds", "numberOfScaffolds")),
        "scaffold_n50": safe_int(get_first(asm, "scaffold_n50", "scaffoldN50")),
        "scaffold_l50": safe_int(get_first(asm, "scaffold_l50", "scaffoldL50")),
        "n_component_sequences": safe_int(get_first(asm, "number_of_component_sequences", "numberOfComponentSequences")),
        "gc_percent": safe_float(get_first(asm, "gc_percent", "gcPercent")),
        "n_chromosomes": safe_float(get_first(asm, "total_number_of_chromosomes", "totalNumberOfChromosomes")),
        "contig_bp": safe_int(get_first(asm, "total_sequence_length", "totalSequenceLength")),
        "checkm_completeness": percent_to_fraction_strict(get_first(chk, "completeness", "checkmCompleteness")),
        "checkm_contamination": percent_to_fraction_strict(get_first(chk, "contamination", "checkmContamination")),
        "checkm_version": get_first(chk, "checkm_version", "checkmVersion"),
    }


def parse_reports(
    reports: List[Dict[str, Any]],
    *,
    return_spark: bool = False,
    spark=None
) -> Union[pd.DataFrame, "SparkDF"]:
    """
    Parse reports into Pandas or Spark DataFrame.

    Args:
        reports: list of raw report dicts
        return_spark: if True, return Spark DataFrame instead of Pandas
        spark: SparkSession required if return_spark=True

    Returns:
        pd.DataFrame or SparkDF
    """
    rows = [parse_report_to_row(r) for r in reports]
    pdf = pd.DataFrame(rows)

    if pdf.empty:
        return pdf

    # Ensure all expected columns exist
    for col in EXPECTED_COLS:
        if col not in pdf.columns:
            pdf[col] = None
    pdf = pdf[EXPECTED_COLS]

    if return_spark:
        if spark is None:
            raise ValueError("SparkSession must be provided when return_spark=True")
        return spark.createDataFrame(pdf, schema=CDM_SCHEMA)

    return pdf
