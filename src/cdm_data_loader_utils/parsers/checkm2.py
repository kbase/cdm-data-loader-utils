"""Read in the quality report from checkm2."""

import csv
from pathlib import Path
from typing import Any

# Ensure required columns are present
REQUIRED_COLS = ["Name", "Completeness", "Contamination"]


def get_checkm2_data(tsv_file_path: Path) -> dict[str, Any]:
    """Parse the data from checkm2, if it exists."""
    checkm2_data = {}
    err_list = []
    try:
        # Open and read the TSV file
        with tsv_file_path.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                raise ValueError("file is not in TSV format")

            missing_cols = [col for col in REQUIRED_COLS if col not in reader.fieldnames]
            if missing_cols:
                err_msg = f"checkm2 output is missing the following columns: {', '.join(missing_cols)}"
                raise ValueError(err_msg)

            # Loop through each row and extract the metrics
            for row in reader:
                if not row["Name"]:
                    err_list.append(f"row {reader.line_num} has no Name value")
                    continue
                # Convert completeness and contamination to floats
                checkm2_data[row["Name"]] = {
                    "checkm2_completeness": float(row["Completeness"]) if row["Completeness"] else None,
                    "checkm2_contamination": float(row["Contamination"]) if row["Contamination"] else None,
                }

    except Exception as err:
        err_msg = f"error parsing checkm2_file: {err!s}"
        raise RuntimeError(err_msg) from err

    if err_list:
        err_msg = "\n".join(["errors found in checkm2_file:", *err_list])
        raise RuntimeError(err_msg)

    if not checkm2_data:
        err_msg = "no valid data found in checkm2_file"
        raise RuntimeError(err_msg)

    return checkm2_data
