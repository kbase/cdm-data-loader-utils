"""Read in the output of bbmap's `stats.sh` command."""

import json
from pathlib import Path
from typing import Any


def get_bbmap_stats(stats_file: Path) -> dict[str, dict[str, Any]]:
    """Parse the output of bbmap's `stats.sh` command.

    :param results_dir: path to `stats.json`, the output of bbmap's stats.sh command
    :type results_dir: Path
    :return: dictionary, indexed by file path and name, of contigset stats
    :rtype: dict[str, Any]
    """
    try:
        # the file is invalid json, so requires special parsing
        stats_text = stats_file.read_text()
    except Exception as err:
        err_msg = f"error parsing stats_file: {err!s}"
        raise RuntimeError(err_msg) from err

    entry_stats = None
    stats_text_edited = stats_text.replace("}\n{", "},\n{")
    try:
        entry_stats = json.loads("[" + stats_text_edited + "]")
    except Exception as err:
        err_msg = f"error parsing stats_file: {err!s}"
        raise RuntimeError(err_msg) from err

    if not isinstance(entry_stats, list):
        err_msg = "stats_file is not in the correct format"
        raise TypeError(err_msg)

    invalid_entries = [entry for entry in entry_stats if not isinstance(entry, dict)]
    if invalid_entries:
        err_msg = "stats_file is not in the correct format: stats must be in dictionary form"
        raise TypeError(err_msg)

    # Parse the output of BBMap's stats.sh."""
    contigset_stats = {}
    err_list = []
    for entry in entry_stats:
        # minimal error checking
        if "filename" not in entry:
            err_list.append(f"invalid entry format: {json.dumps(entry, indent=None)}")
            continue
        # index by file name
        filename = Path(entry["filename"])
        contigset_stats[filename.name] = entry

    if err_list:
        err_msg = "\n".join(["format error:", *err_list])
        raise RuntimeError(err_msg)

    if not contigset_stats:
        err_msg = "no valid data found in stats_file"
        raise RuntimeError(err_msg)

    return contigset_stats
