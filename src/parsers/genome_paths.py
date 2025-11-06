"""Parse the genome_paths_file."""

import json
from pathlib import Path
from typing import Any

VALID_FILE_TYPES = ["gff", "fna", "protein"]


def get_genome_paths(genome_paths_file: Path) -> dict[str, dict[str, Any]]:
    """Read the genome paths file and retrieve the list of files from it.

    :param genome_paths_file: path to the genome paths JSON file
    :type genome_paths_file: Path
    :return: dictionary, indexed by contigset ID, of file types and paths
    :rtype: dict[str, dict[str, Any]]
    """
    try:
        with genome_paths_file.open() as json_file:
            genome_paths = json.load(json_file)
    except Exception as err:
        err_msg = f"error parsing genome_paths_file: {err!s}"
        raise RuntimeError(err_msg) from err

    if not isinstance(genome_paths, dict):
        err_msg = "genome_paths_file is not in the correct format"
        raise TypeError(err_msg)

    if not genome_paths:
        err_msg = "no valid data found in genome_paths_file"
        raise RuntimeError(err_msg)

    entries = {}
    err_list = []
    entry_list = sorted(genome_paths.keys())
    for entry_id in entry_list:
        entry_data = genome_paths[entry_id]
        if not entry_id:
            err_list.append(f"No ID specified for entry {json.dumps(entry_data, indent=None, sort_keys=True)}")
            continue

        if not isinstance(entry_data, dict):
            err_list.append(f"{entry_id}: invalid entry format")
            continue

        if not entry_data:
            err_list.append(f"{entry_id}: no valid file types or paths found")
            continue

        invalid_keys = [k for k in entry_data if k not in VALID_FILE_TYPES]
        if invalid_keys:
            err_list.append(f"{entry_id}: invalid keys: {', '.join(sorted(invalid_keys))}")
            continue

        # entries will be unique so don't have to check whether entry_id already exists
        entries[entry_id] = entry_data

    if not entries:
        err_list.append("No valid entries found in genome_paths_file")

    if err_list:
        err_msg = "Please ensure that the genome_paths_file is in the correct format.\n\n"
        raise RuntimeError(err_msg + "\n".join(err_list))

    return genome_paths
