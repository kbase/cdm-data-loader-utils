"""PDB holdings file parser.

Holdings json files apply one of two schemas:

* ID-only: dict[str, dict[Any]] keyed on PDB IDs with values of empty dicts
* ID/Date: dict[str, str] keyed on PDB IDs with ISO-8601 date string values, e.g., ``"2024-01-15"``

PDB IDs in both schemas are in "extended format" with ``"pdb_"`` followed by 8 lowercase alpha numeric
characters, e.g., ``"pdb_00001abc"``. Classic ID, which only contained 4 alphanumeric characters
after the ``"pdb_"`` prefix are lowercased and zero-padded in extended format:

``"PDB_1ABC"`` (classic) = ``"pdb_00001abc"`` (extended)
"""

import json
from typing import Any

from cdm_data_loaders.pdb.constants import HoldingsFile, HoldingsFileSchemas, PDBRecord


def parse_holdings_file(holdings: HoldingsFile, data: bytes) -> dict[str, PDBRecord]:
    """Parse a PDB holdings file and return data as a dict keyed on PDB ID.

    :param holdings: holdings file type information
    :param data: file data
    :return: parsed PDB record info as a dict keyed on PDB ID
    """
    json_data = json.loads(data)
    if not isinstance(json_data, dict):
        msg = "Holdings file expected to contain a JSON dict"
        raise TypeError(msg)
    return _parse_pdb_record(holdings, json_data)


def _parse_pdb_record(holdings: HoldingsFile, records: dict[str, Any]) -> dict[str, PDBRecord]:
    """Parse a key-value pair from a holdings file entry into a PDBRecord."""
    if holdings.schema == HoldingsFileSchemas.ID_DATE:
        # extract data as a dict of id/updated-date pairs
        return {key.lower(): PDBRecord(id=key.lower(), last_modified=val) for key, val in records.items()}
    if holdings.schema == HoldingsFileSchemas.ID_ONLY:
        # extract keys (ids) from data
        return {key.lower(): PDBRecord(id=key.lower()) for key in records}
    msg = f"Invalid holdings file schema: {holdings.schema}"
    raise ValueError(msg)
