"""Constants, enums, classes, and globals used in PDB piplines."""

import re
from dataclasses import dataclass
from enum import StrEnum, auto
from pathlib import Path, PurePosixPath

DEFAULT_DESTINATION_BUCKET: PurePosixPath = PurePosixPath("cdm-lake")
DEFAULT_DESTINATION_PREFIX: PurePosixPath = PurePosixPath("tenant-general-warehouse") / "refdata" / "datasets" / "pdb"
DEFAULT_HOLDINGS_SNAPSHOT: PurePosixPath = PurePosixPath("current_holdings_snapshot.json.gz")

HOLDINGS_BASE_URL = PurePosixPath("files-beta.rcsb.org") / "pub" / "wwpdb" / "pdb" / "holdings"


class HoldingsFileTypes(StrEnum):
    """Types of holdings files used by PDB."""

    CURRENT = auto()
    LAST_MODIFIED = auto()
    REMOVED = auto()


class HoldingsFileSchemas(StrEnum):
    """Schema for extracted holdings file data."""

    ID_ONLY = auto()
    ID_DATE = auto()


@dataclass(frozen=True)
class HoldingsFile:
    """Filename and schema for a holdings file."""

    filename: Path
    schema: HoldingsFileSchemas


@dataclass(frozen=True)
class ManifestData:
    """Data needed to build the manifest files."""

    new: list[str]
    updated: list[str]
    removed: list[str]
    missing_dates: list[str]


HOLDINGS_FILES: dict[HoldingsFileTypes, HoldingsFile] = {
    HoldingsFileTypes.CURRENT: HoldingsFile(
        filename=Path("current_file_holdings.json.gz"),
        schema=HoldingsFileSchemas.ID_ONLY,
    ),
    HoldingsFileTypes.LAST_MODIFIED: HoldingsFile(
        filename=Path("released_structures_last_modified_dates.json.gz"),
        schema=HoldingsFileSchemas.ID_DATE,
    ),
    HoldingsFileTypes.REMOVED: HoldingsFile(
        filename=Path("all_removed_entries.json.gz"),
        schema=HoldingsFileSchemas.ID_ONLY,
    ),
}

# Extended PDB ID: "pdb_" followed by exactly 8 lower-case alphanumeric characters.
# Classic PDB IDs are [0-9A-Z]{4}; in extended format they are zero-padded to 8
# chars and lowercased, giving [0-9a-z]{8}
PDB_ID_RE = re.compile(r"^pdb_[0-9a-z]{8}$", re.IGNORECASE)
PDB_ID_SEARCH_RE = re.compile(r"pdb_[0-9a-z]{8}", re.IGNORECASE)


@dataclass
class PDBRecord:
    """A single PDB entry as represented in the holdings inventory.

    :param id: extended PDB ID, e.g., ``"pdb_00001abc"``
    :param last_modified: ISO-8601 date of last modification, e.g. ``"2024-01-15"``
    """

    id: str
    last_modified: str = ""
