"""Parser for UniProt release notes."""

# The UniProt consortium European Bioinformatics Institute (EBI), SIB Swiss
# Institute of Bioinformatics and Protein Information Resource (PIR),
# is pleased to announce UniProt Knowledgebase (UniProtKB) Release
# 2025_03 (18-Jun-2025). UniProt (Universal Protein Resource) is a
# comprehensive catalog of information on proteins.

# UniProtKB Release 2025_03 consists of 253,635,358 entries (UniProtKB/Swiss-Prot:
# 573,661 entries and UniProtKB/TrEMBL: 253,061,697 entries)
# UniRef100 Release 2025_03 consists of 465,330,530 entries
# UniRef90 Release 2025_03 consists of 208,005,650 entries
# UniRef50 Release 2025_03 consists of 70,198,728 entries
# UniParc Release 2025_03 consists of 982,121,738 entries, where 915,805,719 are active and 66,316,019 inactive
# UniProt databases can be accessed from the web at http://www.uniprot.org and
# downloaded from http://www.uniprot.org/downloads. Detailed release
# statistics for TrEMBL and Swiss-Prot sections of the UniProt Knowledgebase
# can be viewed at http://www.ebi.ac.uk/uniprot/TrEMBLstats/ and
# http://web.expasy.org/docs/relnotes/relstat.html respectively.

import datetime as dt
import re
from pathlib import Path
from typing import Any

from cdm_data_loader_utils.utils.cdm_logger import log_and_die

RELEASE_VERSION_DATE: re.Pattern[str] = re.compile(
    r"is pleased to announce UniProt Knowledgebase \(UniProtKB\) Release\s+(\w+) \((\d{1,2}-[a-zA-Z]+-\d{4})\)\."
)

UNIPROT_TREMBL_STATS: re.Pattern[str] = re.compile(
    r"UniProtKB Release \w+ consists of ([\d,]+) entries \(UniProtKB/Swiss-Prot:\n([\d,]+) entries and UniProtKB/TrEMBL: ([\d,]+) entries\)"
)

RELEASE_STATS: re.Pattern[str] = re.compile(r"(\w+) Release .*? consists of ([\d,]+) entries")

DATE_FORMAT = "%d-%b-%Y"


def parse_relnotes(relnotes_path: Path) -> dict[str, Any]:
    """Open and read the release notes file, returning it as a text string.

    :param relnotes_path: path to the release notes file
    :type relnotes_path: Path
    :return: string
    :rtype: str
    """
    rel_text = relnotes_path.read_text()
    return parse(rel_text)


def parse(relnotes: str) -> dict[str, Any]:
    """Parse the release notes for a UniProt release.

    :param relnotes: contents of the release notes file as a string
    :type relnotes: str
    :return: key-value pairs with vital release stats
    :rtype: dict[str, Any]
    """
    errors = []
    stats = {}

    relnotes_parts = relnotes.strip().split("\n\n", 1)
    if len(relnotes_parts) != 2:
        log_and_die("Could not find double line break. Relnotes file format may have changed", RuntimeError)

    (intro_str, stats_str) = relnotes_parts

    # remove line breaks for ease of parsing
    intro_str = intro_str.replace("\n", " ")

    rv = re.search(RELEASE_VERSION_DATE, intro_str)
    if not rv:
        errors.append("Could not find text matching the release version date regex.")
    else:
        stats["version"] = rv.groups()[0]
        stats["date_published"] = dt.datetime.strptime(rv.groups()[1], DATE_FORMAT)  # noqa: DTZ007

    # parse the stats section
    uniprot_trembl = re.search(UNIPROT_TREMBL_STATS, stats_str)
    if not uniprot_trembl:
        errors.append("Could not find text matching the UniProt/TrEMBL stats regex.")
    else:
        numbers = uniprot_trembl.groups()
        stats["UniProtKB"] = numbers[0]
        stats["UniProtKB/Swiss-Prot"] = numbers[1]
        stats["UniProtKB/TrEMBL"] = numbers[2]

    all_releases = re.findall(RELEASE_STATS, stats_str)
    if not all_releases:
        errors.append("Could not find text matching the release stats regex.")
    else:
        for release in all_releases:
            (db, number) = release
            if db not in stats:
                stats[db] = number

    # make sure that we have all the UniRef stats
    errors.extend([f"No stats for UniRef{n} found." for n in ["50", "90", "100"] if f"UniRef{n}" not in stats])

    if errors:
        log_and_die("\n".join(errors), RuntimeError)

    return stats
