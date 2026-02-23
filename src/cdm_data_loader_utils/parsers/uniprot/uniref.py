"""UniRef XML File Parser.

Core parser for UniRef XML entries.

Captures the following data about each UniRef cluster:
- ID and name
- members, including representative member and seed
- sequence of the representative member
- whether this is UniRef50, 90, or 100.

Data goes into the Cluster, ClusterMember and Entity tables.

"""

import datetime
from pathlib import Path
from typing import Any

from lxml.etree import Element, tounicode

from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger

UNIREF_URL = "http://uniprot.org/uniref"
UNIREF_URL_BRKT = f"{{{UNIREF_URL}}}"
NS = "ns"
UNIREF_NS = {NS: UNIREF_URL, "": UNIREF_URL}
UNIREF = "UniRef"
UNIREF_VARIANTS = [100, 90, 50]
ENTITY_ID = "entity_id"

PREFIX_TRANSLATION = {
    "UniProtKB ID": "uniprot_name",  # this is actually the uniprot name
    "UniProtKB accession": "uniprot",
    "UniParc ID": "uniparc",
    "UniRef90 ID": "uniref",
    "UniRef50 ID": "uniref",
    "UniRef100 ID": "uniref",
    "NCBI taxonomy": "NCBITaxon",
}

logger = get_cdm_logger()


def dump_xml_element(element: Element) -> str:
    """Dump an XML element as a string.

    :param element: the XML element
    :type element: Element
    :return: the element as a string
    :rtype: str
    """
    unicode_version = tounicode(element)
    return (
        unicode_version.replace(' xmlns="https://uniprot.org/uniref"', "")
        .replace(' xmlns="http://uniprot.org/uniref"', "")
        .replace(' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"', "")
    )


def generate_dbxref(db: str, acc: str) -> str:
    """Generate a database reference that uses BioRegistry prefixes."""
    return f"{PREFIX_TRANSLATION.get(db, db)}:{acc}"


def extract_cluster(entry: Element) -> tuple[dict[str, Any], dict[str, Any]]:
    """Extract the UniRef cluster ID and name."""
    entry_id = entry.get("id")
    if not entry_id:
        return ({}, {})
    cluster_id = f"uniref:{entry_id}"
    cluster = {
        "cluster_id": cluster_id,
        "name": entry.findtext("./ns:name", namespaces=UNIREF_NS),
        "cluster_type": "Protein",
        "description": None,
        # fill in protocol afterwards
    }
    entity = {
        "entity_id": cluster_id,
        "entity_type": "Cluster",
        "data_source": UNIREF,
        "data_source_entity_id": entry_id,
        "data_source_updated": entry.get("updated"),
        # fill in other columns later
    }

    return (cluster, entity)


def extract_cross_refs(
    dbref: Element,
    cluster_id: str,
    is_representative: bool,  # noqa: FBT001
) -> dict[str, Any]:
    """Extract UniRef90/50/UniParc cross references from a single <dbReference> element."""
    entity_db = dbref.get("type")
    entity_xref = dbref.get("id")
    if not entity_xref or not entity_db:
        # should never happen!
        logger.warning("Incomplete dbReference in %s", cluster_id)
        return {}

    entity_id = generate_dbxref(entity_db, entity_xref)

    # if this is a UniProtKB ID, check whether there is a UniProtKB accession
    if entity_db == "UniProtKB ID":
        uniprot_id = None
        uniprot_acc = dbref.find("ns:property[@type='UniProtKB accession']", UNIREF_NS)
        if uniprot_acc is not None:
            uniprot_id = uniprot_acc.get("value")
        if uniprot_id:
            entity_id = f"uniprot:{uniprot_id}"
        else:
            logger.warning("Could not find UniProt accession for UniProtKB ID %s", entity_xref)

    # is this a seed member?
    is_seed_prop = dbref.find("ns:property[@type='isSeed']", UNIREF_NS)
    is_seed = is_seed_prop is not None and is_seed_prop.get("value") == "true"

    return {
        "entity_id": entity_id,
        "cluster_id": cluster_id,
        "is_representative": is_representative,
        "is_seed": is_seed,
    }


def parse_uniref_entry(
    entry: Element, timestamp: datetime.datetime, uniref_variant: str, file_path: str | Path
) -> dict[str, list[dict[str, str]]]:
    """Parse a single UniRef <entry> element into CDM-friendly row tuples."""
    # Cluster basic info
    (cluster, entity) = extract_cluster(entry)
    cluster_id = cluster.get("cluster_id")
    if not cluster_id:
        logger.debug("No cluster ID found in entry")
        return {
            "_parse_error": [
                {"xml": dump_xml_element(entry), "error": "No cluster ID found in entry", "source_file": str(file_path)}
            ]
        }

    # add the current timestamp
    entity["updated"] = timestamp
    # add in the protocol (which version of UniRef this is)
    cluster["protocol"] = uniref_variant

    cluster_members = []
    # Cross references from representative and members
    repr_db = entry.find("ns:representativeMember/ns:dbReference", UNIREF_NS)
    if repr_db is not None:
        cluster_members.append(extract_cross_refs(repr_db, cluster_id, is_representative=True))

    cluster_members.extend(
        [
            extract_cross_refs(mem, cluster_id, is_representative=False)
            for mem in entry.findall("ns:member/ns:dbReference", UNIREF_NS)
        ]
    )

    # save a record of which entry is in which file
    entity_data_source = {
        "entity_id": cluster_id,
        "data_source": entity.get("data_source"),
        "source_file": str(file_path),
    }

    return {
        "entity": [entity],
        "entity_x_source_file": [entity_data_source],
        "cluster": [cluster],
        "clustermember": [c for c in cluster_members if c],
    }
