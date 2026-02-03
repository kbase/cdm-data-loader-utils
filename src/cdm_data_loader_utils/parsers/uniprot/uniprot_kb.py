"""
UniProt XML Parser Core.

Contains the core logic for parsing a UniProt XML file that can be successfully validated against the UniProt XSD.
"""

import datetime
from pathlib import Path
from typing import Any

from lxml.etree import Element, tounicode

from cdm_data_loader_utils.utils.cdm_logger import get_cdm_logger
from cdm_data_loader_utils.utils.helpers import _ensembl_type
from cdm_data_loader_utils.utils.xml_utils import get_text

CONTENT = "content"
DB = "db"
DESCRIPTION = "description"
ENTITY_ID = "entity_id"
HAS_HOST = "RO:0002454"
HOST_OF = "RO:0002453"
IN_TAXON = "RO:0002162"
KEY = "key"
NAME = "name"
NCBI_TAXON = "NCBITaxon"
XREF = "xref"

# XML namespace mapping for UniProt entries
HTTPS_NS = "https://uniprot.org/uniprot"
NS = {"ns": HTTPS_NS, "": HTTPS_NS}
APP_NAME = "uniprotkb_importer"
ENTRY_XML_TAG = f"{{{HTTPS_NS}}}entry"


# CURIE prefixes
PREFIX_TRANSLATION: dict[str, str] = {
    "UniProtKB": "UniProt",
    "UniProtKB/Swiss-Prot": "UniProt",
    "UniProtKB/TrEMBL": "UniProt",
    "NCBI Taxonomy": NCBI_TAXON,
    "GeneID": "NCBIGene",
    "pubmed": "PMID",
}

logger = get_cdm_logger()


# TODO: FIXME!!!
def build_datasource_record(xml_url: str) -> dict:
    """Build a provenance record for the UniProt datasource."""
    return {
        NAME: "UniProt import",
        "source": "UniProt",
        "url": xml_url,
        "accessed": datetime.datetime.now(datetime.UTC).isoformat(),
        "version": 115,
    }


def parse_identifiers(entry: Element) -> list[dict[str, Any]]:
    """Parse the UniProt accessions from an entry element.

    :param entry: entry XML element
    :type entry: Element
    :return: list of identifiers in "db" / "xref" format
    :rtype: list[dict[str, str]]
    """
    return [
        {DB: "UniProt", XREF: acc.text, DESCRIPTION: "UniProt accession"}
        for acc in entry.findall("./ns:accession", NS)
        if get_text(acc)
    ]


def parse_names(entry: Element) -> list[dict[str, str]]:
    """Parse the names from an entry element.

    :param entry: entry XML element
    :type entry: Element
    :return: list of names with a description of the name type
    :rtype: list[dict[str, str]]
    """
    name_list = [
        {NAME: n.text, DESCRIPTION: "UniProt entry name"} for n in entry.findall("./ns:name", NS) if get_text(n)
    ]
    for name_type in ["recommended", "alternative", "submitted"]:
        for name_length in ["full", "short"]:
            name_list.extend(
                [
                    {
                        NAME: n.text,
                        DESCRIPTION: f"UniProt {name_type} {name_length} name",
                    }
                    for n in entry.findall(f"ns:protein/ns:{name_type}Name/ns:{name_length}Name", NS)
                    if get_text(n)
                ]
            )

    # gene names
    for n in entry.findall("./ns:gene/ns:name", NS):
        name_list.extend([{NAME: n.text, DESCRIPTION: f"UniProt gene name, {n.get('type')}"}])

    return name_list


def parse_cross_references(entry: Element) -> list[dict[str, Any]]:
    """Parse the dbReferences from an entry element.

    schema for dbReferences:

    simple example:
    <dbReference type="<db>" id="<local identifier>" />

    dbref with extra info:
    <dbReference type="<db>" id="<local identifier>">
        <molecule id="<uniprot ID>"/>  # 0+
        <property type="..." value="..."/>  # 0+
    </dbReference>

    :param entry: entry XML element
    :type entry: Element
    :return: list of cross references in "db" / "xref" / "description" / "relationship" format
    :rtype: list[dict[str, str]]
    """
    refs = []
    entry_acc = None
    for dbxref in entry.findall("ns:dbReference", NS):
        db: str | None = dbxref.get("type")
        xref: str | None = dbxref.get("id")
        if not db or not xref:
            continue
        molecules = dbxref.findall("ns:molecule[@id]", NS)
        if len(molecules) > 1:
            if not entry_acc:
                entry_acc = entry.find("ns:accession", NS).text
            logger.warning("Found many molecules for entry %s", entry_acc)
        description_suffix = f" for UniProt:{molecules[0].get('id')}" if len(molecules) > 0 else ""

        prop_list = [{"type": p.get("type"), "value": p.get("value")} for p in dbxref.findall("ns:property", NS)]
        if prop_list and db in ("Ensembl", "EMBL", "RefSeq", "MANE-Select"):
            if db == "Ensembl":
                refs.extend(parse_ensembl_dbxref(xref, prop_list, description_suffix))
                continue

            if db == "EMBL":
                # embl / genbank
                refs.extend(parse_embl_dbxref(xref, prop_list, description_suffix))
                continue

            if db == "RefSeq":
                refs.extend(parse_refseq_dbxref(xref, prop_list, description_suffix))
                continue

            if db == "MANE-Select":
                refs.extend(parse_mane_dbxref(xref, prop_list, description_suffix))
                continue

        if db == "GO":
            refs.extend([{DB: "GO", XREF: xref.removeprefix("GO:")}])
            continue

        if description_suffix:
            refs.extend([{DB: db, XREF: xref, DESCRIPTION: f"{db} ID{description_suffix}"}])
            continue

        # default parse
        refs.extend([{DB: db, XREF: xref}])
    return refs


def parse_ensembl_dbxref(xref: str, prop_list: list[dict[str, str]], description_suffix: str) -> list[dict[str, str]]:
    """Parse an Ensembl dbxref."""
    likely_type = _ensembl_type(xref)
    return [
        {DB: "ensembl", XREF: xref, DESCRIPTION: f"Ensembl {likely_type} ID{description_suffix}"},
        *[
            {DB: "ensembl", XREF: p["value"], DESCRIPTION: f"Ensembl {p['type']}{description_suffix}"}
            for p in prop_list
        ],
    ]


def parse_refseq_dbxref(xref: str, prop_list: list[dict[str, str]], description_suffix: str) -> list[dict[str, str]]:
    """Parse a RefSeq dbxref."""
    return [
        {DB: "refseq", XREF: xref, DESCRIPTION: f"RefSeq protein sequence ID{description_suffix}"},
        # add in refs for each of the properties
        *[{DB: "refseq", XREF: p["value"], DESCRIPTION: f"RefSeq {p['type']}{description_suffix}"} for p in prop_list],
    ]


def parse_embl_dbxref(xref: str, prop_list: list[dict[str, str]], description_suffix: str) -> list[dict[str, str]]:
    """Parse an EMBL / GenBan dbxref."""
    mol_types = [p["value"] for p in prop_list if p["type"] == "molecule type"]
    description = f"EMBL/GenBank {mol_types[0]} ID" if mol_types and len(mol_types) == 1 else "EMBL/GenBank ID"
    return [
        {DB: "genbank", XREF: xref, DESCRIPTION: f"{description}{description_suffix}"},
        *[
            {DB: "genbank", XREF: p["value"], DESCRIPTION: f"EMBL/GenBank {p['type']}{description_suffix}"}
            for p in prop_list
            if p["type"] not in ("molecule type", "status")
        ],
    ]


def parse_mane_dbxref(xref: str, prop_list: list[dict[str, str]], description_suffix: str) -> list[dict[str, str]]:
    """Parse an MANE-select dbxref."""
    likely_type = _ensembl_type(xref)
    refs = [
        {
            DB: "ensembl",
            XREF: xref,
            DESCRIPTION: f"Ensembl {likely_type} ID{description_suffix}, via MANE-Select",
        }
    ]
    for p in prop_list:
        if p["type"].startswith("RefSeq "):
            refs.append(
                {
                    DB: "refseq",
                    XREF: p["value"],
                    DESCRIPTION: f"{p['type']}{description_suffix}, via MANE-Select",
                }
            )
        elif p["value"].startswith("ENS"):
            refs.append(
                {
                    DB: "ensembl",
                    XREF: p["value"],
                    DESCRIPTION: f"Ensembl {p['type']}{description_suffix}, via MANE-Select",
                }
            )
    return refs


def parse_protein_info(entry: Element) -> dict[str, Any]:
    """Extract information for the protein table."""
    protein: dict = {}
    # <proteinExistence type="inferred from homology"/>
    # valid values:
    # <xs:enumeration value="evidence at protein level"/>
    # <xs:enumeration value="evidence at transcript level"/>
    # <xs:enumeration value="inferred from homology"/>
    # <xs:enumeration value="predicted"/>
    # <xs:enumeration value="uncertain"/>
    protein_existence = entry.find("ns:proteinExistence", NS)
    if protein_existence is not None:
        protein["evidence_for_existence"] = protein_existence.get("type")

    # <sequence length="122" mass="14969" checksum="C5E63C34B941711C" modified="2009-05-05" version="1">AMINO_ACID_SEQUENCE</sequence>
    seq_elem = entry.find("ns:sequence", NS)
    if seq_elem is not None:
        protein.update(
            {
                "length": int(seq_elem.get("length") or 0),
                "hash": seq_elem.get("checksum"),
            }
        )
    sequence = entry.findtext("sequence", namespaces=NS)
    if sequence:
        protein["sequence"] = sequence
    return protein


def parse_organism(entry: Element) -> list[dict[str, str]]:
    """Parse the organism and organism host statements."""
    taxon_ref = entry.find("ns:organism/ns:dbReference[@type='NCBI Taxonomy']", NS)
    if taxon_ref is not None and taxon_ref.get("id"):
        return [
            {
                DB: NCBI_TAXON,
                XREF: taxon_ref.get("id"),
                DESCRIPTION: "UniProt taxon designation",
                "relationship": "RO:0002162: in taxon",
            }
        ]
    return []


def parse_references(entry: Element) -> dict[str, list[dict[str, str]]]:
    """Parse the publications for an entry and extract any that have DOIs or PMIDs.

    :param entry: the UniProt entry being parsed
    :type entry: Element
    :return: list of publications associated with the entry
    :rtype: list[dict[str, str]]
    """
    publications: set[str] = set()
    all_publications = []
    for reference in entry.findall("ns:reference", NS):
        # save the full version for later parsing if needed
        all_publications.extend([{KEY: reference.get("key"), CONTENT: dump_xml_element(reference)}])
        citation_refs = [
            {DB: ref.get("type").lower(), XREF: ref.get("id")} for ref in reference.findall("citation/dbReference", NS)
        ]
        if not citation_refs:
            continue

        all_ref_types = {c[DB] for c in citation_refs}
        ref_db_priority_order = ["doi", "pmcid", "pmid", "pubmed", "agricola"]
        unexpected_types = all_ref_types - set(ref_db_priority_order)
        if unexpected_types:
            # oh no!!
            logger.warning("Unexpected dbxref types in publications: %s", ", ".join(sorted(unexpected_types)))

        ref_found = False
        for ref_type in ref_db_priority_order:
            refs = [ref[XREF] for ref in citation_refs if ref[DB] == ref_type]
            if refs:
                publications.add(f"{PREFIX_TRANSLATION.get(ref_type, ref_type.upper())}:{refs[0]}")
                ref_found = True
                break

        # oh no!!
        if not ref_found:
            logger.warning("Could not find priority ref type for reference: %s", dump_xml_element(reference))

    return {"all_xml": all_publications, "publication": [{"publication_id": pub} for pub in publications]}


def dump_evidence_xml(entry: Element, uniprot_id: str) -> list[dict[str, str]]:
    """Dump evidence as XML strings for potential later parsing.

    :param entry: entry XML element
    :type entry: Element
    :param uniprot_id: uniprot ID for the entry
    :type uniprot_id: str
    :return: list of dictionaries containing evidence XML
    :rtype: list[dict[str, str]]
    """
    return [
        {ENTITY_ID: uniprot_id, KEY: evidence.get("key"), CONTENT: dump_xml_element(evidence)}
        for evidence in entry.findall("ns:evidence", NS)
    ]


def dump_comment_xml(entry: Element, uniprot_id: str) -> list[dict[str, str]]:
    """Dump comments as XML strings for potential later parsing.

    :param entry: entry XML element
    :type entry: Element
    :param uniprot_id: uniprot ID for the entry
    :type uniprot_id: str
    :return: list of dictionaries containing comment XML
    :rtype: list[dict[str, str]]
    """
    return [{ENTITY_ID: uniprot_id, CONTENT: dump_xml_element(comment)} for comment in entry.findall("ns:comment", NS)]


def dump_xml_element(element: Element) -> str:
    """Dump an XML element as a string.

    :param element: the XML element
    :type element: Element
    :return: the element as a string
    :rtype: str
    """
    unicode_version = tounicode(element)
    return unicode_version.replace(' xmlns="https://uniprot.org/uniprot"', "").replace(
        ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"', ""
    )


def parse_uniprot_entry(
    entry: Element,
    current_timestamp: datetime.datetime,
    file_path: str | Path,
) -> dict[str, Any]:
    """Parse a UniProt entry to yield the appropriate CDM schema tables.

    :param entry: an entry from a UniProt XML file
    :type entry: Element
    :param current_timestamp: current timestamp
    :type current_timestamp: datetime.datetime
    :param file_path: path of the file currently being parsed
    :type file_path: str | Path
    :return: dictionary of parsed data, indexed by table name
    :rtype: dict[str, Any]
    """
    try:
        first_acc = entry.findtext("ns:accession", namespaces=NS)

        if not first_acc or not first_acc.strip():
            logger.error("No accession found for entry")
            return {
                "_parse_error": [
                    {"xml": dump_xml_element(entry), "error": "No accession found", "source_file": str(file_path)}
                ]
            }

        uniprot_id = f"uniprot:{first_acc.strip()}"
        entity = {
            ENTITY_ID: uniprot_id,
            "entity_type": "protein",
            "data_source_entity_id": first_acc.strip(),
            "data_source_created": entry.attrib.get("created"),
            "data_source_modified": entry.attrib.get("modified"),
            "data_source_entity_version": entry.attrib.get("version"),
            "data_source_id": None,  # SwissProt / TrEMBL
            # TODO: add in "created" field when merging datasets?
            "created": None,
            "data_source": f"UniProt/{entry.attrib.get('dataset')}",
            "updated": current_timestamp,
        }

        ref_data = parse_references(entry)
        reference_xml = [{ENTITY_ID: uniprot_id, **ref} for ref in ref_data["all_xml"]]

        # collect raw XML for entities that we will parse later
        comment_xml = dump_comment_xml(entry, uniprot_id)
        evidence_xml = dump_evidence_xml(entry, uniprot_id)

        return {
            "entity": [entity],
            "identifier": [
                {ENTITY_ID: uniprot_id, **e}
                for e in [*parse_identifiers(entry), *parse_cross_references(entry), *parse_organism(entry)]
            ],
            "name": [{ENTITY_ID: uniprot_id, **e} for e in parse_names(entry)],
            "protein": [{"protein_id": uniprot_id, **parse_protein_info(entry)}],
            "entity_x_publication": [{ENTITY_ID: uniprot_id, **e} for e in ref_data["publication"]],
            "entity_x_source_file": [
                {ENTITY_ID: uniprot_id, "data_source": entity["data_source"], "source_file": str(file_path)}
            ],
            # non-schema tables
            "_evidence_xml": evidence_xml,
            "_comment_xml": comment_xml,
            "_reference_xml": reference_xml,
        }
    except Exception as e:
        logger.exception("Parse error while processing entry")
        return {"_parse_error": [{"xml": dump_xml_element(entry), "error": str(e), "source_file": str(file_path)}]}
