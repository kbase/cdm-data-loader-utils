"""
UniProt XML Delta Lake Ingestion Pipeline
=========================================

This script parses UniProt XML (.xml.gz) file and ingests the data into structured Delta Lake tables.

Typical usage:
--------------
Use it in Berdle as:
python3 src/parsers/uniprot.py \
    --xml-url "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz" \
    --output-dir "./output" \
    --namespace "uniprot_db" \
    --batch-size 5000

Use it in the local computer as: 
python3 uniprot.py \
  --xml-url "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz" \
  --output-dir "./output_archaea" \
  --namespace "uniprot_archaea_db" \
  --batch-size 5000

python -m cdm_data_loader_utils.parsers.uniprot \
  --xml-url "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/uniprot_sprot_archaea.xml.gz" \
  --output-dir "./output" \
  --namespace "uniprot_db" \
  --batch-size 5000 \


Arguments:
----------
--xml-url:      URL to the UniProt XML .gz file 
--output-dir:   Output directory for Delta tables and logs (default: './output')
--namespace:    Delta Lake database name (default: 'uniprot_db')
--target-date:  Process entries modified/updated since specific date 
--batch-size:   Number of UniProt entries to process per write batch (default: 5000)

Functionality:
--------------
- Downloads the XML file if not present locally 
- Parses UniProt entries in a memory-efficient streaming fashion
- Maps parsed data into standardized CDM tables
- Writes all tables as Delta Lake tables, supporting incremental import
- Supports overwrite of previous imports and incremental updates by unique entity_id

Typical scenario:
-----------------
- Large-scale UniProt batch parsing and warehouse ingestion
- Efficient data lake ingestion

"""

import datetime
import gzip
import json
import logging
import os
import uuid
import xml.etree.ElementTree as ET

import click
import requests
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.functions import col, split
from pyspark.sql.types import ArrayType, StringType, StructField, StructType

from cdm_data_loader_utils.parsers.shared_identifiers import parse_identifiers_generic
from cdm_data_loader_utils.parsers.xml_utils import clean_dict, find_all_text, get_attr, get_text, parse_db_references

# ---------------------------------------------------------------------
#                              Logging
# ---------------------------------------------------------------------
logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
)


# ---------------------------------------------------------------------
# XML namespace mapping for UniProt entries (used for all XPath queries)
# ---------------------------------------------------------------------
NS = {"ns": "https://uniprot.org/uniprot"}


# ---------------------------------------------------------------------
# Stable ID namespace (UUIDv5)
# ---------------------------------------------------------------------
CDM_UUID_NAMESPACE = uuid.UUID("2d3f6e2a-4d7b-4a8c-9c5a-0e0f7b7d9b3a")


# ---------------------------------------------------------------------
# CURIE prefixes
# ---------------------------------------------------------------------
PREFIX_TRANSLATION: dict[str, str] = {
    "UniProtKB": "UniProt",
    "UniProtKB/Swiss-Prot": "UniProt",
    "UniProtKB/TrEMBL": "UniProt",
    "UniParc": "UniParc",
    "RefSeq": "RefSeq",
    "EMBL": "EMBL",
    "PDB": "PDB",
    "ChEBI": "ChEBI",
    "Rhea": "Rhea",
    "NCBI Taxonomy": "NCBITaxon",
    "GeneID": "NCBIGene",
    "Ensembl": "Ensembl",
    "GO": "GO",
}


# ================================ HELPERS =================================
def delta_table_path(output_dir: str, namespace: str, table: str) -> str:
    return os.path.abspath(os.path.join(output_dir, namespace, table))


def build_datasource_record(xml_url: str) -> dict:
    """Build a provenance record for the UniProt datasource."""
    return {
        "name": "UniProt import",
        "source": "UniProt",
        "url": xml_url,
        "accessed": datetime.datetime.now(datetime.UTC).isoformat(),
        "version": 115,
    }


def save_datasource_record(xml_url: str, output_dir: str) -> dict:
    """Generate and save the datasource provenance record as a JSON file."""
    datasource = build_datasource_record(xml_url)

    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "datasource.json")

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(datasource, f, indent=2)

    logger.info("Saved datasource record to %s", output_path)
    return datasource


def download_file(
    url: str,
    output_path: str,
    chunk_size: int = 1024 * 1024,
    overwrite: bool = False,
) -> None:
    """Download URL -> output_path (streaming)"""
    if os.path.exists(output_path) and not overwrite:
        logger.info("File already exists, skip download: %s", output_path)
        return

    tmp_path = output_path + ".part"
    if os.path.exists(tmp_path):
        try:
            os.remove(tmp_path)
        except Exception:
            pass

    try:
        logger.info("Downloading %s -> %s", url, output_path)
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(tmp_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
        os.replace(tmp_path, output_path)
        logger.info("Download complete: %s", output_path)
    except Exception:
        logger.exception("Failed to download %s", url)
        try:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        except Exception:
            logger.exception("Failed to remove partial download: %s", tmp_path)
        raise


def prepare_local_xml(xml_url: str, output_dir: str, overwrite: bool = False) -> str:
    os.makedirs(output_dir, exist_ok=True)
    local_path = os.path.join(output_dir, os.path.basename(xml_url))
    download_file(xml_url, local_path, overwrite=overwrite)
    return local_path


def stream_uniprot_xml(filepath: str):
    """Stream gzipped UniProt XML entries."""
    logger.info("Streaming UniProt XML from: %s", filepath)
    with gzip.open(filepath, "rb") as f:
        for _, elem in ET.iterparse(f, events=("end",)):
            if elem.tag.endswith("entry"):
                yield elem
                elem.clear()


def get_spark_session(namespace: str) -> SparkSession:
    """Initialize SparkSession with Delta Lake support, and ensure the target database exists."""
    builder = (
        SparkSession.builder.appName("UniProtDeltaIngestion")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config(
            "spark.sql.catalog.spark_catalog",
            "org.apache.spark.sql.delta.catalog.DeltaCatalog",
        )
        .config("spark.databricks.delta.schema.autoMerge.enabled", "true")
    )
    spark = configure_spark_with_delta_pip(builder).getOrCreate()
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")
    return spark


def normalize_prefix(db_type: str) -> str:
    """Map UniProt dbReference @type to a normalized CURIE prefix."""
    return PREFIX_TRANSLATION.get(db_type, db_type.replace(" ", ""))


def make_curie(db_type: str, db_id: str) -> str:
    """Create CURIE with normalized prefix."""
    return f"{normalize_prefix(db_type)}:{db_id}"


# ================================ STABLE ID =================================
def stable_cdm_id_from_uniprot_accession(accession: str, prefix: str = "cdm_prot_") -> str:
    u = uuid.uuid5(CDM_UUID_NAMESPACE, f"UniProt:{accession}")
    return f"{prefix}{u}"


def load_existing_maps(
    spark: SparkSession,
    output_dir: str,
    namespace: str,
) -> tuple[dict[str, str], dict[str, str]]:
    """
    Returns:
      accession_to_entity_id: accession -> entity_id  (from identifiers)
      entity_id_to_created:   entity_id -> created    (from entities)
    """
    accession_to_entity_id: dict[str, str] = {}
    entity_id_to_created: dict[str, str] = {}

    id_path = os.path.join(output_dir, namespace, "identifiers")
    ent_path = os.path.join(output_dir, namespace, "entities")

    if os.path.exists(id_path):
        try:
            df = (
                spark.read.format("delta")
                .load(id_path)
                .filter(col("identifier").startswith("UniProt:"))
                .select(
                    split(col("identifier"), ":").getItem(1).alias("accession"),
                    col("entity_id"),
                )
            )
            for row in df.toLocalIterator():
                acc = row["accession"]
                eid = row["entity_id"]
                if acc and eid:
                    accession_to_entity_id[acc] = eid
            logger.info(
                "Loaded %d accession->entity_id from %s",
                len(accession_to_entity_id),
                id_path,
            )
        except Exception:
            logger.exception("Couldn't load identifiers from %s", id_path)

    if os.path.exists(ent_path):
        try:
            df = spark.read.format("delta").load(ent_path).select("entity_id", "created")
            for row in df.toLocalIterator():
                if row["entity_id"] and row["created"]:
                    entity_id_to_created[row["entity_id"]] = row["created"]
            logger.info(
                "Loaded %d entity_id->created from %s",
                len(entity_id_to_created),
                ent_path,
            )
        except Exception:
            logger.exception("Couldn't load entities from %s", ent_path)

    return accession_to_entity_id, entity_id_to_created


# ================================ PARSERS =================================
def parse_identifiers(entry, cdm_id: str) -> list[dict]:
    out = parse_identifiers_generic(entry=entry, xpath="ns:accession", prefix="UniProt", ns=NS)
    for row in out:
        row["entity_id"] = cdm_id
        row.setdefault("source", "UniProt")
        row.setdefault("description", "UniProt accession")
    return out


def _make_name_record(cdm_id: str, name_text: str, description: str) -> dict:
    return {
        "entity_id": cdm_id,
        "name": name_text,
        "description": description,
        "source": "UniProt",
    }


def parse_names(entry, cdm_id: str) -> list[dict]:
    names: list[dict] = []

    for txt in find_all_text(entry, "ns:name", NS):
        names.append(_make_name_record(cdm_id, txt, "UniProt entry name"))

    protein = entry.find("ns:protein", NS)
    if protein is not None:
        for tag_name, logical_type in [
            ("recommendedName", "recommended"),
            ("alternativeName", "alternative"),
        ]:
            for name_block in protein.findall(f"ns:{tag_name}", NS):
                for xml_tag, length_label in [
                    ("fullName", "full"),
                    ("shortName", "short"),
                ]:
                    elem = name_block.find(f"ns:{xml_tag}", NS)
                    text = get_text(elem)
                    if text:
                        names.append(
                            _make_name_record(
                                cdm_id,
                                text,
                                f"UniProt {logical_type} {length_label} name",
                            )
                        )
    return names


def parse_protein_info(entry, cdm_id: str) -> dict | None:
    protein_info: dict = {}

    protein = entry.find("ns:protein", NS)
    if protein is not None:
        ec_paths = ["ns:recommendedName/ns:ecNumber", "ns:alternativeName/ns:ecNumber"]
        ec_numbers: list[str] = []
        for path in ec_paths:
            ec_numbers.extend(find_all_text(protein, path, NS))
        if ec_numbers:
            protein_info["ec_numbers"] = ";".join(ec_numbers)

    protein_existence = entry.find("ns:proteinExistence", NS)
    if protein_existence is not None:
        protein_info["protein_id"] = cdm_id
        protein_info["evidence_for_existence"] = get_attr(protein_existence, "type")

    seq_elem = entry.find("ns:sequence", NS)
    if seq_elem is not None:
        protein_info.update(
            clean_dict({
                "length": get_attr(seq_elem, "length"),
                "mass": get_attr(seq_elem, "mass"),
                "checksum": get_attr(seq_elem, "checksum"),
                "modified": get_attr(seq_elem, "modified"),
                "sequence_version": get_attr(seq_elem, "version"),
                "sequence": get_text(seq_elem),
            })
        )

    entry_modified = get_attr(entry, "modified") or get_attr(entry, "updated")
    if entry_modified:
        protein_info["entry_modified"] = entry_modified

    return protein_info if protein_info else None


def parse_evidence_map(entry) -> dict[str, dict]:
    evidence_map: dict[str, dict] = {}

    for ev in entry.findall("ns:evidence", NS):
        key = get_attr(ev, "key")
        if not key:
            continue

        evidence_type = get_attr(ev, "type")
        pubs: list[str] = []
        others: list[str] = []

        source = ev.find("ns:source", NS)
        if source is not None:
            raw_pubs, raw_others = parse_db_references(source, NS)

            normalized_pubs: list[str] = []
            for p in raw_pubs:
                up = p.upper()
                if up.startswith("PUBMED:"):
                    _, acc = p.split(":", 1)
                    normalized_pubs.append(f"PMID:{acc}")
                else:
                    normalized_pubs.append(p)

            pubs = normalized_pubs
            others = raw_others

        evidence_map[key] = clean_dict({
            "evidence_type": evidence_type,
            "publications": pubs or None,
            "supporting_objects": others or None,
        })

    return evidence_map


def _make_association(
    cdm_id: str,
    obj: str,
    predicate: str | None = None,
    evidence_key: str | None = None,
    evidence_map: dict | None = None,
) -> dict:
    assoc = {
        "subject": cdm_id,
        "object": obj,
        "predicate": predicate,
        "evidence_type": None,
        "supporting_objects": None,
        "publications": None,
    }
    if evidence_key and evidence_map and evidence_key in evidence_map:
        assoc.update(evidence_map[evidence_key])
    return clean_dict(assoc)


def parse_reaction_association(reaction, cdm_id: str, evidence_map: dict[str, dict]) -> list[dict]:
    associations: list[dict] = []
    for dbref in reaction.findall("ns:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        if not db_type or not db_id:
            continue

        assoc = {
            "subject": cdm_id,
            "predicate": "catalyzes",
            "object": make_curie(db_type, db_id),
            "evidence_type": None,
            "supporting_objects": None,
            "publications": None,
        }
        evidence_key = reaction.get("evidence")
        if evidence_key and evidence_key in evidence_map:
            assoc.update(evidence_map[evidence_key])
        associations.append(clean_dict(assoc))
    return associations


def parse_cofactor_association(cofactor, cdm_id: str) -> list[dict]:
    associations: list[dict] = []
    for dbref in cofactor.findall("ns:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        if not db_type or not db_id:
            continue
        associations.append(
            clean_dict({
                "subject": cdm_id,
                "predicate": "requires_cofactor",
                "object": make_curie(db_type, db_id),
                "evidence_type": None,
                "supporting_objects": None,
                "publications": None,
            })
        )
    return associations


def parse_associations(entry, cdm_id: str, evidence_map: dict[str, dict]) -> list[dict]:
    """
    Only keep:
      - taxonomy association
      - catalytic activity / cofactor associations
    """
    associations: list[dict] = []

    # Taxonomy association
    organism = entry.find("ns:organism", NS)
    if organism is not None:
        taxon_ref = organism.find('ns:dbReference[@type="NCBI Taxonomy"]', NS)
        if taxon_ref is not None:
            tax_id = taxon_ref.get("id")
            if tax_id:
                associations.append(_make_association(cdm_id, f"NCBITaxon:{tax_id}", predicate="in_taxon"))

    # Catalytic activity / cofactor
    for comment in entry.findall("ns:comment", NS):
        comment_type = comment.get("type")
        if comment_type == "catalytic activity":
            for reaction in comment.findall("ns:reaction", NS):
                associations.extend(parse_reaction_association(reaction, cdm_id, evidence_map))
        elif comment_type == "cofactor":
            for cofactor in comment.findall("ns:cofactor", NS):
                associations.extend(parse_cofactor_association(cofactor, cdm_id))

    return associations


def parse_cross_references(entry, cdm_id: str) -> list[dict]:
    """Generic <dbReference> -> cross_references table."""
    rows: list[dict] = []

    for dbref in entry.findall("ns:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        if not db_type or not db_id:
            continue

        xref_type = normalize_prefix(db_type)

        if ":" in db_id:
            xref = db_id
        else:
            xref = f"{xref_type}:{db_id}"

        rows.append(
            clean_dict({
                "entity_id": cdm_id,
                "xref_type": xref_type,
                "xref_value": db_id,
                "xref": xref,
            })
        )

    return rows


def parse_publications(entry) -> list[str]:
    publications: list[str] = []
    for reference in entry.findall("ns:reference", NS):
        citation = reference.find("ns:citation", NS)
        if citation is None:
            continue

        raw_pubs, _ = parse_db_references(citation, NS)
        for p in raw_pubs:
            up = p.upper()
            if up.startswith("PUBMED:"):
                _, acc = p.split(":", 1)
                publications.append(f"PMID:{acc}")
            elif up.startswith("DOI:"):
                _, acc = p.split(":", 1)
                publications.append(f"DOI:{acc}")

    return list(dict.fromkeys(publications))


def parse_uniprot_entry(
    entry,
    cdm_id: str,
    current_timestamp: str,
    datasource_name: str = "UniProt import",
    prev_created: str | None = None,
) -> dict:
    entity_created = prev_created or current_timestamp
    entity_updated = current_timestamp

    uniprot_created = entry.attrib.get("created")
    uniprot_modified = entry.attrib.get("modified") or entry.attrib.get("updated")
    uniprot_version = entry.attrib.get("version")

    entity = {
        "entity_id": cdm_id,
        "entity_type": "protein",
        "data_source": datasource_name,
        "created": entity_created,
        "updated": entity_updated,
        "version": uniprot_version,
        "uniprot_created": uniprot_created,
        "uniprot_modified": uniprot_modified,
    }

    evidence_map = parse_evidence_map(entry)

    return {
        "entity": entity,
        "identifiers": parse_identifiers(entry, cdm_id),
        "names": parse_names(entry, cdm_id),
        "protein": parse_protein_info(entry, cdm_id),
        "associations": parse_associations(entry, cdm_id, evidence_map),
        "cross_references": parse_cross_references(entry, cdm_id),
        "publications": parse_publications(entry),
    }


# ================================ SCHEMA =================================
schema_entities = StructType([
    StructField("entity_id", StringType(), False),
    StructField("entity_type", StringType(), False),
    StructField("data_source", StringType(), False),
    StructField("created", StringType(), True),
    StructField("updated", StringType(), True),
    StructField("version", StringType(), True),
    StructField("uniprot_created", StringType(), True),
    StructField("uniprot_modified", StringType(), True),
])

schema_identifiers = StructType([
    StructField("entity_id", StringType(), False),
    StructField("identifier", StringType(), False),
    StructField("source", StringType(), True),
    StructField("description", StringType(), True),
])

schema_proteins = StructType([
    StructField("protein_id", StringType(), False),
    StructField("ec_numbers", StringType(), True),
    StructField("evidence_for_existence", StringType(), True),
    StructField("length", StringType(), True),
    StructField("mass", StringType(), True),
    StructField("checksum", StringType(), True),
    StructField("modified", StringType(), True),
    StructField("sequence_version", StringType(), True),
    StructField("sequence", StringType(), True),
    StructField("entry_modified", StringType(), True),
])

schema_names = StructType([
    StructField("entity_id", StringType(), False),
    StructField("name", StringType(), False),
    StructField("description", StringType(), True),
    StructField("source", StringType(), True),
])

schema_associations = StructType([
    StructField("subject", StringType(), True),
    StructField("object", StringType(), True),
    StructField("predicate", StringType(), True),
    StructField("evidence_type", StringType(), True),
    StructField("supporting_objects", ArrayType(StringType()), True),
    StructField("publications", ArrayType(StringType()), True),
])

schema_cross_references = StructType([
    StructField("entity_id", StringType(), False),
    StructField("xref_type", StringType(), True),
    StructField("xref_value", StringType(), True),
    StructField("xref", StringType(), True),
])

schema_publications = StructType([
    StructField("entity_id", StringType(), False),
    StructField("publication", StringType(), True),
])


# ================================ DELTA WRITE =================================
def ensure_tables_registered(spark: SparkSession, output_dir: str, namespace: str, table_names: list[str]) -> None:
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")
    for tbl in table_names:
        # delta_dir = os.path.abspath(os.path.join(output_dir, namespace, tbl))
        delta_dir = delta_table_path(output_dir, namespace, tbl)
        spark.sql(
            f"""
            CREATE TABLE IF NOT EXISTS {namespace}.{tbl}
            USING DELTA
            LOCATION '{delta_dir}'
            """
        )


def save_batches_to_delta(
    spark: SparkSession,
    tables: dict[str, tuple[list, StructType]],
    output_dir: str,
    namespace: str,
    mode: str = "append",
) -> None:
    for table_name, (records, schema) in tables.items():
        if not records:
            continue

        # delta_dir = os.path.abspath(os.path.join(output_dir, namespace, table_name))
        delta_dir = delta_table_path(output_dir, namespace, table_name)
        df = spark.createDataFrame(records, schema)
        writer = df.write.format("delta").mode(mode)

        if mode == "append":
            writer = writer.option("mergeSchema", "true")
        if mode == "overwrite":
            writer = writer.option("overwriteSchema", "true")

        writer.save(delta_dir)


## =============================== MAIN PARSING LOOP =================================
def parse_entries(
    local_xml_path: str,
    target_date: str | None,
    batch_size: int,
    spark: SparkSession,
    tables: dict[str, tuple[list, StructType]],
    output_dir: str,
    namespace: str,
    current_timestamp: str,
    accession_to_entity_id: dict[str, str],
    entity_id_to_created: dict[str, str],
    mode: str,
) -> tuple[int, int]:
    target_date_dt = None
    if target_date:
        try:
            target_date_dt = datetime.datetime.strptime(target_date, "%Y-%m-%d")
            logger.info("Target date filter enabled: >= %s", target_date)
        except Exception:
            logger.warning("Invalid target date provided: %s (ignored)", target_date)
            target_date_dt = None

    entry_count, skipped = 0, 0

    for entry_elem in stream_uniprot_xml(local_xml_path):
        try:
            mod_date = entry_elem.attrib.get("modified") or entry_elem.attrib.get("updated")

            if target_date_dt and mod_date:
                try:
                    entry_date_dt = datetime.datetime.strptime(mod_date[:10], "%Y-%m-%d")
                    if entry_date_dt < target_date_dt:
                        skipped += 1
                        continue
                except Exception:
                    skipped += 1
                    continue

            main_accession_elem = entry_elem.find("ns:accession", NS)
            if main_accession_elem is None or not main_accession_elem.text:
                skipped += 1
                continue

            accession = main_accession_elem.text.strip()

            cdm_id = accession_to_entity_id.get(accession) or stable_cdm_id_from_uniprot_accession(accession)
            prev_created = entity_id_to_created.get(cdm_id)

            record = parse_uniprot_entry(entry_elem, cdm_id, current_timestamp, prev_created=prev_created)

            tables["entities"][0].append(record["entity"])
            tables["identifiers"][0].extend(record["identifiers"])
            tables["names"][0].extend(record["names"])

            if record["protein"]:
                tables["proteins"][0].append(record["protein"])

            tables["associations"][0].extend(record["associations"])
            tables["cross_references"][0].extend(record["cross_references"])

            for pub in record["publications"]:
                tables["publications"][0].append({
                    "entity_id": cdm_id,
                    "publication": pub,
                })

            entry_count += 1

            if entry_count % batch_size == 0:
                save_batches_to_delta(spark, tables, output_dir, namespace, mode=mode)
                for v in tables.values():
                    v[0].clear()
                logger.info("Processed and saved %d entries...", entry_count)

        except Exception:
            logger.exception("Error parsing UniProt entry, skipping")
            skipped += 1

    save_batches_to_delta(spark, tables, output_dir, namespace, mode=mode)
    return entry_count, skipped


def ingest_uniprot(
    xml_url: str,
    output_dir: str,
    namespace: str,
    target_date: str | None = None,
    batch_size: int = 5000,
    mode: str = "append",
    overwrite_download: bool = False,
) -> None:
    current_timestamp = datetime.datetime.now(datetime.UTC).isoformat()

    local_xml_path = prepare_local_xml(xml_url, output_dir, overwrite=overwrite_download)
    save_datasource_record(xml_url, output_dir)

    spark = get_spark_session(namespace)
    if mode == "append":
        accession_to_entity_id, entity_id_to_created = load_existing_maps(spark, output_dir, namespace)
    else:
        accession_to_entity_id, entity_id_to_created = {}, {}

    # accession_to_entity_id, entity_id_to_created = load_existing_maps(spark, output_dir, namespace)

    entities: list[dict] = []
    identifiers: list[dict] = []
    names: list[dict] = []
    proteins: list[dict] = []
    associations: list[dict] = []
    cross_references: list[dict] = []
    publications: list[dict] = []

    tables: dict[str, tuple[list, StructType]] = {
        "entities": (entities, schema_entities),
        "identifiers": (identifiers, schema_identifiers),
        "names": (names, schema_names),
        "proteins": (proteins, schema_proteins),
        "associations": (associations, schema_associations),
        "cross_references": (cross_references, schema_cross_references),
        "publications": (publications, schema_publications),
    }

    ensure_tables_registered(
        spark,
        output_dir,
        namespace,
        [
            "entities",
            "identifiers",
            "names",
            "proteins",
            "associations",
            "cross_references",
            "publications",
        ],
    )

    logger.info(
        "Starting UniProt ingestion: xml=%s | namespace=%s | mode=%s | batch_size=%d",
        xml_url,
        namespace,
        mode,
        batch_size,
    )

    entry_count, skipped = parse_entries(
        local_xml_path=local_xml_path,
        target_date=target_date,
        batch_size=batch_size,
        spark=spark,
        tables=tables,
        output_dir=output_dir,
        namespace=namespace,
        current_timestamp=current_timestamp,
        accession_to_entity_id=accession_to_entity_id,
        entity_id_to_created=entity_id_to_created,
        mode=mode,
    )

    logger.info("Completed parsing UniProt XML. processed=%d skipped=%d", entry_count, skipped)

    logger.info("Verifying Delta tables in namespace `%s`", namespace)
    spark.sql(f"SHOW TABLES IN {namespace}").show(truncate=False)

    for tbl in [
        "entities",
        "identifiers",
        "names",
        "proteins",
        "associations",
        "cross_references",
        "publications",
    ]:
        logger.info("Verifying table: %s.%s", namespace, tbl)
        spark.sql(f"SELECT COUNT(*) AS row_count FROM {namespace}.{tbl}").show(truncate=False)
        spark.sql(f"SELECT * FROM {namespace}.{tbl} LIMIT 5").show(truncate=False)

    spark.stop()
    logger.info("Done")


# ================================ CLI =================================
@click.command()
@click.option("--xml-url", required=True, help="URL to UniProt XML (.xml.gz)")
@click.option(
    "--output-dir",
    default="output",
    show_default=True,
    help="Output directory for Delta tables",
)
@click.option(
    "--namespace",
    default="uniprot_db",
    show_default=True,
    help="Delta Lake database name",
)
@click.option(
    "--target-date",
    default=None,
    help="Only process entries modified/updated since this date (YYYY-MM-DD)",
)
@click.option(
    "--batch-size",
    default=5000,
    show_default=True,
    help="Batch size for writing Delta tables",
)
@click.option(
    "--mode",
    type=click.Choice(["append", "overwrite"]),
    default="append",
    show_default=True,
)
@click.option(
    "--overwrite-download",
    is_flag=True,
    help="Force re-download XML even if file exists",
)
def main(xml_url, output_dir, namespace, target_date, batch_size, mode, overwrite_download):
    ingest_uniprot(
        xml_url=xml_url,
        output_dir=output_dir,
        namespace=namespace,
        target_date=target_date,
        batch_size=int(batch_size),
        mode=mode,
        overwrite_download=overwrite_download,
    )


if __name__ == "__main__":
    main()
