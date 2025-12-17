"""UniProt XML Delta Lake Ingestion Pipeline.

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
from typing import Dict, List, Optional, Tuple

import click
import requests
from delta import configure_spark_with_delta_pip
from pyspark.sql import SparkSession
from pyspark.sql.types import (
    ArrayType,
    StringType,
    StructField,
    StructType,
)

from cdm_data_loader_utils.parsers.shared_identifiers import parse_identifiers_generic
from cdm_data_loader_utils.parsers.xml_utils import (
    clean_dict,
    find_all_text,
    get_attr,
    get_text,
    parse_db_references,
)

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


# ================================ HELPERS =================================
def build_datasource_record(xml_url: str) -> dict:
    """Build a provenance record for the UniProt datasource."""
    return {
        "name": "UniProt import",
        "source": "UniProt",
        "url": xml_url,
        "accessed": datetime.datetime.now(datetime.timezone.utc).isoformat(),
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


def download_file(url: str, output_path: str, chunk_size: int = 8192, overwrite: bool = False) -> None:
    """Download a file from a given URL to a local output path."""
    if os.path.exists(output_path) and not overwrite:
        logger.info("File already exists, skip download: %s", output_path)
        return

    try:
        logger.info("Downloading %s -> %s", url, output_path)
        with requests.get(url, stream=True, timeout=60) as response:
            response.raise_for_status()
            with open(output_path, "wb") as f:
                for chunk in response.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)
        logger.info("Download complete: %s", output_path)
    except Exception:
        logger.exception("Failed to download %s", url)

        # remove incomplete file
        try:
            if os.path.exists(output_path):
                os.remove(output_path)
        except Exception:
            logger.exception("Failed to remove incomplete download: %s", output_path)
        raise


def prepare_local_xml(xml_url: str, output_dir: str) -> str:
    """Download UniProt XML (.xml.gz) to local output directory if needed."""
    os.makedirs(output_dir, exist_ok=True)
    local_xml_path = os.path.join(output_dir, os.path.basename(xml_url))
    download_file(xml_url, local_xml_path)
    return local_xml_path


def stream_uniprot_xml(filepath: str):
    """
    Stream and parse UniProt XML entries from a local gzipped file.
    Yields each <entry> element as soon as it is parsed to avoid loading entire XML.
    """
    logger.info("Streaming UniProt XML from: %s", filepath)
    with gzip.open(filepath, "rb") as f:
        context = ET.iterparse(f, events=("end",))
        for _, element in context:
            if element.tag.endswith("entry"):
                yield element
                element.clear()


def get_spark_session(namespace: str) -> SparkSession:
    """Initialize SparkSession with Delta Lake support, and ensure the target database exists."""
    builder = (
        SparkSession.builder.appName("UniProtDeltaIngestion")
        .config("spark.sql.extensions", "io.delta.sql.DeltaSparkSessionExtension")
        .config(
            "spark.sql.catalog.spark_catalog",
            "org.apache.spark.sql.delta.catalog.DeltaCatalog",
        )
    )
    spark = configure_spark_with_delta_pip(builder).getOrCreate()
    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")
    return spark


# ================================ STABLE ID =================================
def stable_cdm_id_from_accession(accession: str, prefix: str = "cdm_prot_") -> str:
    """
    Deterministic/stable CDM ID using UUIDv5 on accession.
    This prevents entity_id duplication across runs.
    """
    u = uuid.uuid5(CDM_UUID_NAMESPACE, accession)
    return f"{prefix}{u}"


def load_existing_identifiers(spark: SparkSession, output_dir: str, namespace: str) -> Dict[str, str]:
    """
    Load existing identifiers table to map UniProt accession -> existing entity_id.
    Path: <output_dir>/<namespace>/identifiers
    """
    access_to_cdm_id: Dict[str, str] = {}
    id_path = os.path.join(output_dir, namespace, "identifiers")

    if not os.path.exists(id_path):
        logger.info("No previous identifiers delta at %s", id_path)
        return access_to_cdm_id

    try:
        df = spark.read.format("delta").load(id_path).select("identifier", "entity_id")
        for row in df.collect():
            ident = row["identifier"]
            if not ident or ":" not in ident:
                continue
            prefix, acc = ident.split(":", 1)
            if prefix.upper() == "UNIPROT" and acc:
                access_to_cdm_id[acc] = row["entity_id"]
        logger.info(
            "Loaded %d existing accession->entity_id mappings from %s",
            len(access_to_cdm_id),
            id_path,
        )
    except Exception:
        logger.exception("Couldn't load identifiers table from %s", id_path)

    return access_to_cdm_id


# ================================ PARSERS =================================
def parse_identifiers(entry, cdm_id: str) -> List[dict]:
    out = parse_identifiers_generic(entry=entry, xpath="ns:accession", prefix="UniProt", ns=NS)
    for row in out:
        row["entity_id"] = cdm_id
    return out


def _make_name_record(cdm_id: str, name_text: str, description: str) -> dict:
    return {
        "entity_id": cdm_id,
        "name": name_text,
        "description": description,
        "source": "UniProt",
    }


def parse_names(entry, cdm_id: str) -> List[dict]:
    names: List[dict] = []

    top_level_names = find_all_text(entry, "ns:name", NS)
    for txt in top_level_names:
        names.append(_make_name_record(cdm_id, txt, "UniProt protein name"))

    protein = entry.find("ns:protein", NS)
    if protein is not None:
        name_blocks_spec = [
            ("recommendedName", "recommended"),
            ("alternativeName", "alternative"),
        ]
        length_spec = [("fullName", "full"), ("shortName", "short")]

        for tag_name, logical_type in name_blocks_spec:
            for name_block in protein.findall(f"ns:{tag_name}", NS):
                for xml_tag, length_label in length_spec:
                    elem = name_block.find(f"ns:{xml_tag}", NS)
                    text = get_text(elem)
                    if not text:
                        continue
                    desc = f"UniProt {logical_type} {length_label} name"
                    names.append(_make_name_record(cdm_id, text, desc))

    return names


def parse_protein_info(entry, cdm_id: str) -> Optional[dict]:
    protein_info: dict = {}

    protein = entry.find("ns:protein", NS)
    if protein is not None:
        ec_paths = ["ns:recommendedName/ns:ecNumber", "ns:alternativeName/ns:ecNumber"]
        ec_numbers: List[str] = []
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
            clean_dict(
                {
                    "length": get_attr(seq_elem, "length"),
                    "mass": get_attr(seq_elem, "mass"),
                    "checksum": get_attr(seq_elem, "checksum"),
                    "modified": get_attr(seq_elem, "modified"),
                    "sequence_version": get_attr(seq_elem, "version"),
                    "sequence": get_text(seq_elem),
                }
            )
        )

    entry_modified = get_attr(entry, "modified") or get_attr(entry, "updated")
    if entry_modified:
        protein_info["entry_modified"] = entry_modified

    return protein_info if protein_info else None


def parse_evidence_map(entry) -> Dict[str, dict]:
    evidence_map: Dict[str, dict] = {}

    for ev in entry.findall("ns:evidence", NS):
        key = get_attr(ev, "key")
        if not key:
            continue

        evidence_type = get_attr(ev, "type")
        pubs: List[str] = []
        others: List[str] = []

        source = ev.find("ns:source", NS)
        if source is not None:
            raw_pubs, raw_others = parse_db_references(source, NS)

            normalized_pubs: List[str] = []
            for p in raw_pubs:
                up = p.upper()
                if up.startswith("PUBMED:"):
                    _, acc = p.split(":", 1)
                    normalized_pubs.append(f"PMID:{acc}")
                else:
                    normalized_pubs.append(p)

            pubs = normalized_pubs
            others = raw_others

        evidence_map[key] = clean_dict(
            {
                "evidence_type": evidence_type,
                "publications": pubs or None,
                "supporting_objects": others or None,
            }
        )

    return evidence_map


def _make_association(
    cdm_id: str,
    obj: str,
    predicate: Optional[str] = None,
    evidence_key: Optional[str] = None,
    evidence_map: Optional[dict] = None,
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


def parse_reaction_association(reaction, cdm_id: str, evidence_map: Dict[str, dict]) -> List[dict]:
    associations: List[dict] = []
    for dbref in reaction.findall("ns:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        if not db_type or not db_id:
            continue
        assoc = {
            "subject": cdm_id,
            "predicate": "catalyzes",
            "object": f"{db_type}:{db_id}",
            "evidence_type": None,
            "supporting_objects": None,
            "publications": None,
        }
        evidence_key = reaction.get("evidence")
        if evidence_key and evidence_key in evidence_map:
            assoc.update(evidence_map[evidence_key])
        associations.append(clean_dict(assoc))
    return associations


def parse_cofactor_association(cofactor, cdm_id: str) -> List[dict]:
    associations: List[dict] = []
    for dbref in cofactor.findall("ns:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        if not db_type or not db_id:
            continue
        assoc = {
            "subject": cdm_id,
            "predicate": "requires_cofactor",
            "object": f"{db_type}:{db_id}",
            "evidence_type": None,
            "supporting_objects": None,
            "publications": None,
        }
        associations.append(clean_dict(assoc))
    return associations


def parse_associations(entry, cdm_id: str, evidence_map: Dict[str, dict]) -> List[dict]:
    """
    removed annotation GO-2x
    We only keep:
      - taxonomy association
      - generic dbReference cross-references
      - catalytic activity / cofactor associations
    """
    associations: List[dict] = []

    # Taxonomy association
    organism = entry.find("ns:organism", NS)
    if organism is not None:
        taxon_ref = organism.find('ns:dbReference[@type="NCBI Taxonomy"]', NS)
        if taxon_ref is not None:
            tax_id = taxon_ref.get("id")
            if tax_id:
                associations.append(_make_association(cdm_id, f"NCBITaxon:{tax_id}"))

    # Generic dbReference cross-references
    for dbref in entry.findall("ns:dbReference", NS):
        db_type = dbref.get("type")
        db_id = dbref.get("id")
        if not db_type or not db_id:
            continue
        evidence_key = dbref.get("evidence")
        associations.append(
            _make_association(
                cdm_id,
                f"{db_type}:{db_id}",
                evidence_key=evidence_key,
                evidence_map=evidence_map,
            )
        )

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


def parse_publications(entry) -> List[str]:
    publications: List[str] = []
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

    # dedup while preserving order
    return list(dict.fromkeys(publications))


def parse_uniprot_entry(
    entry,
    cdm_id: str,
    current_timestamp: str,
    datasource_name: str = "UniProt import",
    prev_created: Optional[str] = None,
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
        "publications": parse_publications(entry),
    }


# ================================ SCHEMA =================================
schema_entities = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("entity_type", StringType(), False),
        StructField("data_source", StringType(), False),
        StructField("created", StringType(), True),
        StructField("updated", StringType(), True),
        StructField("version", StringType(), True),
        StructField("uniprot_created", StringType(), True),
        StructField("uniprot_modified", StringType(), True),
    ]
)

schema_identifiers = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("identifier", StringType(), False),
        StructField("source", StringType(), True),
        StructField("description", StringType(), True),
    ]
)

schema_proteins = StructType(
    [
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
    ]
)

schema_names = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("name", StringType(), False),
        StructField("description", StringType(), True),
        StructField("source", StringType(), True),
    ]
)

schema_associations = StructType(
    [
        StructField("subject", StringType(), True),
        StructField("object", StringType(), True),
        StructField("predicate", StringType(), True),
        StructField("evidence_type", StringType(), True),
        StructField("supporting_objects", ArrayType(StringType()), True),
        StructField("publications", ArrayType(StringType()), True),
    ]
)

schema_publications = StructType(
    [
        StructField("entity_id", StringType(), False),
        StructField("publication", StringType(), True),
    ]
)


# ================================ DELTA WRITE =================================
def save_batches_to_delta(
    spark: SparkSession,
    tables: Dict[str, Tuple[list, StructType]],
    output_dir: str,
    namespace: str,
    mode: str = "append",
) -> None:
    """
    Persist batches of parsed records for each CDM table into Delta Lake format,
    and register them into Spark SQL as external Delta tables under <namespace>.

    Paths:
      <output_dir>/<namespace>/<table_name>
    """

    spark.sql(f"CREATE DATABASE IF NOT EXISTS {namespace}")

    for table_name, (records, schema) in tables.items():
        if not records:
            continue

        delta_dir = os.path.abspath(os.path.join(output_dir, namespace, table_name))

        logger.debug(
            "Saving table %s.%s to %s (mode=%s, records=%d)",
            namespace,
            table_name,
            delta_dir,
            mode,
            len(records),
        )

        try:
            df = spark.createDataFrame(records, schema)

            writer = df.write.format("delta").mode(mode)
            if mode == "overwrite":
                writer = writer.option("overwriteSchema", "true")
            else:
                writer = writer.option("mergeSchema", "true")

            writer.save(delta_dir)

            spark.sql(
                f"""
                CREATE TABLE IF NOT EXISTS {namespace}.{table_name}
                USING DELTA
                LOCATION '{delta_dir}'
                """
            )
        except Exception:
            logger.exception("Failed to save '%s' to Delta at %s", table_name, delta_dir)
            raise


# ================================ MAIN PARSE LOOP =================================
def parse_entries(
    local_xml_path: str,
    target_date: Optional[str],
    batch_size: int,
    spark: SparkSession,
    tables: Dict[str, Tuple[list, StructType]],
    output_dir: str,
    namespace: str,
    current_timestamp: str,
    accession_to_entity_id: Dict[str, str],
    mode: str,
) -> Tuple[int, int]:
    """
    Parse UniProt XML entries and write to Delta in batches.
    Returns (processed_entry_count, skipped_entry_count).
    """

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

            # ----------------------------
            # Prevent entity_id duplication:
            # Reuse from existing identifiers table if present
            # Else stable UUIDv5(accession)
            # ----------------------------
            cdm_id = accession_to_entity_id.get(accession)
            prev_created = None

            if not cdm_id:
                cdm_id = stable_cdm_id_from_accession(accession)
            else:
                prev_created = None

            record = parse_uniprot_entry(entry_elem, cdm_id, current_timestamp, prev_created=prev_created)

            tables["entities"][0].append(record["entity"])
            tables["identifiers"][0].extend(record["identifiers"])
            tables["names"][0].extend(record["names"])

            if record["protein"]:
                tables["proteins"][0].append(record["protein"])

            tables["associations"][0].extend(record["associations"])
            tables["publications"][0].extend(
                {"entity_id": record["entity"]["entity_id"], "publication": pub} for pub in record["publications"]
            )

            entry_count += 1

            if entry_count % batch_size == 0:
                save_batches_to_delta(spark, tables, output_dir, namespace, mode=mode)
                for v in tables.values():
                    v[0].clear()
                logger.info("Processed and saved %d entries...", entry_count)

        except Exception:
            logger.exception("Error parsing UniProt entry, skipping")
            skipped += 1
            continue

    save_batches_to_delta(spark, tables, output_dir, namespace, mode=mode)
    return entry_count, skipped


def ingest_uniprot(
    xml_url: str,
    output_dir: str,
    namespace: str,
    target_date: Optional[str] = None,
    batch_size: int = 5000,
    mode: str = "append",
) -> None:
    """
    Ingest UniProt XML into Delta tables under <namespace>.

    mode:
      - append: incremental append (mergeSchema enabled)
      - overwrite: replace existing Delta data (overwriteSchema enabled)
    """
    current_timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()

    local_xml_path = prepare_local_xml(xml_url, output_dir)
    save_datasource_record(xml_url, output_dir)

    spark = get_spark_session(namespace)

    # Load existing accession -> entity_id mappings
    accession_to_entity_id = load_existing_identifiers(spark, output_dir, namespace)

    # In-memory batch buffers
    entities, identifiers, names, proteins, associations, publications = (
        [],
        [],
        [],
        [],
        [],
        [],
    )

    tables = {
        "entities": (entities, schema_entities),
        "identifiers": (identifiers, schema_identifiers),
        "names": (names, schema_names),
        "proteins": (proteins, schema_proteins),
        "associations": (associations, schema_associations),
        "publications": (publications, schema_publications),
    }

    logger.info(
        "Starting UniProt ingestion: xml=%s | namespace=%s | mode=%s | batch_size=%d",
        xml_url,
        namespace,
        mode,
        batch_size,
    )

    # ------------------------------------------------------------------
    # Main parsing + write loop
    # ------------------------------------------------------------------
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
        mode=mode,
    )

    logger.info(
        "Completed parsing UniProt XML. processed=%d skipped=%d",
        entry_count,
        skipped,
    )

    # ------------------------------------------------------------------
    # Verification
    # ------------------------------------------------------------------
    logger.info("Verifying Delta tables in namespace `%s`", namespace)

    spark.sql(f"SHOW TABLES IN {namespace}").show(truncate=False)

    table_names = [
        "entities",
        "identifiers",
        "names",
        "proteins",
        "associations",
        "publications",
    ]

    for tbl in table_names:
        logger.info("Verifying table: %s.%s", namespace, tbl)

        # Row count
        spark.sql(f"SELECT COUNT(*) AS row_count FROM {namespace}.{tbl}").show(truncate=False)

        # Small preview
        spark.sql(f"SELECT * FROM {namespace}.{tbl} LIMIT 5").show(truncate=False)

    spark.stop()

    logger.info(
        "All Delta tables successfully created and registered in Spark SQL under `%s`.",
        namespace,
    )


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
    help="Delta write mode",
)
def main(xml_url, output_dir, namespace, target_date, batch_size, mode):
    ingest_uniprot(
        xml_url=xml_url,
        output_dir=output_dir,
        namespace=namespace,
        target_date=target_date,
        batch_size=int(batch_size),
        mode=mode,
    )


if __name__ == "__main__":
    main()
