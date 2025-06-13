import requests
import os
import xml.etree.ElementTree as ET
import uuid
from datetime import datetime, timezone
import json
import glob
import pandas as pd
import time

"""
pipeline:
Paging to capture cluster id → batch download XML → parsing XML structure → output in CDM format
"""

def fetch_uniref100_ids(
        uniref_level="100",
        out_tsv="uniref100_ids_part.tsv",
        max_count=1000,
        page_size=100,
        max_retry=5
        ):
    
    """
    Fetch UniRef100 cluster IDs using the UniProt REST API with automatic paging.

    Args:
    uniref_level (str): UniRef level, default "100".
    out_tsv (str): Output TSV filename.
    max_ids (int): Maximum number of IDs to fetch.
    page_size (int): Number of IDs per page.
    max_retry (int): Maximum retry count per request.

    Returns:
    list: List of UniRef100 cluster IDs.

    Purpose: 
    - Every time when I ran the fetch_uniref100_ids function, I am grabbing the IDs of all the latest version of UniRef100 clusters in real time via UniProt's REST API.
    - Whenever UniProt official data modifies or deletes any clusters, the next run will get a list of the latest IDs.

    """
    
    # Base URL for UniRef search API
    url = "https://rest.uniprot.org/uniref/search" 
    
    parameters = {
        "query": f"identity:{uniref_level}.0",  ## Query condition, identity level
        "format": "tsv", ## The output format is a tsv table
        "fields": "id",  ## Only ID 
        "size": page_size 
        }
    
    ## Stores all cluster IDs, use cursor parameters for API paging
    ids = []
    cursor = None
    
    ## Automatically page through to get all IDs
    while True:
        if cursor:
            ## If it's not the first page, add the cursor parameter to specify the next page
            parameters["cursor"] = cursor 
        else:
            ## Make sure the first page doesn't have a cursor parameter
            parameters.pop("cursor", None) 

        ## Automatically retry when encountering network anomalies and occasional API errors
        for retry in range(max_retry):
            try:
                r = requests.get(url, params=parameters, timeout=10) 
            except Exception as e:
                print(f"Exception: {e}, retry {retry+1}") ## avoid network and web service issues
                time.sleep(3)
                continue
            if r.ok:
                break
            print(f"Error: {r.status_code}, retry {retry+1}")
            time.sleep(3)
        else:
            print("Failed after retries, aborting fetch.")
            break
        
        ## Parses the contents returned by the API, ensure one cluster id per line 
        ## original text content returned by the API and split each line by removing the whitespace
        lines = r.text.strip().split("\n") 
        if lines and lines[0].startswith("Cluster ID"):
            lines = lines[1:] ## skip the first line of Cluster ID header

        ## Ensure that only real cluster IDs are added to the final result
        ids.extend([line.strip() for line in lines if line.strip()])

        ## Determines if the number of IDs collected reaches maximum number. 
        ## If reached, truncate to the specified maximum number and end the while loop.
        if len(ids) >= max_count:
            ids = ids[:max_count]
            break

        ## Automate multiple pages of data and crawling if the API has the next page
        next_link = r.links.get("next", {}).get("url", None)
        if not next_link:
            break
        
        ## Check if the URL of the next page contains the cursor parameter. 
        ## If it does, extract the cursor value and use it as a paging parameter in the next request.
        if "cursor=" in next_link:
            cursor = next_link.split("cursor=")[-1]
        else: 
            print("Cursor noe found in the next link")

        time.sleep(3)
    
    ## Write all collected Cluster IDs to out_tsv
    with open(out_tsv, "w") as fout:
        fout.write("Cluster ID\n")
        for cid in ids:
            fout.write(f"{cid}\n")
    print(f"Written {len(ids)} IDs to {out_tsv}")

    return ids


def batch_download_uniref_xml(ids, out_dir=".", max_retry=3):
    """
    batch_download_uniref_xml will download only what is not available locally based on the latest ID list, and skip what already exists locally with content. Newly added data will be downloaded and old data will not be downloaded again.

    """
    for idx, cid in enumerate(ids, 1):
        filename = os.path.join(out_dir, f"{cid}.xml")
        if os.path.exists(filename) and os.path.getsize(filename) > 0:
            continue

        url = f"https://rest.uniprot.org/uniref/{cid}.xml"

        for retry in range(max_retry):
            try:
                r = requests.get(url, timeout=10)
                if r.ok:
                    with open(filename, "wb") as fout:
                        fout.write(r.content)
                    break
                else:
                    print(f"[{idx}/{len(ids)}] Failed {cid}, status:{r.status_code}, retry {retry+1}")
            except Exception as e:
                print(f"[{idx}/{len(ids)}] Exception: {e}, retry {retry+1}")
            time.sleep(3)
        else:
            print(f"[{idx}/{len(ids)}] Totally failed to download {cid}")


"""
For each UniRef cluster <entry>, the core information falls into three main categories:
- Cluster-level attributes (e.g. cluster id, update time, annotations, public species, number of members, etc., directly under the <entry> tag)

- representativeMember, usually only one, contains all the attributes of the representative sequence and representative protein of the cluster.

- Member list, there can be many members, each member is a <member> block, which contains various attributes of the specific protein.


"""

def parse_representative(entry, ns):
    """

    Accurately extracts all information in <representativeMember> with clear field mapping and robust exception handling.

    entry:
        entry: the currently parsed <entry> element
        ns: namespace dictionary (typically {“u”: “http://uniprot.org/uniref”})

    output:
        rep_info: dict with all normalized fields representing members

    """

    # Initialization returns a dictionary, all target fields are filled with None by default
    rep_info = {
        "type": None, "id": None, "source_organism": None, "taxon": None,
        "protein_name": None, "uniparc_id": None, "uniref90_id": None,
        "uniref50_id": None, "uniprotkb_accession": None, "length": None,
        "is_seed": None, "sequence": None, "seq_length": None, "checksum": None
    }

    # Mapping table of UniRef property names to CDM standard fields
    field_map = {
        "source organism": "source_organism",
        "NCBI taxonomy": "taxon",
        "protein name": "protein_name",
        "UniParc ID": "uniparc_id",
        "UniRef90 ID": "uniref90_id",
        "UniRef50 ID": "uniref50_id",
        "UniProtKB accession": "uniprotkb_accession",
        "length": "length",
        "isSeed": "is_seed"
    }

    # Locate the <dbReference> node under <representativeMember>
    rep = entry.find("u:representativeMember/u:dbReference", ns)
    if rep is not None:
        # 1. Populate type and id (directly from attributes)
        rep_info["type"] = rep.attrib.get("type")
        rep_info["id"] = rep.attrib.get("id")

        # 2. Iterate over all <property>, use field_map for field alignment
        for prop in rep.findall("u:property", ns):
            ptype = prop.attrib.get("type")
            value = prop.attrib.get("value")
            if ptype in field_map:
                key = field_map[ptype]

                # Doing conversions based on field types
                if key == "taxon" or key == "length":
                    try:
                        rep_info[key] = int(value)     
                    except Exception:
                        rep_info[key] = value          
                elif key == "is_seed":
                    rep_info[key] = value.lower() == "true"   
                else:
                    rep_info[key] = value                    

        # 3. Parses <sequence> information on behalf of member
        seq_elem = entry.find("u:representativeMember/u:sequence", ns)
        if seq_elem is not None:
            rep_info["sequence"] = seq_elem.text                             
            try:
                rep_info["seq_length"] = int(seq_elem.attrib.get("length", 0))  
            except Exception:
                rep_info["seq_length"] = seq_elem.attrib.get("length")          
            rep_info["checksum"] = seq_elem.attrib.get("checksum")          

    return rep_info


def parse_members(entry, ns):
    # Initialize Return List
    members = []

    # Defining the mapping of XML attribute names to CDM field names
    field_map = {
        "source organism": "source_organism",
        "NCBI taxonomy": "taxon",
        "protein name": "protein_name",
        "UniParc ID": "uniparc_id",
        "UniRef90 ID": "uniref90_id",
        "UniRef50 ID": "uniref50_id",
        "UniProtKB accession": "uniprotkb_accession",
        "length": "length",
        "isSeed": "is_seed"
    }

    # Iterate through all <member><dbReference> nodes
    for member in entry.findall("u:member/u:dbReference", ns):
        # Initialize the information structure of a single member
        mem_info = {
            "type": member.attrib.get("type"),
            "id": member.attrib.get("id"),
            "source_organism": None, "taxon": None, "protein_name": None,
            "uniparc_id": None, "uniref90_id": None, "uniref50_id": None,
            "uniprotkb_accession": None, "length": None, "is_seed": None
        }

        # Iterate over <property> tags to extract field information
        for prop in member.findall("u:property", ns):
            ptype = prop.attrib.get("type")
            value = prop.attrib.get("value")
            if ptype in field_map:
                key = field_map[ptype]

                if key in ("taxon", "length"):
                    try:
                        mem_info[key] = int(value)
                    except ValueError:
                        mem_info[key] = value
                elif key == "is_seed":
                    mem_info[key] = value.lower() == "true"
                else:
                    mem_info[key] = value
            else:
                # If an unexpected field appears, print a debug message
                print(f"Unexpected property: {ptype} = {value}")

        members.append(mem_info)
    return members


def get_entry_property(entry, ns, property_name, as_int=False):
    # Iterate over all <property> tags under the entry node
    for prop in entry.findall("u:property", ns):
        # If the type attribute of property is equal to property_name
        if prop.attrib.get("type") == property_name:
            value = prop.attrib.get("value")
            # If as_int=True, convert the string value to int
            if as_int:
                try:
                    return int(value)
                except Exception:
                    return value
            return value
    return None


def parse_uniref_xml_to_cdm(xml_file):
    ns = {"u": "http://uniprot.org/uniref"}
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
    except Exception as e:
        print(f"Error parsing {xml_file}: {e}")
        return []

    cdm_results = []
    for entry in root.findall("u:entry", ns):
        cluster_id = entry.attrib.get("id")
        updated = entry.attrib.get("updated")
        uniref_level = cluster_id.split("_")[0].replace("UniRef", "") if cluster_id else None
        name_elem = entry.find("u:name", ns)
        functional_annotation = name_elem.text if name_elem is not None else None
        member_count = get_entry_property(entry, ns, "member count", as_int=True)
        common_taxon = get_entry_property(entry, ns, "common taxon")
        common_taxon_id = get_entry_property(entry, ns, "common taxon ID", as_int=True)
        rep_info = parse_representative(entry, ns)
        members = parse_members(entry, ns)
        now = datetime.now(timezone.utc).isoformat()

        cdm = {
            "entity_id": f"CDM:{uuid.uuid4()}",
            "entity_type": "protein_cluster",
            "source": "UniRef",
            "uniref_level": uniref_level,
            "cluster_id": cluster_id,
            "updated": updated,
            "member_count": member_count,
            "functional_annotation": functional_annotation,
            "common_taxon": common_taxon,
            "common_taxon_id": common_taxon_id,
            "representative": rep_info,
            "members": members,
            "created": now,
            "cdm_version": "1.0",  # 可以考虑作为常量
        }
        cdm_results.append(cdm)
    return cdm_results
    

def flatten_cdm_entry(cdm):
    """
    Flat CDM structure, one row per member, retain clustering attributes + all attributes for that member
    """
    basic = {
        "entity_id": cdm["entity_id"],
        "entity_type": cdm["entity_type"],
        "source": cdm["source"],
        "uniref_level": cdm["uniref_level"],
        "cluster_id": cdm["cluster_id"],
        "updated": cdm["updated"],
        "member_count": cdm["member_count"],
        "functional_annotation": cdm["functional_annotation"],
        "common_taxon": cdm["common_taxon"],
        "common_taxon_id": cdm["common_taxon_id"],
        "created": cdm["created"],
        "cdm_version": cdm["cdm_version"],
    }

    rep = cdm.get("representative", {})
    rep_flat = {f"rep_{k}": v for k, v in rep.items()}
    records = []

    for i in cdm.get("members", []):
        mem_flat = {f"mem_{k}": v for k, v in i.items()}
        row = {**basic, **rep_flat, **mem_flat}
        records.append(row)
    if not records:
        row = {**basic, **rep_flat}
        records.append(row)
    return records


def main():
    tsv_path = "uniref100_ids_part.tsv"
    max_count = 1000
    page_size = 100
    
    ## Get a list of UniRef100 group IDs
    if not os.path.exists(tsv_path):
        fetch_uniref100_ids(out_tsv=tsv_path, max_ids=max_count, page_size=page_size)
    
    ## Read the list of IDs and download all XMLs
    with open(tsv_path) as fin:
        ids = [line.strip() for line in fin if line.strip() and not line.startswith("Cluster ID")]
    batch_download_uniref_xml(ids)
    
    ## Parses all XML to CDM
    all_cdm_entries = []
    xml_files = glob.glob("UniRef100_*.xml")
    
    for xml_file in xml_files:
        try:
            all_cdm_entries.extend(parse_uniref_xml_to_cdm(xml_file))
            print(f"Parsed {xml_file}")
        except Exception as e:
            print(f"Error parsing {xml_file}: {e}")
    
    ## Export as jsonl
    with open("all_uniref100_cdm.jsonl", "w") as fout:
        for cdm in all_cdm_entries:
            fout.write(json.dumps(cdm, ensure_ascii=False) + "\n")
    print("All clusters written to all_uniref100_cdm.jsonl")
    
    ## Spread export to csv and parquet
    all_flat = []
    for cdm in all_cdm_entries:
        all_flat.extend(flatten_cdm_entry(cdm))

    df_flat = pd.DataFrame(all_flat)
    df_flat.to_csv("all_uniref100_flat.csv", index=False)
    df_flat.to_parquet("all_uniref100_flat.parquet", index=False)
    
    ## Exporting an example and full array
    if all_cdm_entries:
        with open("uniref100_cdm_example.json", "w") as fout:
            json.dump(all_cdm_entries[0], fout, indent=2, ensure_ascii=False)
        print("example cluster save as uniref100_cdm_example.json")

        with open("all_uniref100_cdm_array.json", "w") as fout:
            json.dump(all_cdm_entries, fout, indent=2, ensure_ascii=False)
        print("all clusters as JSON array save as all_uniref100_cdm_array.json")

if __name__ == "__main__":
    main()
