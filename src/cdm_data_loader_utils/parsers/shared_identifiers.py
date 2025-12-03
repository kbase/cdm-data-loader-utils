from xml_utils import get_text

def parse_identifiers_generic(entry, xpath, prefix, ns):
    result = []
    for node in entry.findall(xpath, ns):
        text = get_text(node)
        if not text:
            continue
        result.append({
            "identifier": f"{prefix}:{text}",
            "source": prefix,
            "description": f"{prefix} accession"
        })
    return result