"""
Prefix Mapping Validation Script

python3 src/cdm_data_loader_utils/parsers/prefix_mapping.py
"""

import json
from pathlib import Path
from collections import Counter


REGISTRY_PATH = Path("registry.json")
MAPPING_PATH = Path("dataset/uniprot_prefix_remapping.json")


# ------------------------------------------------------------------
# Loaders
# ------------------------------------------------------------------


def load_registry(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def load_mapping(path: Path) -> list:
    with open(path) as f:
        return json.load(f)


# ------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------


def analyze_mapping_structure(mapping: list) -> None:
    print(f"Total mapping entries: {len(mapping)}")

    all_keys = set()
    for row in mapping:
        all_keys.update(row.keys())

    print("Keys present in mapping file:")
    print(all_keys)


def analyze_status_distribution(mapping: list) -> None:
    status_counts = Counter(row.get("_status") for row in mapping)

    for status, count in sorted(status_counts.items(), key=lambda x: str(x[0])):
        print(f"{status}: {count}")


def extract_canonical_targets(mapping: list) -> set:
    """
    Collect all canonical prefixes referenced in the remapping file.
    """

    canonical_targets = set()

    for row in mapping:
        # Single match case
        match = row.get("match")
        if match:
            canonical_targets.add(match.strip().lower())

        # Multiple matches case
        matches = row.get("matches")
        if matches:
            for m in matches:
                canonical_targets.add(m.strip().lower())

    return canonical_targets


def validate_against_registry(canonical_targets: set, registry: dict) -> None:
    registry_prefixes = set(registry.keys())
    invalid = canonical_targets - registry_prefixes

    if invalid:
        print(f"\n{len(invalid)} targets that do not appear in egistry")
        for item in sorted(invalid):
            print(f"{item}")
    else:
        print("\nAll referenced canonical targets are present in the registry.")


def analyze_comments(mapping: list) -> None:
    unique_comments = set()

    for row in mapping:
        comments = row.get("comment")
        if comments:
            for c in comments:
                unique_comments.add(c.strip())

    for comment in sorted(unique_comments):
        print(f"{comment}")


def inspect_none_status(mapping: list) -> None:
    none_rows = [row for row in mapping if row.get("_status") is None]
    print(f"\nNumber of entries without status: {len(none_rows)}")

    if not none_rows:
        return

    for row in none_rows[:10]:
        print(f"{row.get('__prefix')}")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    print("\nBioregistry JSON:")
    registry = load_registry(REGISTRY_PATH)

    print("\nUniProt prefix remapping file:")
    mapping = load_mapping(MAPPING_PATH)

    analyze_mapping_structure(mapping)
    analyze_status_distribution(mapping)

    canonical_targets = extract_canonical_targets(mapping)
    validate_against_registry(canonical_targets, registry)

    analyze_comments(mapping)
    inspect_none_status(mapping)


if __name__ == "__main__":
    main()
