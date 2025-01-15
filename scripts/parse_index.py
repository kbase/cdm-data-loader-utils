"""Script to parse a file containing paths to genome data."""

import json
import sys


def main() -> None:
    """Extracts values associated with a specified key from a JSON file and prints them.

    Usage:
        parse_json.py <input_json_file> <key>

    Arguments:
        input_json_file: Path to the input JSON file.
        key: The key whose values need to be extracted.

    Output:
        Prints the values associated with the specified key, one per line.
    """
    # Check if the correct number of arguments is provided
    if len(sys.argv) < 3:
        print("Usage: parse_json.py <input_json_file> <key>")
        sys.exit(1)

    input_file = sys.argv[1]  # Path to the JSON file
    target_key = sys.argv[2]  # Key to extract values for

    try:
        # Attempt to open and read the JSON file
        with open(input_file, "r") as file:
            data = json.load(file)
    except Exception as e:
        print(f"Error reading JSON file: {e}")
        sys.exit(1)

    # throw error if top level data structure is incorrect
    if not isinstance(data, dict):
        print("Error: expected JSON file to be a dictionary")
        sys.exit(1)

    # iterate over dictionary values and extract the values from dicts containing the target key
    extracted_paths = [
        entry[target_key]
        for entry in data.values()
        if isinstance(entry, dict) and target_key in entry
    ]

    # Print each extracted path on a new line
    for path in extracted_paths:
        print(path)


if __name__ == "__main__":
    main()
