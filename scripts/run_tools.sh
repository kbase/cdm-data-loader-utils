#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <input_json_file> <output_directory>"
    exit 1
fi

INPUT_JSON=$1       # Path to the input JSON file
OUTPUT_DIR=$2       # Directory where the results will be stored
PYTHON_SCRIPT="scripts/parse_index.py"  # Python script to parse the JSON file

# Ensure the output directory exists
mkdir -p "$OUTPUT_DIR"

# Check if Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script $PYTHON_SCRIPT not found."
    exit 1
fi

# Parse the JSON file and capture the output
PARSED_PATHS=$(python3 "$PYTHON_SCRIPT" "$INPUT_JSON" fna)
if [ $? -ne 0 ]; then
    echo "Error parsing JSON file."
    exit 1
fi

# Convert the parsed paths into a comma-separated list
IFS=$'\n' read -r -d '' -a PATHS_ARRAY < <(printf '%s\0' "$PARSED_PATHS")
PATHS_COMMA=$(IFS=,; echo "${PATHS_ARRAY[*]}")

# Check if there are paths to process
if [ -z "$PATHS_COMMA" ]; then
    echo "No paths found in the JSON file."
    exit 1
fi

# ensure conda environment is active and checkm2 db is downloaded
conda env create -f $BASE_DIR/env.yml
conda activate genome_loader_env
checkm2 database --download

# Generate output file paths
STATS_OUTPUT="$OUTPUT_DIR/stats.json"
CHECKM2_OUTPUT="$OUTPUT_DIR/checkm2/"

mkdir -p "$CHECKM2_OUTPUT"

# Run statswrapper.sh and save the results to the output directory
echo "Running statswrapper.sh..."
statswrapper.sh in="$PATHS_COMMA" --format=8 > "$STATS_OUTPUT"
if [ $? -ne 0 ]; then
    echo "Error running statswrapper.sh."
    exit 1
fi

# Run checkm2 and save the results to the output directory
echo "Running checkm2..."
checkm2 predict --threads $NUM_THREADS_CHECKM2_RUN --input "${PATHS_ARRAY[@]}" -o "$CHECKM2_OUTPUT"
if [ $? -ne 0 ]; then
    echo "Error running checkm2."
    exit 1
fi

echo "Process completed successfully. Results saved to $OUTPUT_DIR."
