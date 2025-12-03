#!/bin/bash
set -x

# Get the parent directory of the dir where the run_tests script is located
SCRIPT_DIR="$(dirname "$(dirname "$(readlink -f "$0")")")"
# Change the current working directory to SCRIPT_DIR
cd "$SCRIPT_DIR"

uv run pytest --cov=src --cov-report=xml --cov-append tests/
