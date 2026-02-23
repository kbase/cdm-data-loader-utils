#!/bin/bash
set -x

# Get the parent directory of the dir where the run_tests script is located
SCRIPT_DIR="$(dirname "$(dirname "$(readlink -f "$0")")")"
# Change the current working directory to SCRIPT_DIR
cd "$SCRIPT_DIR"

# use the system packages in this virtual environment
uv venv --system-site-packages
source $(SCRIPT_DIR)/.venv/bin/activate

# run the tests using the active venv and with the dev dependencies installed.
uv run --active --frozen --group dev pytest --cov=src --cov-report=xml
