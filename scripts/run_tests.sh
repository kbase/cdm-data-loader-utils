#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(dirname "$(dirname "$(readlink -f "$0")")")"
cd "$SCRIPT_DIR"

SYSTEM_PYTHON=/opt/conda/bin/python
VENV="$SCRIPT_DIR/docker_env"
CONSTRAINTS="$(mktemp)"
trap 'rm -f "$CONSTRAINTS"' EXIT

# ensure that these package(s) match the system versions
"$SYSTEM_PYTHON" -I - > "$CONSTRAINTS" <<'PY'
from importlib.metadata import version

for name in ["pyspark"]:
    print(f"{name}=={version(name)}")
PY

# use system python but install own dependencies
uv venv --python "$SYSTEM_PYTHON" "$VENV"

uv pip install --python "$VENV/bin/python" --constraint "$CONSTRAINTS" -e . --group dev

"$VENV/bin/python" -m pytest -m "not requires_ceph" --cov=src --cov-report=xml "$@"
