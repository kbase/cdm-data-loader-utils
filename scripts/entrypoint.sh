#!/usr/bin/env bash
set -euo pipefail

# Ensure at least one argument is provided
if [ "$#" -eq 0 ]; then
  echo "Usage: $0 {uniref|uniprot} [args...]"
  exit 1
fi

cmd="$1"
shift

case "$cmd" in
  uniref)
    # Run the uniref pipeline with any additional arguments via tini
    exec /usr/bin/tini -- uv run uniref_pipeline "$@"
    ;;
  uniprot)
    # Run the uniprot pipeline with any additional arguments via tini
    exec /usr/bin/tini -- uv run uniprot_pipeline "$@"
    ;;
  *)
    echo "Error: unknown command '$cmd'; valid commands are 'uniref' or 'uniprot'." >&2
    exit 1
    ;;
esac
