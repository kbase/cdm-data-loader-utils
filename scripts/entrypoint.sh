#!/usr/bin/env bash
set -euo pipefail

# Ensure at least one argument is provided
if [ "$#" -eq 0 ]; then
  echo "Usage: $0 {uniref|uniprot|test} [args...]"
  exit 1
fi

cmd="$1"
shift

case "$cmd" in
  uniref)
    # Run the uniref pipeline with any additional arguments via tini
    exec /usr/bin/tini -- uv run --no-sync uniref_pipeline "$@"
    ;;
  uniprot)
    # Run the uniprot pipeline with any additional arguments via tini
    exec /usr/bin/tini -- uv run --no-sync uniprot_pipeline "$@"
    ;;
  test)
    # run the tests
    exec /usr/bin/tini -- uv run --no-sync pytest -m "not requires_spark"
    ;;
  bash)
    exec /usr/bin/tini -- /bin/bash
    ;;
  *)
    echo "Error: unknown command '$cmd'; valid commands are 'uniref' or 'uniprot'." >&2
    exit 1
    ;;
esac
