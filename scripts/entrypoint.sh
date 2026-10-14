#!/usr/bin/env bash
set -euo pipefail

# Each pipeline name must match a uv script name, as it is run as `uv run <name>`.
PIPELINES=(all_the_bacteria ncbi_ftp_sync ncbi_rest_api uniprot uniref)
OTHER=(xml_split test ceph_integration_test bash)
VALID_COMMANDS=("${PIPELINES[@]}" "${OTHER[@]}")

DLT_CONF_SRC_DIR=/input_dir/dlt_conf
DLT_CONF_DEST_DIR=.dlt

usage() {
  local joined
  joined=$(IFS='|'; echo "${VALID_COMMANDS[*]}")
  echo "Usage: $0 {${joined}} [args...]" >&2
}

# Return 0 if the first argument matches any of the remaining arguments.
contains() {
  local needle="$1"
  shift
  local item
  for item in "$@"; do
    if [ "$item" = "$needle" ]; then
      return 0
    fi
  done
  return 1
}

# Copy any dlt config files from the input dir into ./.dlt; exit on any failure.
copy_dlt_config() {
  local file
  for file in config.toml secrets.toml; do
    if [ -e "${DLT_CONF_SRC_DIR}/${file}" ]; then
      mkdir -p "${DLT_CONF_DEST_DIR}" || {
        echo "Error: could not create '${DLT_CONF_DEST_DIR}' in $(pwd)." >&2
        exit 1
      }
      cp "${DLT_CONF_SRC_DIR}/${file}" "${DLT_CONF_DEST_DIR}/${file}" || {
        echo "Error: could not copy '${DLT_CONF_SRC_DIR}/${file}' to '${DLT_CONF_DEST_DIR}/${file}'." >&2
        exit 1
      }
    fi
  done
}

if [ "$#" -eq 0 ]; then
  usage
  exit 1
fi

cmd="$1"
shift

if contains "$cmd" "${PIPELINES[@]}"; then
  copy_dlt_config
  exec /usr/bin/tini -- uv run --no-sync "$cmd" "$@"
fi

case "$cmd" in
  xml_split)
    exec /usr/bin/tini -- xml_file_splitter "$@"
    ;;
  test)
    exec /usr/bin/tini -- uv run --no-sync pytest -m "not requires_spark and not requires_ceph"
    ;;
  ceph_integration_test)
    # run the CEPH-backed integration tests (requires a running CEPH instance)
    exec /usr/bin/tini -- uv run --no-sync pytest -m "requires_ceph" -v "$@"
    ;;
  bash)
    exec /usr/bin/tini -- /bin/bash
    ;;
  *)
    echo "Error: unknown command '$cmd'." >&2
    usage
    exit 1
    ;;
esac
