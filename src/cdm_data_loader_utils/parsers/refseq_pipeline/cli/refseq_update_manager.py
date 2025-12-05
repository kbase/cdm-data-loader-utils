import logging
from datetime import datetime
from pathlib import Path

import click

from cdm_data_loader_utils.parsers.refseq_pipeline.core.config import REFSEQ_ASSEMBLY_SUMMARY_URL
from cdm_data_loader_utils.parsers.refseq_pipeline.core.refseq_io import (
    download_text,
    normalize_multiline_text,
    text_sha256,
)

"""

PYTHONPATH=src/parsers \
python src/parsers/refseq_pipeline/cli/refseq_update_manager.py \
  --database refseq_api \
  --hash-table assembly_hashes \
  --kind md5 \
  --chunk-size 500 \
  --max-workers 12 \
  --fast-mode \
  --data-dir delta_data/refseq

"""


from cdm_data_loader_utils.parsers.refseq_pipeline.cli import detect_updates, snapshot_hashes

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def get_project_root() -> Path:
    try:
        return Path(__file__).resolve().parents[2]
    except NameError:
        return Path.cwd()


def get_latest_local_index(bronze_dir: Path) -> Path | None:
    files = sorted(bronze_dir.glob("assembly_summary_refseq.*.tsv"))
    return files[-1] if files else None


def parse_tag_from_filename(path: Path) -> str | None:
    name = path.name  # assembly_summary_refseq.20250930.tsv
    parts = name.split(".")
    if len(parts) >= 3:
        return parts[1]
    return None


@click.command()
@click.option(
    "--database",
    default="refseq_api",
    show_default=True,
    help="Spark SQL database where hash snapshots are stored.",
)
@click.option(
    "--hash-table",
    default="assembly_hashes",
    show_default=True,
    help="Delta table name for hash snapshots.",
)
@click.option(
    "--kind",
    type=click.Choice(["annotation", "md5"]),
    default="md5",
    show_default=True,
    help="Hashing strategy passed to snapshot_hashes.",
)
@click.option(
    "--tag",
    default=None,
    help="Tag for the new snapshot (default: today's date in YYYYMMDD).",
)
@click.option(
    "--force-snapshot/--no-force-snapshot",
    default=False,
    show_default=True,
    help="Force creating a new snapshot even if the index TSV did not change.",
)
@click.option(
    "--chunk-size",
    default=500,
    show_default=True,
    help="Chunk size passed to snapshot_hashes.",
)
@click.option(
    "--max-workers",
    default=12,
    show_default=True,
    help="Max worker threads passed to snapshot_hashes.",
)
@click.option(
    "--fast-mode/--no-fast-mode",
    default=True,
    show_default=True,
    help="Whether to use fast-mode in snapshot_hashes.",
)
@click.option(
    "--data-dir",
    default=None,
    help="Base directory for Delta data (default: <project_root>/delta_data/refseq).",
)
@click.option(
    "--output-dir",
    default=None,
    help="Directory to save CSV diff from detect_updates "
    "(default: <project_root>/diff_results/<old_tag>_to_<new_tag>).",
)
def main(
    database: str,
    hash_table: str,
    kind: str,
    tag: str | None,
    force_snapshot: bool,
    chunk_size: int,
    max_workers: int,
    fast_mode: bool,
    data_dir: str | None,
    output_dir: str | None,
):
    project_root = get_project_root()
    logger.info(f"[update] Project root resolved to: {project_root}")

    bronze_dir = project_root / "bronze" / "refseq" / "indexes"
    bronze_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"[update] Using bronze index dir: {bronze_dir}")

    if data_dir:
        data_dir_path = Path(data_dir)
    else:
        data_dir_path = project_root / "delta_data" / "refseq"
    logger.info(f"[update] Using Delta data dir: {data_dir_path}")

    logger.info(f"[update] Downloading latest assembly index from: {REFSEQ_ASSEMBLY_SUMMARY_URL}")
    raw_txt = download_text(REFSEQ_ASSEMBLY_SUMMARY_URL)
    latest_txt = normalize_multiline_text(raw_txt)
    latest_sha = text_sha256(latest_txt)
    logger.info(f"[update] Latest remote index SHA256: {latest_sha}")

    last_local = get_latest_local_index(bronze_dir)
    old_tag: str | None = None
    old_sha: str | None = None

    if last_local and last_local.exists():
        logger.info(f"[update] Found existing local index: {last_local}")
        old_tag = parse_tag_from_filename(last_local)
        try:
            old_txt = normalize_multiline_text(last_local.read_text(encoding="utf-8"))
            old_sha = text_sha256(old_txt)
            logger.info(f"[update] Last local index tag={old_tag}, SHA256={old_sha}")
        except Exception as e:
            logger.warning(f"[update] Failed to read previous index {last_local}: {e}")
    else:
        logger.info("[update] No previous local index found. This appears to be the first run.")

    if old_sha is not None and old_sha == latest_sha and not force_snapshot:
        logger.info("[update] Remote index is identical to the last local index.")
        logger.info("[update] No snapshot will be created (use --force-snapshot to override).")
        return

    new_tag = tag or datetime.utcnow().strftime("%Y%m%d")
    logger.info(f"[update] New snapshot tag will be: {new_tag}")

    new_index_path = bronze_dir / f"assembly_summary_refseq.{new_tag}.tsv"

    new_index_path.write_text(latest_txt, encoding="utf-8")
    logger.info(f"[update] Saved latest index to: {new_index_path}")

    logger.info("[update] Running snapshot_hashes to create new snapshot...")

    snapshot_args = [
        "--database",
        database,
        "--hash-table",
        hash_table,
        "--tag",
        new_tag,
        "--kind",
        kind,
        "--index-path",
        str(new_index_path),
        "--chunk-size",
        str(chunk_size),
        "--max-workers",
        str(max_workers),
        "--data-dir",
        str(data_dir_path),
    ]
    # fast-mode flag
    if fast_mode:
        snapshot_args.append("--fast-mode")
    else:
        snapshot_args.append("--no-fast-mode")

    logger.info(f"[update] snapshot_hashes args: {snapshot_args}")
    snapshot_hashes.main(args=snapshot_args, standalone_mode=False)

    if old_tag:
        if output_dir:
            output_dir_path = Path(output_dir)
        else:
            output_dir_path = project_root / "diff_results" / f"{old_tag}_to_{new_tag}"

        output_dir_path.mkdir(parents=True, exist_ok=True)

        logger.info("[update] Running detect_updates to compute snapshot diff")
        detect_args = [
            "--database",
            database,
            "--table",
            hash_table,
            "--old-tag",
            old_tag,
            "--new-tag",
            new_tag,
            "--output",
            str(output_dir_path),
        ]

        logger.info(f"[update] detect_updates args: {detect_args}")
        detect_updates.main(args=detect_args, standalone_mode=False)
        logger.info(f"[update] Diff CSV should be written under: {output_dir_path}")
    else:
        logger.info("[update] No previous tag found. Skipping detect_updates diff.")

    logger.info("[update] RefSeq update manager completed successfully.")


if __name__ == "__main__":
    main()
