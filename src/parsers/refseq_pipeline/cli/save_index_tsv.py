import os
import click
from datetime import datetime
from refseq_pipeline.core.config import REFSEQ_ASM_REPORTS
from refseq_pipeline.core.refseq_io import download_text


def save_assembly_index(destination_path: str) -> bool:
    """
    Download the RefSeq assembly index file and save it to the given path.

    Args:
        destination_path (str): Full output file path where the index will be saved.

    Returns:
        bool: True if saved successfully, False otherwise.
    """
    try:
        # Download the index text from RefSeq FTP
        txt = download_text(REFSEQ_ASM_REPORTS)

        # Save the text content to a local file
        with open(destination_path, "w", encoding="utf-8") as f:
            f.write(txt)

        print(f"[index] saved to: {destination_path}")
        return True
    
    except Exception as e:
        print(f"[index] failed to download or save: {e}")
        return False


@click.command()

@click.option(
    "--out-dir",
    required=True,
    help="Directory to save the downloaded RefSeq assembly index."
    )

@click.option(
    "--tag",
    default=datetime.today().strftime("%Y%m%d"),
    show_default=True,
    help="Filename tag used to version the saved index."
    )


def main(out_dir: str, tag: str):
    """
    CLI entrypoint to download and save RefSeq assembly index.

    Example usage:
        python -m refseq_pipeline.cli.save_index_tsv \
            --out-dir bronze/refseq/indexes \
            --tag 20250930
    """
    # Ensure output directory exists
    os.makedirs(out_dir, exist_ok=True)

    # Construct full path with tag
    out_path = os.path.join(out_dir, f"assembly_summary_refseq.{tag}.tsv")

    # Call the core function to download & save
    success = save_assembly_index(out_path)

    # Exit with error code if failed (useful in pipelines)
    if not success:
        raise SystemExit("[fatal] Assembly index save failed.")


if __name__ == "__main__":
    main()
