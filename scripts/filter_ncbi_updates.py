"""Generate lists of NCBI assembly IDs to be updated, filtered by RefSeq membership and a list of assemblies of interest.

Update dates for the NCBI FTP site can be found in the file:
https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ftp_mod_times.txt.gz

Assemblies whose update date is later than the input date are selected.

By default, all RefSeq IDs are considered "assemblies of interest", and only Genbank assembly IDs are filtered. To
filter RefSeq IDs as well as Genbank, pass the flag --filter_refseq.
"""

import math
import re
import sys
from collections.abc import Iterable
from datetime import date, datetime
from enum import StrEnum
from pathlib import Path
from typing import Final

import pandas as pd
from pydantic import Field, FilePath, field_validator
from pydantic_settings import BaseSettings, CliApp, CliImplicitFlag, SettingsConfigDict

ASSEMBLY_DIR_RE: re.Pattern[str] = re.compile(r"^(GC[AF]_\d{9}\.\d+)_(.+)$")
BATCH_SIZE: Final[int] = 1000


class Source(StrEnum):
    """Data source."""

    REFSEQ = "refseq"
    GENBANK = "genbank"


class Period(StrEnum):
    """Positioning in the fourth dimension."""

    BEFORE = "before"
    AFTER = "after"


Splits = dict[tuple[Source, Period], pd.DataFrame]


class Config(BaseSettings):
    """Generate lists of NCBI assembly IDs to be updated, filtered by modification date and a list of interesting IDs.

    Files are output in batched sets of 1000, suitable for feeding to ncbi `datasets` or other commands.
    """

    model_config = SettingsConfigDict(
        cli_prog_name="filter_ncbi_updates",
        cli_kebab_case=True,
        cli_enforce_required=True,
        extra="forbid",
    )

    ftp_mod_times_path: FilePath = Field(
        default=Path("ftp_mod_times.txt"),
        description=(
            "NCBI FTP directory update times, downloaded from"
            "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/ftp_mod_times.txt.gz"
            "This is a TSV file without a header in the format"
            "assembly_dir<TAB>modification_date"
        ),
    )
    assemblies_of_interest_path: FilePath = Field(
        default=Path("assemblies.txt"),
        description="Single-column file of assembly IDs (no header).",
    )
    mod_date: date = Field(
        description="Modification date threshold, format YYYY/MM/DD.",
    )
    output_dir: Path = Field(
        description="Directory to write output files to.",
    )
    filter_refseq: CliImplicitFlag[bool] = Field(
        default=False,
        description=(
            "If set, filter RefSeq assembly IDs against the assembly list. "
            "By default, all RefSeq IDs are kept and only Genbank IDs are "
            "filtered against the assemblies of interest list."
        ),
    )
    batch_size: int = Field(
        default=BATCH_SIZE,
        description="Number of assembly IDs per output batch file.",
    )

    @field_validator("mod_date", mode="before")
    @classmethod
    def _parse_mod_date(cls, value: str):
        """Modification date parser."""
        if isinstance(value, date):
            return value
        return datetime.strptime(value, "%Y/%m/%d").date()  # noqa: DTZ007

    @property
    def derived_dir(self) -> Path:
        """Derived directory."""
        return self.output_dir / "derived"

    @property
    def mod_date_tag(self) -> str:
        """Modification date in fs-friendly format."""
        return self.mod_date.strftime("%Y%m%d")

    def source_dir(self, source: Source) -> Path:
        """Source directory."""
        return self.output_dir / source.value

    def ensure_dirs(self) -> None:
        """Ensure that all subdirectories required by the script have been created."""
        self.derived_dir.mkdir(parents=True, exist_ok=True)
        for source in Source:
            self.source_dir(source).mkdir(parents=True, exist_ok=True)

    def cli_cmd(self) -> None:
        """Entry point invoked by CliApp.run()."""
        run_pipeline(self)


def parse_assembly_dir(assembly_dir: str) -> tuple[str, str]:
    """Split 'GCF_000123456.1_SomeName' -> ('GCF_000123456.1', 'SomeName')."""
    match = ASSEMBLY_DIR_RE.match(assembly_dir)
    if not match:
        msg = f"assembly_dir '{assembly_dir}' does not match expected pattern"
        raise ValueError(msg)
    return match.group(1), match.group(2)


def infer_source(assembly_id: str) -> Source:
    """Infer whether this is a RefSeq or Genbank assembly ID."""
    if assembly_id.startswith("GCF_"):
        return Source.REFSEQ
    if assembly_id.startswith("GCA_"):
        return Source.GENBANK
    msg = f"Unrecognized assembly_id prefix: {assembly_id}"
    raise ValueError(msg)


def load_ftp_mod_times(path: Path) -> pd.DataFrame:
    """Load and fully parse the raw ftp_mod_times.txt(.gz) into a tidy DataFrame."""
    raw = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=["assembly_dir", "modification_date"],
        dtype=str,
    )

    ids_and_names = raw["assembly_dir"].apply(parse_assembly_dir)
    raw["assembly_id"] = ids_and_names.apply(lambda x: x[0])
    raw["assembly_name"] = ids_and_names.apply(lambda x: x[1])
    raw["source"] = raw["assembly_id"].apply(infer_source)
    raw["modification_date"] = pd.to_datetime(raw["modification_date"], format="%Y/%m/%d").dt.date

    return raw[["assembly_dir", "assembly_id", "assembly_name", "source", "modification_date"]]


def load_interesting_assembly_ids(path: Path) -> set[str]:
    """Read in the text file containing the IDs of the assemblies of interest."""
    df = pd.read_csv(path, sep="\t", header=None, names=["assembly_id"], dtype=str)
    return set(df["assembly_id"])


def build_interesting_membership_mask(
    parsed_df: pd.DataFrame,
    interesting_ids: set[str],
    filter_refseq: bool,  # noqa: FBT001
) -> pd.Series:
    """Decide which rows of parsed_df should be included in the filtered subset.

    - GenBank rows are always required to be in interesting_ids.
    - RefSeq rows are required to be in interesting_ids only if filter_refseq is True;
      otherwise all RefSeq rows pass through unconditionally.
    """
    in_interesting = parsed_df["assembly_id"].isin(interesting_ids)

    if filter_refseq:
        return in_interesting

    is_refseq = parsed_df["source"] == Source.REFSEQ
    return is_refseq | in_interesting


def filter_to_interesting(
    parsed_df: pd.DataFrame,
    interesting_ids: set[str],
    filter_refseq: bool = False,  # noqa: FBT001, FBT002
) -> pd.DataFrame:
    """Keep interesting (and, optionally, all RefSeq) assemblies.

    Assemblies are deduplicated on the most recent modification date per assembly_id.
    """
    mask = build_interesting_membership_mask(parsed_df, interesting_ids, filter_refseq)

    subset = parsed_df[mask].copy()
    sorted_subset = subset.sort_values("modification_date", ascending=False)
    deduped_subset = sorted_subset.drop_duplicates(subset="assembly_id", keep="first")
    return deduped_subset.sort_values("assembly_id")


def split_by_period(interesting_df: pd.DataFrame, mod_date: date) -> dict[Period, pd.DataFrame]:
    """Split into rows on/after mod_date ("after") vs strictly before it."""
    is_after = interesting_df["modification_date"] >= mod_date
    return {
        Period.AFTER: interesting_df[is_after].sort_values("assembly_id"),
        Period.BEFORE: interesting_df[~is_after].sort_values("assembly_id"),
    }


def split_by_source(df: pd.DataFrame) -> dict[Source, pd.DataFrame]:
    """Split dataframes according to the contents of the source column."""
    return {source: df[df["source"] == source] for source in Source}


def build_splits(interesting_df: pd.DataFrame, mod_date: date) -> Splits:
    """Produce the four (source, period) DataFrames used downstream."""
    by_period = split_by_period(interesting_df, mod_date)
    splits: Splits = {}
    for period, period_df in by_period.items():
        for source, source_df in split_by_source(period_df).items():
            splits[(source, period)] = source_df
    return splits


def save_tsv(df: pd.DataFrame, path: Path) -> None:
    """Save output as TSV."""
    df.to_csv(path, sep="\t", index=False)
    print(f"  wrote {len(df)} rows -> {path}")


def save_split_tables(splits: Splits, config: Config) -> None:
    """Save the split data files as TSV."""
    for (source, period), df in splits.items():
        path = config.derived_dir / f"{source.value}_{period.value}_{config.mod_date_tag}.tsv"
        save_tsv(df, path)


def batch(items: list[str], size: int) -> Iterable[list[str]]:
    """Batch up items into batches of size `size`."""
    for i in range(0, len(items), size):
        yield items[i : i + size]


def batch_count(n_items: int, batch_size: int) -> int:
    """Count 'em batches and weep!"""
    return math.ceil(n_items / batch_size) if n_items else 0


def compute_batch_width(splits: Splits, batch_size: int) -> int:
    """Zero-padding width wide enough for the largest batch count across all splits."""
    max_batches = max((batch_count(len(df), batch_size) for df in splits.values()), default=0)
    return max(len(str(max_batches)), 1)


def write_id_batches(ids: list[str], out_dir: Path, prefix: str, batch_size: int, width: int) -> int:
    """Spew out chunks of ID lists."""
    n_written = 0
    for i, chunk in enumerate(batch(ids, batch_size), start=1):
        file_num = str(i).zfill(width)
        out_path = out_dir / f"{prefix}-{file_num}.txt"
        out_path.write_text("\n".join(chunk) + "\n")
        n_written += 1
    return n_written


def write_batched_ids(splits: Splits, config: Config) -> None:
    """Write out the batched IDs into sequentially numbered files."""
    width = compute_batch_width(splits, config.batch_size)
    for (source, period), df in splits.items():
        ids = df["assembly_id"].tolist()
        prefix = f"{period.value}_{config.mod_date_tag}"
        n_written = write_id_batches(ids, config.source_dir(source), prefix, config.batch_size, width)
        print(f"  wrote {n_written} batch file(s) for {source.value}/{period.value}")


def run_pipeline(config: Config) -> None:
    """Run the ID list generation pipeline."""
    config.ensure_dirs()

    print(f"Parsing {config.ftp_mod_times_path} ...")
    parsed_df = load_ftp_mod_times(config.ftp_mod_times_path)
    save_tsv(parsed_df, config.derived_dir / "ftp_mod_times_parsed.tsv")

    print(f"Reading {config.assemblies_of_interest_path} ...")
    interesting_ids = load_interesting_assembly_ids(config.assemblies_of_interest_path)
    print(f"  found {len(interesting_ids)} unique assembly IDs")

    mode = "filtering RefSeq and GenBank" if config.filter_refseq else "keeping all RefSeq, filtering GenBank only"
    print(f"  NCBI filtering mode: {mode}")

    interesting_df = filter_to_interesting(parsed_df, interesting_ids, filter_refseq=config.filter_refseq)
    save_tsv(interesting_df, config.derived_dir / "interesting_assemblies.tsv")

    splits = build_splits(interesting_df, config.mod_date)
    save_split_tables(splits, config)
    write_batched_ids(splits, config)

    print("Done!")


if __name__ == "__main__":
    CliApp.run(Config, cli_args=sys.argv[1:]).cli_cmd()
