from pathlib import Path

import click
from refseq_pipeline.core.spark_delta import build_spark


@click.command()
@click.option("--database", required=True, help="Name of Spark SQL database.")
@click.option("--table", required=True, help="Name of the Delta table to register.")
@click.option("--path", required=True, type=click.Path(exists=False), help="Path to the Delta table.")
def main(database, table, path):
    abs_path = Path(path).expanduser().resolve()
    print(f"[debug] Received path: {abs_path}")
    print(f"[debug] Exists on disk? {abs_path.exists()}")

    if not abs_path.exists():
        raise click.ClickException(f"[error] Path does not exist: {abs_path}")

    spark = build_spark(database)
    spark.sql(f"CREATE TABLE IF NOT EXISTS {database}.{table} USING DELTA LOCATION '{abs_path}'")
    print(f"[success] Table registered: {database}.{table} → {abs_path}")


if __name__ == "__main__":
    main()


"""
python -m refseq_pipeline.cli.register_table \
  --database refseq_api \
  --table assembly_stats \
  --path "/Users/yuewang/Documents/RefSeq/delta_data/refseq_api/assembly_stats"

"""
