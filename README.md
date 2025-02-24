# cdm-data-loader-utils
Repo for CDM input data loading and wrangling


## Installation

The data loader utils run on python 3.12 and above.

To install, run

```sh
> pip install -r requirements.txt
```

### Tests

To run the tests, execute the command:

```sh
> pytest tests/
```

To generate coverage for the tests, run

```sh
> pytest --cov=src --cov-report=xml tests/
```

The standard python `coverage` package is used and coverage can be generate as html or other formats by changing the parameters.


## Loading genomes, contigs, and features

The [genome loader](src/parsers/genome_loader.py) can be used to load and integrate data from related GFF and FASTA files. Currently, the loader requires a GFF file and two FASTA files (one for amino acid seqs, one for nucleic acid seqs) for each genome. The list of files to be processed should be specified in the genome paths file, which has the following format:

```json
{
    "FW305-3-2-15-C-TSA1.1": {
        "fna": "tests/data/FW305-3-2-15-C-TSA1/FW305-3-2-15-C-TSA1_scaffolds.fna",
        "gff": "tests/data/FW305-3-2-15-C-TSA1/FW305-3-2-15-C-TSA1_genes.gff",
        "protein": "tests/data/FW305-3-2-15-C-TSA1/FW305-3-2-15-C-TSA1_genes.faa"
    },
    "FW305-C-112.1": {
        "fna": "tests/data/FW305-C-112.1/FW305-C-112.1_scaffolds.fna",
        "gff": "tests/data/FW305-C-112.1/FW305-C-112.1_genes.gff",
        "protein": "tests/data/FW305-C-112.1/FW305-C-112.1_genes.faa"
    }
}
```

## Running bbmap stats and checkm2 on genome or contigset files

[run_tools.sh](scripts/run_tools.sh) runs the stats script from [bbmap](https://sourceforge.net/projects/bbmap/) and [checkm2](https://github.com/chklovski/CheckM2) on files with the suffix "fna". These tools can be installed using conda:

```sh
conda env create -f env.yml
conda activate genome_loader_env
# download the checkm2 database
checkm2 database --download
```

Run the stats and checkm2 tools with the following command:

```sh
bash scripts/run_tools.sh path/to/genome_paths_file.json output_dir
```
where `path/to/genome_paths_file.json` specifies the path to the genome paths file (format specified above) and `output_dir` is the directory for the results.
