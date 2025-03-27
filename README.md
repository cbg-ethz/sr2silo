# sr2silo
## Wrangele BAM nucleotide alignments to cleartext alignments
<picture>
  <source
    media="(prefers-color-scheme: light)"
    srcset="resources/graphics/logo.svg">
  <source
    media="(prefers-color-scheme: dark)"
    srcset="resources/graphics/logo_dark_mode.svg">
  <img alt="Logo" src="resources/logo.svg" width="15%" />
</picture>

[![Project Status: POC – This project is currently under active development.](https://www.repostatus.org/badges/latest/concept.svg)](https://www.repostatus.org/#concept)
[![CI/CD](https://github.com/gordonkoehn/UsefulGnom/actions/workflows/test.yml/badge.svg)](https://github.com/gordonkoehn/UsefulGnom/actions/workflows/test.yml)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC.svg)](https://docs.pytest.org/en/stable/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)
[![Pyright](https://img.shields.io/badge/type%20checked-pyright-blue.svg)](https://github.com/microsoft/pyright)

### General Use: Convert Nucleotide Alignment Reads - CIGAR in .BAM to Cleartext JSON
sr2silo can convert millions of Short-Read nucleotide read in the form of a .bam CIGAR
alignments to cleartext alignments. Further, it will gracefully extract insertions
and deletions. Optionally, sr2silo can translate and align each read using [diamond / blastX](https://github.com/bbuchfink/diamond). And again handle insertions and deletions.

Your input `.bam/.sam` with one line as:
````
294	163	NC_045512.2	79	60	31S220M	=	197	400	CTCTTGTAGAT	FGGGHHHHLMM	...
````

sr2silo outputs per read a JSON (mock output):

```
{
  "metadata":{
    "read_id":"AV233803:AV044:2411515907:1:10805:5199:3294",
      ...
    },
    "nucleotideInsertions":{
                            "main":[10 : ACTG]
                            },
    "aminoAcidInsertions":{
                            "E":[],
                            ...
                            "ORF1a":[2323 : TG, 2389 : CA],
                            ...
                            "S":[23 : A]
                            },
    "alignedNucleotideSequences":
                                {
                                  "main":"NNNNNNNNNNNNNNNNNNCGGTTTCGTCCGTGTTGCAGCCG...GTGTCAACATCTTAAAGATGGCACTTGTGNNNNNNNNNNNNNNNNNNNNNNNN"
                                  },
    "unalignedNucleotideSequences":{
                                  "main":"CGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTTGTCCGGGTGTGA...TACAGGTTCGCGACGTGCTCGTGTGAAAGATGGCACTTGTG"
                                  },
    "alignedAminoAcidSequences":{
                "E":"",
                ...
                "ORF1a":"...NMESLVPGFNEKTHVQLSLPVLQVRVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN...",
                ...
                "S":""}
      }
```

The total output is handled in an `.ndjson.zst`.

### Resource Requirements

When running sr2silo, particularly the `import-to-loculus` command, be aware of memory and storage requirements:

- Standard configuration uses 8GB RAM and one CPU core
- Processing batches of 100k reads requires ~3GB RAM plus ~3GB for Diamond
- Temporary storage needs (especially on clusters) can reach 30-50GB

For detailed information about resource requirements, especially for cluster environments, please refer to the [Resource Requirements documentation](docs/usage/resource_requirements.md).

### Wrangling Short-Read Genomic Alignments for SILO Database

Originally this was started for wargeling short-read genomic alignments for from wastewater-sampling, into a format for easy import into [Loculus](https://github.com/loculus-project/loculus) and its sequence database SILO.

sr2silo is designed to process a nucliotide alignments from `.bam` files with metadata, translate and align reads in amino acids, gracefully handling all insertions and deletions and upload the results to the backend [LAPIS-SILO](https://github.com/GenSpectrum/LAPIS-SILO).

For the V-Pipe to Silo implementation we carry through the following metadata:
```
  "metadata":{
    "read_id":"AV233803:AV044:2411515907:1:10805:5199:3294",
    "sample_id":"A1_05_2024_10_08",
    "batch_id":"20241024_2411515907",
    "sampling_date":"2024-10-08",
    "sequencing_date":"2024-10-24",
    "location_name":"Lugano (TI)",
    "read_length":"250","primer_protocol":"v532",
    "location_code":"05",
    "flow_cell_serial_number":"2411515907"
    "sequencing_well_position":"A1",
    "primer_protocol_name":"SARS-CoV-2 ARTIC V5.3.2",
    "nextclade_reference":"sars-cov-2"
    }
```

### Setting up the repository

To build the package and maintain dependencies, we use [Poetry](https://python-poetry.org/).
In particular, it's good to install it and become familiar with its basic functionalities by reading the documentation.

### Installation

The project uses a modular environment system to separate core functionality, development requirements, and workflow dependencies. Environment files are located in the `environments/` directory:

#### Core Environment Setup

For basic usage of sr2silo:
```bash
make setup
```
This creates the core conda environment with essential dependencies and installs the package using Poetry.

#### Development Environment

For development work:
```bash
make setup-dev
```
This command sets up the development environment with Poetry.
#### Workflow Environment

For working with the snakemake workflow:
```bash
make setup-workflow
```
This creates an environment specifically configured for running the sr2silo in snakemake workflows.

#### All Environments

You can set up all environments at once:
```bash
make setup-all
```

### Additional Setup for Development

After setting up the development environment:
```bash
conda activate sr2silo-dev
poetry install --with dev
poetry run pre-commit install
```

### Run Tests

```bash
make test
```
or
```bash
conda activate sr2silo-dev
pytest
```

### Run CLI

The sr2silo CLI has two main commands:

1. `run` - Not yet implemented command for future functionality
2. `import-to-loculus` - Convert BAM alignments to SILO format and optionally upload

#### Basic Usage

The main command you'll use is `import-to-loculus`:

```bash
sr2silo import-to-loculus \
    --input-file INPUT.bam \
    --sample-id SAMPLE_ID \
    --batch-id BATCH_ID \
    --timeline-file TIMELINE.tsv \
    --primer-file PRIMERS.yaml \
    --output-fp OUTPUT.ndjson \
    --reference sars-cov-2
```

#### Required Arguments

- `--input-file, -i`: Path to the input BAM alignment file
- `--sample-id, -s`: Sample ID to use for metadata
- `--batch-id, -b`: Batch ID to use for metadata
- `--timeline-file, -t`: Path to the timeline metadata file
- `--primer-file, -p`: Path to the primers configuration file
- `--output-fp, -o`: Path for the output file (will be auto-suffixed with .ndjson.zst)

#### Optional Arguments

- `--reference, -r`: Reference genome to use (default: "sars-cov-2")
- `--upload/--no-upload`: Whether to upload results to S3 and submit to SILO (default: no-upload)

#### Example Usage

Here's a complete example with sample data:

```bash
sr2silo import-to-loculus \
    --input-file ./data/sample/alignments/REF_aln_trim.bam \
    --sample-id "A1_05_2024_10_08" \
    --batch-id "20241024_2411515907" \
    --timeline-file ./data/timeline.tsv \
    --primer-file ./data/primers.yaml \
    --output-fp ./results/output.ndjson \
    --reference sars-cov-2
```

To also upload the results to SILO, add the `--upload` flag:

```bash
sr2silo import-to-loculus \
    # ...same arguments as above... \
    --upload
```

### Tool Sections
The code quality checks run on GitHub can be seen in
 - ``.github/workflows/test.yml`` for the python package CI/CD,

We are using:

  * [Ruff](https://github.com/charliermarsh/ruff) to lint the code.
  * [Black](https://github.com/psf/black) to format the code.
  * [Pyright](https://github.com/microsoft/pyright) to check the types.
  * [Pytest](https://docs.pytest.org/) to run the unit tests code and workflows.
  * [Interrogate](https://interrogate.readthedocs.io/) to check the documentation.


## Contributing

This project welcomes contributions and suggestions. For details, visit the repository's [Contributor License Agreement (CLA)](https://cla.opensource.microsoft.com) and [Code of Conduct](https://opensource.microsoft.com/codeofconduct/) pages.
