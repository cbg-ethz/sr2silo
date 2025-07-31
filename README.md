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
[![CI/CD](https://github.com/cbg-ethz/sr2silo/actions/workflows/test.yml/badge.svg)](https://github.com/cbg-ethz/sr2silo/actions/workflows/test.yml)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC.svg)](https://docs.pytest.org/en/stable/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)
[![Pyright](https://img.shields.io/badge/type%20checked-pyright-blue.svg)](https://github.com/microsoft/pyright)

### General Use: Convert Nucleotide Alignment Reads - CIGAR in .BAM to Cleartext JSON
sr2silo can convert millions of Short-Read nucleotide reads in the form of .bam CIGAR
alignments to cleartext alignments compatible with LAPIS-SILO v0.8.0+. It gracefully extracts insertions
and deletions. Optionally, sr2silo can translate and align each read using [diamond / blastX](https://github.com/bbuchfink/diamond), handling insertions and deletions in amino acid sequences as well.

Your input `.bam/.sam` with one line as:
````
294	163	NC_045512.2	79	60	31S220M	=	197	400	CTCTTGTAGAT	FGGGHHHHLMM	...
````

sr2silo outputs per read a JSON (compatible with LAPIS-SILO v0.8.0+):

```json
{
  "read_id": "AV233803:AV044:2411515907:1:10805:5199:3294",
  "sample_id": "A1_05_2024_10_08",
  "batch_id": "20241024_2411515907",
  "sampling_date": "2024-10-08",
  "location_name": "Lugano (TI)",
  "read_length": "250",
  "location_code": "05",
  "main": {
    "sequence": "CGGTTTCGTCCGTGTTGCAGCCG...GTGTCAACATCTTAAAGATGGCACTTGTG",
    "insertions": ["10:ACTG", "456:TACG"],
    "offset": 4545
  },
  "unaligned_main": "CGGTTTCGTCCGTGTTGCAGCCGATCATCTAGGT...TACAGGTTCGCGACGTGCTCGTGTGAAAGATGGCACTTGTG",
  "S": {
    "sequence": "MESLVPGFNEKTHVQLSLPVLQVRVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGV",
    "insertions": ["23:A", "145:KLM"],
    "offset": 78
  },
  "ORF1a": {
    "sequence": "XXXMESLVPGFNEKTHVQLSLPVLQVRVRGFGDSVEEVLSEARQHLKDGTCGLV",
    "insertions": ["2323:TG", "2389:CA"],
    "offset": 678
  },
  "E": null,
  "M": null,
  "N": null,
  "ORF1b": null,
  "ORF3a": null,
  "ORF6": null,
  "ORF7a": null,
  "ORF7b": null,
  "ORF8": null,
  "ORF10": null
}
```

The total output is handled in an `.ndjson.zst`.

### Resource Requirements

When running sr2silo, particularly the `process-from-vpipe` command, be aware of memory and storage requirements:

- Standard configuration uses 8GB RAM and one CPU core
- Processing batches of 100k reads requires ~3GB RAM plus ~3GB for Diamond
- Temporary storage needs (especially on clusters) can reach 30-50GB

For detailed information about resource requirements, especially for cluster environments, please refer to the [Resource Requirements documentation](docs/usage/resource_requirements.md).

### Wrangling Short-Read Genomic Alignments for SILO Database

Originally this was started for wrangling short-read genomic alignments from wastewater-sampling, into a format for easy import into [Loculus](https://github.com/loculus-project/loculus) and its sequence database SILO.

sr2silo is designed to process nucleotide alignments from `.bam` files with metadata, translate and align reads in amino acids, gracefully handling all insertions and deletions and upload the results to the backend [LAPIS-SILO](https://github.com/GenSpectrum/LAPIS-SILO) v0.8.0+.

**New Output Format for LAPIS-SILO v0.8.0+:**
- Metadata fields are now at the root level (no nested "metadata" object)
- Genomic segments use a structured format with `sequence`, `insertions`, and `offset` fields
- The main nucleotide segment is required and contains the primary alignment
- Gene segments (S, ORF1a, etc.) contain amino acid sequences or `null` if empty
- Insertions use the format `"position:sequence"` (e.g., `"123:ACGT"`)
- Unaligned sequences are prefixed with `unaligned_` (e.g., `unaligned_main`)

For the V-Pipe to Silo implementation we include the following metadata fields at the root level:
```json
{
  "read_id": "AV233803:AV044:2411515907:1:10805:5199:3294",
  "sample_id": "A1_05_2024_10_08",
  "batch_id": "20241024_2411515907",
  "sampling_date": "2024-10-08",
  "location_name": "Lugano (TI)",
  "read_length": "250",
  "location_code": "05"
}
```

### Setting up the repository

To build the package and maintain dependencies, we use [Poetry](https://python-poetry.org/).
In particular, it's good to install it and become familiar with its basic functionalities by reading the documentation.

### Installation

sr2silo can be installed either from Bioconda or from source.

#### Install from Bioconda

The easiest way to install sr2silo is through the Bioconda channel:

```bash
# Add necessary channels if you haven't already
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Install sr2silo
conda install sr2silo
```

#### Install from Source

For development purposes or to install the latest version, you can install from source using Poetry:

The project uses a modular environment system to separate core functionality, development requirements, and workflow dependencies. Environment files are located in the `environments/` directory:

##### Core Environment Setup

For basic usage of sr2silo:
```bash
make setup
```
This creates the core conda environment with essential dependencies and installs the package using Poetry.

##### Development Environment

For development work:
```bash
make setup-dev
```
This command sets up the development environment with Poetry.
##### Workflow Environment

For working with the snakemake workflow:
```bash
make setup-workflow
```
This creates an environment specifically configured for running the sr2silo in snakemake workflows.

##### All Environments

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

### Usage

sr2silo follows a two-step workflow:

1. **Process data:** `sr2silo process-from-vpipe --help`
2. **Submit to Loculus:** `sr2silo submit-to-loculus --help`

```bash
# Example: Process V-Pipe data
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --lapis-url https://lapis.cov-spectrum.org/open/v2 \
    --output-fp output.ndjson

# Example: Submit to Loculus (use environment variables for credentials)
export KEYCLOAK_TOKEN_URL=https://auth.example.com/token
export SUBMISSION_URL=https://api.example.com/submit
export GROUP_ID=123
export USERNAME=your-username
export PASSWORD=your-password

sr2silo submit-to-loculus --processed-file output.ndjson.zst
```

**Note:** Use environment variables for credentials to avoid exposing sensitive information in command history.

### Environment Variable Configuration

sr2silo supports flexible configuration through environment variables, making it easy to use in different deployment scenarios including conda packages and pip installations.

**Key features:**
- CLI parameters override environment variables
- **Recommended for credentials to avoid exposing sensitive information in command history**

**Common configuration via environment variables:**
```bash
# Authentication credentials (recommended approach for security)
export KEYCLOAK_TOKEN_URL=https://auth.example.com/token
export SUBMISSION_URL=https://backend.example.com/api
export GROUP_ID=123
export USERNAME=your-username
export PASSWORD=your-password

# Run with required CLI arguments (timeline file must be specified)
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file /path/to/timeline.tsv \
    --lapis-url https://lapis.cov-spectrum.org/open/v2 \
    --output-fp output.ndjson

# Submission using environment variables for credentials
sr2silo submit-to-loculus \
    --processed-file output.ndjson.zst
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
