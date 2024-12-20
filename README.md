# sr2silo
<picture>
  <source
    media="(prefers-color-scheme: light)"
    srcset="resources/logo.svg">
  <source
    media="(prefers-color-scheme: dark)"
    srcset="resources/logo_dark_mode.svg">
  <img alt="Logo" src="resources/logo.svg" width="15%" />
</picture>

[![Project Status: WIP – This project is currently under active development.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![CI/CD](https://github.com/gordonkoehn/UsefulGnom/actions/workflows/test.yml/badge.svg)](https://github.com/gordonkoehn/UsefulGnom/actions/workflows/test.yml)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC.svg)](https://docs.pytest.org/en/stable/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)
[![Pyright](https://img.shields.io/badge/type%20checked-pyright-blue.svg)](https://github.com/microsoft/pyright)


### Wrangling Short-Read Genomic Alignments for SILO Database

This project will wrangle short-read genomic alignments, for example from wastewater-sampling, into a format for easy import into the SILO sequencing database.
### Usage of the V-Pipe Docker

The V-Pipe Docker is designed to process a single `.bam` file and upload the results to SILO.

## Project Organization

- `silo-input-transformer`: Is a rust based utility to handle the `fasta` to `ndjson` transformation and is here imported as a git submodule.
- `.github/workflows`: Contains GitHub Actions used for building, testing, and publishing.
install, and whether or not to mount the project directory into the container.
- `.vscode/settings.json`: Contains VSCode settings specific to the project, such as the Python interpreter to use and the maximum line length for auto-formatting.
- `src`: Place new source code here.
- `scripts`: Place new source code here, temporary and intermediate works.
- `tests`: Contains Python-based test cases to validate source code.
- `pyproject.toml`: Contains metadata about the project and configurations for additional tools used to format, lint, type-check, and analyze Python code.

### Setting up the repository

To build the package and maintain dependencies, we use [Poetry](https://python-poetry.org/).
In particular, it's good to install it and become familiar with its basic functionalities by reading the documentation.


### Setting up the Development Environment

1. Create and activate the conda environment from the `environment.yml` file:
  ```bash
  conda env create -f environment.yml
  conda activate sr2silo
  ```


2. Set up the environment with development tools:
  ```bash
  poetry install --with dev
  poetry run pre-commit install
  ```

Then, you will be able to run tests:
```bash
$ poetry run pytest
```
... or check the types:
```bash
$ poetry run pyright
```

Alternatively, you may prefer to work with the right Python environment using:
```bash
$ poetry shell
$ pytest
```

### [WIP]: Run V-Pipe to SILO Transformation
This is currently implemented as script and under heavy development.
To run, we recommend a build as a docker compose as it relies on other RUST components.

#### Configuration

Edit the `docker-compose.env` file in the `docker-compose` directory with the following paths:

```env
SAMPLE_DIR=../../../data/sr2silo/daemon_test/samples/A1_05_2024_10_08/20241024_2411515907/alignments/
SAMPLE_ID=A1_05_2024_10_08
BATCH_ID=20241024_2411515907
TIMELINE_FILE=../../../data/sr2silo/daemon_test/timeline.tsv
NEXTCLADE_REFERENCE=sars-cov2
RESULTS_DIR=./results
```


#### Docker Secrets
To upload the processed outputs S3 storage is required.

For sensitive information like AWS credentials, use Docker secrets. Create the following files in the secrets directory:

- `secrets/aws_access_key_id.txt`:

```YourAWSAccessKeyId```

- `secrets/aws_secret_access_key.txt`:

```YourAWSSecretAccessKey```

- `secrets/aws_default_region.txt`:
```YourAWSRegion```

#### Run Transformation

To process a single sample, run the following command:

```sh
docker-compose --env-file .env up --build
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
