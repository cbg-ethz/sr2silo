[![Project Status: WIP – This project is currently under active development.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![CI/CD](https://github.com/gordonkoehn/UsefulGnom/actions/workflows/test.yml/badge.svg)](https://github.com/gordonkoehn/UsefulGnom/actions/workflows/test.yml)
[![Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC.svg)](https://docs.pytest.org/en/stable/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)
[![Pyright](https://img.shields.io/badge/type%20checked-pyright-blue.svg)](https://github.com/microsoft/pyright)

# sr2silo
### Wrangling Short-Read Genomic Alignments for SILO Database

This project will wrangle short-read genomic alignments, for example from wastewater-sampling, into a format for easy import into the SILO sequencing database.

## Project Organization

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

  Install Nextclade:
  ```bash
  conda install -c bioconda nextclade
  ```
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

#### Tool Sections
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
