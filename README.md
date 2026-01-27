<div align="center">

<picture>
  <source media="(prefers-color-scheme: light)" srcset="resources/graphics/logo.svg">
  <source media="(prefers-color-scheme: dark)" srcset="resources/graphics/logo_dark_mode.svg">
  <img alt="sr2silo logo" src="resources/graphics/logo.svg" width="200px" />
</picture>

# sr2silo

**Convert BAM nucleotide alignments to cleartext alignments for LAPIS-SILO**

[![Status: Public Beta](https://img.shields.io/badge/Status-Public%20Beta-blue)](https://github.com/cbg-ethz/sr2silo)
[![CI/CD](https://github.com/cbg-ethz/sr2silo/actions/workflows/test.yml/badge.svg)](https://github.com/cbg-ethz/sr2silo/actions/workflows/test.yml)
[![Pytest](https://img.shields.io/badge/tested%20with-pytest-0A9EDC.svg)](https://docs.pytest.org/en/stable/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)
[![Pyright](https://img.shields.io/badge/type%20checked-pyright-blue.svg)](https://github.com/microsoft/pyright)

[Documentation](https://cbg-ethz.github.io/sr2silo/) · [Installation](#installation) · [Quick Start](#quick-start)

</div>

---

sr2silo processes short-read nucleotide alignments from `.bam` files, translates and aligns reads in amino acids, and outputs JSON compatible with [LAPIS-SILO](https://github.com/GenSpectrum/LAPIS-SILO) v0.8.0+.

## Installation

```bash
conda install -c bioconda sr2silo
```

## Quick Start

```bash
# Process BAM data
sr2silo process-from-vpipe \
    --input-file input.bam \
    --sample-id SAMPLE_001 \
    --timeline-file timeline.tsv \
    --organism covid \
    --output-fp output.ndjson.zst

# Submit to Loculus
sr2silo submit-to-loculus \
    --processed-file output.ndjson.zst
```

## Documentation

Full documentation is available at the [sr2silo documentation site](https://cbg-ethz.github.io/sr2silo/):

- [Configuration](https://cbg-ethz.github.io/sr2silo/usage/configuration/) - Environment variables and CLI options
- [Multi-Organism Support](https://cbg-ethz.github.io/sr2silo/usage/organisms/) - Supported organisms and adding new ones
- [Deployment](https://cbg-ethz.github.io/sr2silo/usage/deployment/) - Multi-virus cluster deployment
- [API Reference](https://cbg-ethz.github.io/sr2silo/api/loculus/) - Python API documentation

## Development

```bash
make setup-dev
conda activate sr2silo-dev
poetry install --with dev
pytest
```

## License

See [LICENSE](LICENSE) for details.
