# Contributing

We welcome contributions to sr2silo! Whether you're fixing bugs, adding features, improving documentation, or reporting issues, your help is appreciated.

## Ways to Contribute

- **Report bugs**: Open an issue describing the problem and how to reproduce it
- **Suggest features**: Open an issue to discuss new functionality
- **Submit pull requests**: Fix bugs or implement new features
- **Improve documentation**: Help make the docs clearer and more complete

## Development Setup

We use [Poetry](https://python-poetry.org/) for dependency management and packaging. The project provides multiple environment configurations in the `environments/` directory.

### Environment Options

**Core Environment** (basic usage):
```bash
make setup
```

**Development Environment** (recommended for contributors):
```bash
make setup-dev
```

**Workflow Environment** (for Snakemake workflows):
```bash
make setup-workflow
```

**All Environments**:
```bash
make setup-all
```

### After Setup

```bash
conda activate sr2silo-dev
poetry install --with dev
poetry run pre-commit install
```

### Running Tests

```bash
make test
# or directly:
pytest
```

## Getting Started

1. Fork the repository
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/sr2silo.git
   cd sr2silo
   ```
3. Set up the development environment (see above)
4. Create a feature branch following our [branching strategy](contributing/branching-strategy.md)
5. Make your changes and run tests
6. Submit a pull request to the `dev` branch

## Code Guidelines

- Follow existing code style and patterns
- Add tests for new functionality
- Update documentation as needed
- Keep commits focused and write descriptive commit messages

## Questions?

Open an issue on GitHub if you have questions about contributing.
