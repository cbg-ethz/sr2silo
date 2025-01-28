# Makefile for managing the sr2silo installation

# Variables
ENV_NAME = sr2silo
ENV_FILE = environment.yml

# Default target
.PHONY: all
all: install

# Create the Conda environment
.PHONY: create_env
create_env:
	conda env create -f $(ENV_FILE)

# Build the Rust project as a Python extension
.PHONY: build_nextclade
build_nextclade:
	conda run -n $(ENV_NAME) bash -c "cd nextclade/packages/nextclade && maturin develop --release -b pyo3"

# Build the Rust project as a Python extension
.PHONY: build_silo_input_transformer
build_silo_input_transformer:
	conda run -n $(ENV_NAME) bash -c "cd silo_input_transformer && maturin develop --release"

# Install the sr2silo package
.PHONY: install
install: create_env sr2silo
	conda run -n $(ENV_NAME) pip install -e .

# Clean up the Conda environment
.PHONY: clean
clean:
	conda env remove -n $(ENV_NAME)
