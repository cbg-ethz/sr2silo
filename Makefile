.PHONY: create-env
create-env:
	@if conda info --envs | awk '{print $$1}' | grep -wq sr2silo; then \
	    echo "Conda environment sr2silo already exists. Skipping creation."; \
	else \
	    echo "Creating Conda environment..."; \
	    conda env create -f environment.yml; \
	fi

.PHONY: activate-env
activate-env:
	@echo "Activate the environment with: conda activate sr2silo"

.PHONY: install-poetry
install-poetry:
	@echo "Installing Poetry in the sr2silo environment..."
	@conda run -n sr2silo pip install poetry
	@conda run -n sr2silo poetry install

.PHONY: install-snakemake
install-snakemake:
	@echo "Installing Snakemake in the sr2silo environment..."
	@conda run -n sr2silo pip install snakemake

.PHONY: setup
setup: create-env activate-env install-poetry install-diamond
	@echo "Environment setup complete."

.PHONY: setup-with-snakemake
setup-with-snakemake: create-env activate-env install-poetry install-diamond install-snakemake
	@echo "Environment setup with Snakemake complete."

.PHONY: install-diamond
install-diamond:
	@echo "Installing Diamond into the sr2silo environment..."
	conda install -n sr2silo -c bioconda -c conda-forge diamond

.PHONY: clean
clean:
	@echo "Removing Conda environment..."
	conda env remove -n sr2silo
	@echo "Environment removed."
