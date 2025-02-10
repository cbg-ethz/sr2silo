.PHONY: create-env
create-env:
	@echo "Creating Conda environment..."
	conda env create -f environment.yml

.PHONY: activate-env
activate-env:
	activate-env:
	@echo "Activate the environment with: conda activate sr2silo"

.PHONY: install-poetry
install-poetry:
	@echo "Installing Poetry dependencies..."
	poetry install

.PHONY: setup
setup: create-env activate-env install-poetry install-diamond
	@echo "Environment setup complete."

.PHONY: install-diamond
install-diamond:
	@echo "Installing Diamond into the sr2silo environment..."
	conda install -n sr2silo -c bioconda -c conda-forge diamond

.PHONY: clean
clean:
	@echo "Removing Conda environment..."
	conda env remove -n sr2silo
	@echo "Environment removed."
