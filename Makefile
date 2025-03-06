.PHONY: create-env
create-env:
	@if conda info --envs | awk '{print $$1}' | grep -wq sr2silo; then \
		echo "Conda environment sr2silo already exists. Updating environment..."; \
		conda env update -f environment.yml --prune; \
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

.PHONY: setup
setup: create-env activate-env install-poetry
	@echo "Environment setup complete."

.PHONY: clean
clean:
	@echo "Removing Conda environment..."
	@conda env remove -n sr2silo
	@echo "Environment removed."
