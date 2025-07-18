.PHONY: create-env
create-env:
	@if conda info --envs | awk '{print $$1}' | grep -wq sr2silo; then \
		echo "Conda environment sr2silo already exists. Updating environment..."; \
		conda env update -f environments/environment.yml --prune; \
	else \
		echo "Creating Conda environment..."; \
		conda env create -f environments/environment.yml; \
	fi

.PHONY: create-dev-env
create-dev-env:
	@if conda info --envs | awk '{print $$1}' | grep -wq sr2silo-dev; then \
		echo "Conda environment sr2silo-dev already exists. Updating environment..."; \
		conda env update -f environments/dev-environment.yml --prune; \
	else \
		echo "Creating development Conda environment..."; \
		conda env create -f environments/dev-environment.yml; \
	fi

.PHONY: create-workflow-env
create-workflow-env:
	@if conda info --envs | awk '{print $$1}' | grep -wq sr2silo-workflow; then \
		echo "Conda environment sr2silo-workflow already exists. Updating environment..."; \
		conda env update -f workflow/envs/sr2silo.yaml --prune; \
	else \
		echo "Creating workflow Conda environment..."; \
		conda env create -f workflow/envs/sr2silo.yaml; \
	fi

.PHONY: activate-env
activate-env:
	@echo "Activate the environment with: conda activate sr2silo"

.PHONY: activate-dev-env
activate-dev-env:
	@echo "Activate the development environment with: conda activate sr2silo-dev"

.PHONY: activate-workflow-env
activate-workflow-env:
	@echo "Activate the workflow environment with: conda activate sr2silo-workflow"

.PHONY: install
install:
	@echo "Installing the package in the sr2silo environment..."
	@conda run -n sr2silo pip install -e .

.PHONY: install-poetry-dev
install-poetry-dev:
	@echo "Installing Poetry with development dependencies in the sr2silo-dev environment..."
	@conda run -n sr2silo-dev pip install poetry
	@conda run -n sr2silo-dev poetry install --with dev

.PHONY: install-workflow
install-workflow:
	@echo "Installing the package in the sr2silo-workflow environment..."
	@conda run -n sr2silo-workflow pip install -e .

.PHONY: setup
setup: create-env activate-env install
	@echo "Environment setup complete."

.PHONY: setup-dev
setup-dev: create-dev-env activate-dev-env install-poetry-dev
	@echo "Development environment setup complete."

.PHONY: setup-workflow
setup-workflow: create-workflow-env activate-workflow-env install-workflow
	@echo "Workflow environment setup complete."

.PHONY: setup-all
setup-all: setup setup-dev setup-workflow
	@echo "All environments setup complete."

.PHONY: clean
clean:
	@echo "Removing Conda environment..."
	@conda env remove -n sr2silo
	@echo "Environment removed."

.PHONY: clean-dev
clean-dev:
	@echo "Removing development Conda environment..."
	@conda env remove -n sr2silo-dev
	@echo "Development environment removed."

.PHONY: clean-workflow
clean-workflow:
	@echo "Removing workflow Conda environment..."
	@conda env remove -n sr2silo-workflow
	@echo "Workflow environment removed."

.PHONY: clean-all
clean-all: clean clean-dev clean-workflow
	@echo "All environments removed."

.PHONY: test
test:
	@echo "Running tests in the development environment..."
	@conda run -n sr2silo-dev pytest
