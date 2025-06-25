#!/bin/bash
# Simple build script for local conda package

# Check if conda-build is installed
if ! conda list conda-build | grep -q conda-build; then
    echo "conda-build is not installed. Creating and using build environment..."
    # Check if build-env exists, if not create it
    if ! conda env list | grep -q build-env; then
        conda create -n build-env conda-build -c conda-forge -y
    fi
    echo "Activating build environment..."
    eval "$(conda shell.bash hook)"
    conda activate build-env
fi

# Build the package locally for testing
conda build conda-recipe -c conda-forge -c bioconda

# Get the path to the built package
PACKAGE_PATH=$(conda build conda-recipe --output)
echo "Package built at: $PACKAGE_PATH"

# Optional: Install locally for testing
echo "To install locally, run:"
echo "conda install -c local $PACKAGE_PATH"
