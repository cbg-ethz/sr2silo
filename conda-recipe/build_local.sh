#!/bin/bash
# Simple build script for local conda package

# Build the package locally for testing
conda build conda-recipe -c conda-forge -c bioconda

# Get the path to the built package
PACKAGE_PATH=$(conda build conda-recipe --output)
echo "Package built at: $PACKAGE_PATH"

# Optional: Install locally for testing
echo "To install locally, run:"
echo "conda install -c local $PACKAGE_PATH"
