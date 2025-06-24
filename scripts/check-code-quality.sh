#!/bin/bash
# Script to run all code quality checks and fixes
# This script is designed to be used by GitHub Copilot Agent and developers

set -e

echo "🔍 Running code quality checks and fixes..."

# Check if poetry is available
if ! command -v poetry &> /dev/null; then
    echo "❌ Poetry not found. Please install poetry first."
    exit 1
fi

# Install dependencies if needed
echo "📦 Ensuring dependencies are installed..."
poetry install --with dev

# Run pre-commit hooks (if available and working)
if [ -f ".pre-commit-config.yaml" ]; then
    echo "🪝 Running pre-commit hooks..."
    if poetry run pre-commit run --all-files; then
        echo "✅ Pre-commit hooks passed"
    else
        echo "⚠️  Pre-commit hooks found issues or failed (will fix with individual tools)"
    fi
fi

# Format Python code with Black
echo "🖤 Formatting Python code with Black..."
poetry run black .

# Sort imports with isort (using poetry run to ensure it's in the venv)
echo "📦 Sorting imports with isort..."
if poetry run python -c "import isort" 2>/dev/null; then
    poetry run python -m isort .
else
    echo "⚠️  isort not available, skipping import sorting"
fi

# Fix linting issues with Ruff
echo "🦀 Fixing linting issues with Ruff..."
poetry run ruff check --fix .

# Run final checks
echo "✅ Running final checks..."

# Check Black formatting
echo "  - Checking Black formatting..."
poetry run black --check .

# Check Ruff linting
echo "  - Checking Ruff linting..."
poetry run ruff check .

# Check type hints with pyright
echo "  - Checking type hints with pyright..."
if poetry run python -c "import pyright" 2>/dev/null; then
    poetry run pyright
else
    echo "⚠️  Pyright not available, skipping type checking"
fi

# Check documentation coverage with interrogate
echo "  - Checking documentation coverage..."
if poetry run python -c "import interrogate" 2>/dev/null; then
    poetry run interrogate src
else
    echo "⚠️  Interrogate not available, skipping documentation check"
fi

echo "🎉 Core code quality checks completed!"
echo "💡 You can now safely commit your changes."