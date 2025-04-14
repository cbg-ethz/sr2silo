#   -------------------------------------------------------------
#   Copyright (c) Microsoft Corporation. All rights reserved.
#   Licensed under the MIT License. See LICENSE in project root for information.
#   -------------------------------------------------------------
"""sr2silo connects pairs, normalizes reads, and converts BAM to SAM files."""
from __future__ import annotations

import importlib.metadata

import sr2silo.vpipe as vpipe

try:
    __version__ = importlib.metadata.version("sr2silo")
except importlib.metadata.PackageNotFoundError:
    # Package is not installed
    import tomllib  # Use tomllib for reading pyproject.toml (Python 3.11+)

    try:
        with open("pyproject.toml", "rb") as f:
            pyproject_data = tomllib.load(f)
            __version__ = pyproject_data["tool"]["poetry"]["version"]
    except (FileNotFoundError, KeyError):
        __version__ = "0.1.0"  # Fallback if pyproject.toml is missing or invalid

__all__ = ["vpipe"]
