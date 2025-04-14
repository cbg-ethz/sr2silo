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
    __version__ = "0.0.4"  # Fallback to match pyproject.toml

__all__ = ["vpipe"]
