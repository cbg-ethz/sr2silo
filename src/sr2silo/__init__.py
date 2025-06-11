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
    __version__ = "unknown"

__all__ = ["vpipe"]
