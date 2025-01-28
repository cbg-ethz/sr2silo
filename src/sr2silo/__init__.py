#   -------------------------------------------------------------
#   Copyright (c) Microsoft Corporation. All rights reserved.
#   Licensed under the MIT License. See LICENSE in project root for information.
#   -------------------------------------------------------------
"""sr2silo connects pairs, normalizes reads, and converts BAM to SAM files."""
from __future__ import annotations

import sr2silo.process as process
import sr2silo.vpipe as vpipe

__version__ = "0.0.2"

__all__ = ["vpipe", "process"]
