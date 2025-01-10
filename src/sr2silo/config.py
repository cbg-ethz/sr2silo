"""Configuration utilities."""
from __future__ import annotations

import os


def is_ci_environment() -> bool:
    """Check if running in a CI environment."""
    return os.getenv("CI", "false").lower() in ("yes", "true", "t", "1")
