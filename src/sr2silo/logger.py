"""Logging functions for sr2silo."""

from __future__ import annotations

import logging
from contextlib import contextmanager


@contextmanager
def suppress_info_and_below():
    """Suppress INFO and below log messages.

    A context manager that temporarily changes the logging level to WARNING,
    suppressing INFO, DEBUG, and NOTSET level messages within its context.
    The original logging level is restored when exiting the context.
    """
    logger = logging.getLogger()
    original_level = logger.getEffectiveLevel()  # Save current level
    logger.setLevel(logging.WARNING)  # Suppress INFO and below
    try:
        yield
    finally:
        logger.setLevel(original_level)  # Restore original level
