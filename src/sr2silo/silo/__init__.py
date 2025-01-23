"""Implements all SILO and LAPIS specific actions, as
   data wrangling for SILO, establishing a connection to
   LAPIS and submitting sequences to LAPIS."""

from __future__ import annotations

from sr2silo.silo.lapis import LapisClient, Submission
from sr2silo.silo.wrangle import wrangle_for_transformer

__all__ = ["LapisClient", "Submission", "wrangle_for_transformer"]
