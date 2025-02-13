"""Implements all SILO and LAPIS specific actions, as
   data wrangling for SILO, establishing a connection to
   LAPIS and submitting sequences to LAPIS."""

from __future__ import annotations

from sr2silo.silo.lapis import LapisClient, Submission

__all__ = ["LapisClient", "Submission"]
