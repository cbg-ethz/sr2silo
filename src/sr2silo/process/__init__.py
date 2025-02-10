"""Implements Processing.

  translation of alignments, read pairing, and normalization, and conversions.
  """

from __future__ import annotations

from sr2silo.process.convert import (
    bam_to_sam,
    pad_alignment,
)
from sr2silo.process.merge import pair_normalize_reads
from sr2silo.process.translation_aligment import translate

__all__ = [
    "bam_to_sam",
    "pair_normalize_reads",
    "translate",
    "pad_alignment",
]
