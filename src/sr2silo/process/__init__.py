"""Implements Processing.

  translation of alignments, read pairing, and normalization, and conversions.
  """

from __future__ import annotations

from sr2silo.process.convert import bam_to_sam, pad_alignment, sort_and_index_bam
from sr2silo.process.interface import (
    AAInsertion,
    AAInsertionSet,
    AlignedRead,
    Gene,
    NucInsertion,
)
from sr2silo.process.merge import pair_normalize_reads
from sr2silo.process.translation_aligment import translate_nextclade

__all__ = [
    "bam_to_sam",
    "pair_normalize_reads",
    "translate_nextclade",
    "pad_alignment",
    "AAInsertion",
    "AAInsertionSet",
    "AlignedRead",
    "NucInsertion",
    "Gene",
    "sort_and_index_bam",
]
