"""Implements Processing.

  translation of alignments, read pairing, and normalization, and conversions.
  """

from __future__ import annotations

from sr2silo.process.convert import (
    bam_to_sam,
    get_gene_set_from_ref,
    pad_alignment,
    sort_and_index_bam,
)
from sr2silo.process.interface import (
    AAInsertion,
    AAInsertionSet,
    AlignedRead,
    Gene,
    NucInsertion,
)

__all__ = [
    "bam_to_sam",
    "pad_alignment",
    "AAInsertion",
    "AAInsertionSet",
    "AlignedRead",
    "NucInsertion",
    "Gene",
    "sort_and_index_bam",
    "nuc_to_aa_alignment",
    "GeneSet",
    "get_gene_set_from_ref",
]
