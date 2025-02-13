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
    GeneSet,
)
from sr2silo.process.merge import pair_normalize_reads
from sr2silo.process.translate_align import (
                  translate_nextclade,
                  nuc_to_aa_alignment
)
from sr2silo.process.convert import get_gene_set_from_ref

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
    "nuc_to_aa_alignment"
    "GeneSet",
    "get_gene_set_from_ref"
]
