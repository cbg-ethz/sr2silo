"""Implements Processing.

translation of alignments, read pairing, and normalization, and conversions.
"""

from __future__ import annotations

from sr2silo.process.convert import (
    bam_to_fasta,
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
    GeneSet,
    NucInsertion,
)
from sr2silo.process.translate_align import (
    enrich_read_with_metadata,
    nuc_to_aa_alignment,
    parse_translate_align,
    parse_translate_align_in_batches,
)

__all__ = [
    "bam_to_sam",
    "bam_to_fasta",
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
    "parse_translate_align",
    "enrich_read_with_metadata",
    "parse_translate_align_in_batches",
]
