"""Implements Processing.

translation of alignments, read pairing, and normalization, and conversions.
"""

from __future__ import annotations

from sr2silo.process.convert import (
    bam_to_fasta_query,
    bam_to_sam,
    get_gene_set_from_ref,
    sam_to_bam,
    sort_and_index_bam,
    sort_bam_file,
    sort_sam_by_qname,
)
from sr2silo.process.interface import (
    AAInsertion,
    AAInsertionSet,
    AlignedRead,
    Gene,
    GeneSet,
    NucInsertion,
)
from sr2silo.process.merge import paired_end_read_merger
from sr2silo.process.translate_align import (
    curry_read_with_metadata,
    nuc_to_aa_alignment,
    parse_translate_align,
    parse_translate_align_in_batches,
)

__all__ = [
    # from sr2silo.process.convert
    "bam_to_fasta_query",
    "bam_to_sam",
    "get_gene_set_from_ref",
    "get_gene_set_from_ref",
    "sam_to_bam",
    "sort_and_index_bam",
    "sort_bam_file",
    "sort_sam_by_qname",
    # from sr2silo.process.interface
    "AAInsertion",
    "AAInsertionSet",
    "AlignedRead",
    "Gene",
    "GeneSet",
    "Gene",
    "GeneSet",
    "NucInsertion",
    # from sr2silo.process.merge
    "paired_end_read_merger",
    # from sr2silo.process.translate_align
    "nuc_to_aa_alignment",
    "parse_translate_align",
    "parse_translate_align_in_batches",
    "curry_read_with_metadata",
]
