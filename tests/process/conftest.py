"""
Fixtures for the process module.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict

import pytest

import sr2silo.process.translate_align as translate_align
from sr2silo.process.interface import AlignedRead


@pytest.fixture
def aligned_reads() -> Dict[str, AlignedRead]:
    """
    Small mock data with 42 real reads from the combined.bam file.

    Current dataset misses Amino Acid Insertions - i.e. not tested here.
    """

    nuc_ref_fp = Path("resources/sars-cov-2/nuc_reference_genomes.fasta")
    aa_ref_fp = Path("resources/sars-cov-2/aa_reference_genomes.fasta")
    nuc_alignment_fp = Path("tests/data/bam/combined.bam")

    aligned_reads = translate_align.parse_translate_align(
        nuc_ref_fp, aa_ref_fp, nuc_alignment_fp
    )
    return aligned_reads
