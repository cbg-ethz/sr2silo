"""Tests for the translation module."""

from __future__ import annotations

import logging
from pathlib import Path

import sr2silo.process.translate_align as translate_align
from sr2silo.process.translate_align import AlignedRead


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def test_translate():
    """Test the translation function."""

    translate_align.translate(
        [Path("tests/data/merged_expected.fasta")],
        Path("output/"),
        "nextstrain/sars-cov-2/XBB",
    )
    assert True


def test_parse_translate_align():
    """Test the parse_translate_align function.

      Current dataset misses Amino Acid Insertions - i.e. not tested here.
    """

    nuc_ref_fp = Path("resources/sars-cov-2/nuc_reference_genomes.fasta")
    aa_ref_fp = Path("resources/sars-cov-2/aa_reference_genomes.fasta")
    nuc_alignment_fp = Path("tests/data/bam/combined.bam")

    aligned_reads = translate_align.parse_translate_align(nuc_ref_fp, aa_ref_fp, nuc_alignment_fp)

    # load the expected aligned reads
    expected_aligned_reads = []
    with open("tests/data/process/aligned_reads.ndjson") as f:
        for line in f:
            # read in each line
            aligned_read = AlignedRead.from_str(line)

            expected_aligned_reads.append(aligned_read)

    assert aligned_reads == expected_aligned_reads


def test_read_in_AligendReads_nuc_seq():
    """Test the read_in_AlignedReads_nuc_seq function."""
    raise NotImplementedError


def test_read_in_AligendReads_nuc_ins():
    """Test the read_in_AlignedReads_nuc_ins function."""

    raise NotImplementedError


def test_read_in_AligendReads_aa_ins():
    """Test the read_in_AligendReads_aa_ins function."""
    raise NotImplementedError
