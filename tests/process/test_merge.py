"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

from pathlib import Path

from sr2silo.process import paired_end_read_merger


def test_paired_end_read_merger():

    SAM_FP = Path("tests/data/bam/combined.sam")
    REF_GENOME_FP = Path("resources/sars-cov-2/nuc_reference_genomes.fasta")

    output_merged_sam_fp = Path("tests/data/bam/merged.sam")

    paired_end_read_merger(SAM_FP, REF_GENOME_FP, output_merged_sam_fp)

    with open(output_merged_sam_fp, "r") as f:
        lines = f.readlines()
