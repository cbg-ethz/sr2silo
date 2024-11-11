"""
This module contains tests for the conversion functions in the sr2silo package.
"""

from __future__ import annotations

from pathlib import Path

from sr2silo.convert import bam_to_sam

# Define test data paths outside the function
INPUT_BAM_PATH = Path("tests/data/REF_aln_trim_subsample.bam")
EXPECTED_SAM_PATH = Path("tests/data/REF_aln_trim_subsample_expected.sam")


def test_bam_to_sam():
    """Test the bam_to_sam function."""
    sam_data = bam_to_sam(INPUT_BAM_PATH)

    with EXPECTED_SAM_PATH.open() as f:
        expected_data = f.read()

    # compare the sam_data with the expected data
    assert (
        sam_data == expected_data
    ), "The converted SAM data does not match the expected SAM data"
