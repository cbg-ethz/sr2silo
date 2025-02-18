"""
This module contains tests for the conversion functions in the sr2silo package.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict

import pytest

from sr2silo.process import bam_to_sam
from sr2silo.process.convert import (
    normalize_reads,
    pad_alignment,
    sam_to_seq_and_indels,
)
from sr2silo.process.interface import AAInsertion


def test_bam_to_sam(bam_data: Dict):
    """Test the bam_to_sam function."""

    print(bam_data)

    sam_data = bam_to_sam(bam_data["bam_path"])

    # compare the sam_data with the expected data
    assert (
        sam_data == bam_data["sam_data"]
    ), "The converted SAM data does not match the expected SAM data"


def test_normalize_reads(sam_with_insert_data):
    """Test the normalize_reads function normalize_reads"""

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_dir_path = Path(temp_dir)
        output_fasta = temp_dir_path / "REF_aln_trim_subsample.fasta"
        output_insertions = temp_dir_path / "REF_aln_trim_subsample_insertions.fasta"

        # call the function
        normalize_reads(
            sam_with_insert_data["sam_data"], output_fasta, output_insertions
        )

        # check that the output files were created
        assert output_fasta.exists(), "The output fasta file was not created"
        assert output_insertions.exists(), "The output insertions file was not created"

        # check that the output files are not empty
        assert output_fasta.stat().st_size > 0, "The output fasta file is empty"
        assert (
            output_insertions.stat().st_size > 0
        ), "The output insertions file is empty"

        print(f"Output fasta file: {output_fasta}")
        print(f"Output insertions file: {output_insertions}")

        # check that the output files match the expected files

        with output_fasta.open() as f:
            output_fasta_str = f.read()

        assert (
            output_fasta_str == sam_with_insert_data["cleartext"]
        ), "The output fasta data does not match the expected fasta data"

        with output_insertions.open() as f:
            output_insertions_str = f.read()

        assert (
            output_insertions_str == sam_with_insert_data["insertions"]
        ), "The output insertions data does not match the expected insertions data"


@pytest.mark.skip(reason="Not implemented")
def test_sort_bam_file():
    """Test the sort_bam_file function."""

    raise NotImplementedError


@pytest.mark.skip(reason="Not implemented")
def test_create_index():
    """Test the index_bam_file function."""

    raise NotImplementedError


@pytest.mark.skip(reason="Not implemented")
def test_bam_to_fasta():
    """Test the bam_to_fasta function."""

    raise NotImplementedError


def test_pad_alignment():
    """Test the pad_alignment function with various inputs."""

    # Test with string input and default unknown_char 'N'
    seq = "ACTG"
    ref_start = 2
    ref_length = 10  # Expected: "NNACTGNNNN" (2 left, 4 right)
    expected = "NNACTGNNNN"
    result = pad_alignment(seq, ref_start, ref_length)
    assert result == expected, f"Expected {expected}, got {result}"

    # Test with list input, same parameters
    seq_list = ["A", "C", "T", "G"]
    expected = "NNACTGNNNN"
    result = pad_alignment(seq_list, ref_start, ref_length)
    assert result == expected, f"Expected {expected}, got {result}"

    # Test with custom unknown_char '-' (variation of deletion char)
    expected = "--ACTG----"
    result = pad_alignment(seq, ref_start, ref_length, unknown_char="-")
    assert result == expected, f"Expected {expected}, got {result}"


def test_sam_to_seq_and_indels():
    """Test the sam_to_seq_and_indels function, including AAInsertion conversion."""

    # Example from the docstring:
    # sequence = "AGCTTAGCTAGCTT"
    # cigar = "5M1I5M1D3M"
    seq = "AGCTTAGCTAGCTT"
    cigar = "5M1I5M1D3M"

    # Expected:
    # - Cleartext: "AGCTTGCTAGCTT" [matches from 5M, 5M, 3M]
    # - Insertion: one insertion at position 5 with base "A"
    # - Deletion: one deletion at position 10 with length 1
    expected_cleartext = "AGCTTGCTAGCTT"
    expected_deletions = [(10, 1)]

    cleartext, insertions, deletions = sam_to_seq_and_indels(seq, cigar)

    assert (
        cleartext == expected_cleartext
    ), f"Expected cleartext {expected_cleartext}, got {cleartext}"
    assert (
        deletions == expected_deletions
    ), f"Expected deletions {expected_deletions}, got {deletions}"
    assert len(insertions) == 1, f"Expected 1 insertion, got {len(insertions)}"
    # Check that the insertion is an instance of AAInsertion and has correct attributes.
    ins = insertions[0]
    assert isinstance(ins, AAInsertion), "Insertion is not an instance of AAInsertion"
    assert ins.position == 5, f"Expected insertion position 5, got {ins.position}"
    assert ins.sequence == "A", f"Expected insertion sequence 'A', got {ins.sequence}"


def test_get_gene_set_from_reference():
    """Test the get_gene_set_from_reference function."""

    raise NotImplementedError
