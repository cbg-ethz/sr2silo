"""
This module contains tests for the conversion functions in the sr2silo package.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from sr2silo.process import bam_to_sam
from sr2silo.process.convert import normalize_reads


def test_bam_to_sam(bam_data):
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



def test_sort_bam_file():
    """Test the sort_bam_file function."""

    raise NotImplementedError


def test_create_index():
    """Test the index_bam_file function."""

    raise NotImplementedError


def test_bam_to_fasta():
    """Test the bam_to_fasta function."""

    raise NotImplementedError

def test_pad_alignment():
    """Test the pad_alignment function."""

    raise NotImplementedError

def test_sam_to_seq_and_indels():
    """Test the sam_to_seq_and_indels function."""

    raise NotImplementedError

def test_get_gene_set_from_reference():
    """Test the get_gene_set_from_reference function."""

    raise NotImplementedError