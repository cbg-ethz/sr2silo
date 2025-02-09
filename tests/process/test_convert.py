"""
This module contains tests for the conversion functions in the sr2silo package.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from sr2silo.process import bam_to_sam
from sr2silo.process.convert import bam_to_cleartext_alignment, normalize_reads


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


def test_bam_to_cleartext_alignment():
    """Test the bam_to_cleartext_alignment function."""

    ref_seq = Path("resources/sars-cov-2/NC_045512.2.fasta")

    expected_output_file = "tests/data/bam/expected_cleartext.ndjson"

    with open(expected_output_file) as f:
        expected_cleartext = f.read()

    save_to_another_dir = False  # Toggle variable to save to another directory

    with tempfile.TemporaryDirectory() as temp_dir:
        output_file = Path(temp_dir) / "REF_aln_trim_subsample_cleartext.fasta"

        bam_to_cleartext_alignment(
            Path("tests/data/bam/combined.bam"), output_file, ref_seq
        )

        output_file_path = Path(output_file)

        assert output_file_path.exists(), "The output file was not created"
        assert output_file_path.stat().st_size > 0, "The output file is empty"

        with output_file_path.open() as f:
            cleartext = f.read()

        if save_to_another_dir:
            another_dir = Path("tests/data/bam")
            another_dir.mkdir(exist_ok=True)
            another_output_file = another_dir / "REF_aln_trim_subsample_cleartext.fasta"
            another_output_file.write_text(cleartext)
            print(f"Output file also saved to: {another_output_file}")

        assert (
            cleartext == expected_cleartext
        ), "The converted cleartext data does not match the expected cleartext data"
