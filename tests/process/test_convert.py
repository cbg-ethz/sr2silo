"""
This module contains tests for the conversion functions in the sr2silo package.
"""

from __future__ import annotations

from typing import Dict

import pytest

from sr2silo.process import bam_to_sam
from sr2silo.process.convert import (
    bam_to_fastq_handle_indels,
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

    # assert the types of the outputs
    assert isinstance(cleartext, str), "cleartext is not a string"
    assert isinstance(insertions, list), "insertions is not a list"
    assert isinstance(deletions, list), "deletions is not a list"

    # assert the elements of insertiosn and deletions
    for ins in insertions:
        assert isinstance(
            ins, AAInsertion
        ), "insertions contains a non-AAInsertion object"
    for del_ in deletions:
        assert isinstance(del_, tuple), "deletions contains a non-tuple object"

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


def test_get_gene_set_from_ref():
    """Test the get_gene_set_from_ref function using the real AA reference file."""
    from pathlib import Path

    from sr2silo.process.convert import get_gene_set_from_ref

    aa_ref_fp = Path("resources/sars-cov-2/aa_reference_genomes.fasta")
    gene_set = get_gene_set_from_ref(aa_ref_fp)
    gene_names = gene_set.get_gene_name_list()

    # Assert that expected gene names are present and their lengths are > 0.
    for expected in ["E", "S", "ORF1a", "ORF7a"]:
        assert expected in gene_names, f"Expected gene {expected} to be in gene set"
        gene_length = gene_set.get_gene_length(expected)
        assert (
            gene_length > 0
        ), f"Expected gene length for {expected} to be > 0, got {gene_length}"


def test_bam_to_fastq_handle_indels(dummy_alignment, tmp_path):
    """
    Test bam_to_fastq_handle_indels using the dummy AlignmentFile from conf_test.py.
    The dummy read has:
      - query_sequence "ACTG"
      - query_qualities [30, 31, 32, 33]
      - cigartuples: [(0,2), (1,1), (0,1)]
      - reference_start: 100
    Expected FASTQ record:
      @read1
      ACG
      +
      ?@B
      alignment_position:100
    Expected insertion file:
      read1    102    T    A
    """
    # Create temporary files for FASTQ and insertions
    fastq_file = tmp_path / "output.fastq"
    insertions_file = tmp_path / "insertions.txt"

    # dummy_alignment is provided by the fixture
    bam_to_fastq_handle_indels(dummy_alignment, fastq_file, insertions_file)

    expected_fastq = (
        "@read1\n"
        "ACG\n"  # from match of 2 bases ("AC") and then match of 1 ("G")
        "+\n"
        "?@B\n"  # qualities: chr(30+33)="?" , chr(31+33)="@" , chr(33+33)="B"
        "alignment_position:100\n"
    )
    fastq_content = fastq_file.read_text()
    assert (
        fastq_content == expected_fastq
    ), f"FASTQ output mismatch:\nExpected:\n{expected_fastq}\nGot:\n{fastq_content}"

    expected_insertion = "read1\t102\tT\tA\n"
    insertion_content = insertions_file.read_text()
    assert insertion_content == expected_insertion, (
        f"Insertion output mismatch:\nExpected:\n{expected_insertion}\n"
        f"Got:\n{insertion_content}"
    )


def test_get_gene_set_from_ref_malformed_no_sequence(tmp_path):
    """Test get_gene_set_from_ref with a FASTA file that has header(s) but no
    sequence lines."""
    malformed = tmp_path / "malformed.fasta"
    # Write a header without a following sequence line.
    malformed.write_text(">GeneX\n")
    from sr2silo.process.convert import get_gene_set_from_ref

    gene_set = get_gene_set_from_ref(malformed)
    # Expect GeneSet to be empty since no sequence was provided.
    assert (
        gene_set.get_gene_name_list() == []
    ), "Expected empty gene set for header without sequence"

    def test_get_gene_set_from_ref_malformed_blank_lines(tmp_path):
        """Test get_gene_set_from_ref with a FASTA file that has multiple headers
        and blank sequence lines."""

    # Create file with headers with blank sequence lines
    content = """>GeneA
    >GeneB
    AGCTAGCT
    >GeneC

    """
    malformed.write_text(content)
    from sr2silo.process.convert import get_gene_set_from_ref

    gene_set = get_gene_set_from_ref(malformed)
    # Expect only GeneB to be added (GeneA and GeneC have blank sequences)
    gene_names = gene_set.get_gene_name_list()
    assert "GeneB" in gene_names, "Expected GeneB to be parsed"
    assert (
        "GeneA" not in gene_names
    ), "Expected GeneA to be skipped due to missing sequence"
    assert (
        "GeneC" not in gene_names
    ), "Expected GeneC to be skipped due to missing sequence"
