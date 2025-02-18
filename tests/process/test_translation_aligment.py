"""Tests for the translation module."""

from __future__ import annotations

import copy
import logging
from pathlib import Path
from typing import Dict

import pytest

from sr2silo.process.interface import Gene, GeneName, GeneSet, AAInsertionSet,AAInsertion, AASequenceSet, AlignedRead
import sr2silo.process.translate_align as translate_align
from sr2silo.process import convert

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def test_translate():
    """Test the translation function."""

    translate_align.translate_nextclade(
        [Path("tests/data/merged_expected.fasta")],
        Path("output/"),
        "nextstrain/sars-cov-2/XBB",
    )
    assert True


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


def test_parse_translate_align(aligned_reads):
    """Test the parse_translate_align function.

    Current dataset misses Amino Acid Insertions - i.e. not tested here.
    """

    # load the expected aligned reads
    expected_aligned_reads = []
    with open("tests/data/process/aligned_reads.ndjson") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            expected_aligned_reads.append(AlignedRead.from_str(line))

    # Convert actual_aligned_reads from dict to list for comparison.
    actual_aligned_reads = list(aligned_reads.values())

    # Compare each read attribute using their dictionary representations.
    for i, (actual, expected) in enumerate(
        zip(actual_aligned_reads, expected_aligned_reads)
    ):
        for attr in actual.to_dict().keys():
            act_attr = actual.to_dict()[attr]
            exp_attr = expected.to_dict()[attr]
            assert act_attr == exp_attr, (
                f"Mismatch at read {i} for attribute {attr}: "
                f"Expected {exp_attr} but got {act_attr}"
            )


@pytest.mark.skip(reason="Not implemented")
def test_read_in_AligendReads_nuc_seq():
    """Test the read_in_AlignedReads_nuc_seq function."""
    raise NotImplementedError


def test_read_in_AligendReads_nuc_ins(aligned_reads):
    """Test the read_in_AlignedReads_nuc_ins function."""

    fasta_nuc_insertions_file = Path("tests/data/process/nuc_insertions.fasta")

    expected_aligned_reads = translate_align.read_in_AlignedReads_nuc_ins(
        aligned_reads, fasta_nuc_insertions_file
    )

    # make list of aligned reads to compare
    expected_aligned_reads = list(expected_aligned_reads.values())

    actual_aligned_reads = copy.deepcopy(aligned_reads)
    # set the nucitides to empty lists for comparison
    for read in actual_aligned_reads.values():
        read.nucleotide_insertions = []

    # Convert actual_aligned_reads from dict to list for comparison.
    actual_aligned_reads = list(aligned_reads.values())

    # Compare each read's read_id and nucleotide insertions.
    for i, (actual, expected) in enumerate(
        zip(actual_aligned_reads, expected_aligned_reads)
    ):
        assert actual.read_id == expected.read_id, (
            f"Mismatch at read {i} for read_id: "
            f"Expected {expected.read_id} but got {actual.read_id}"
        )
        assert actual.nucleotide_insertions == expected.nucleotide_insertions, (
            f"Mismatch at read {i} for nucleotide_insertions: "
            f"Expected {expected.nucleotide_insertions} "
            f"but got {actual.nucleotide_insertions}"
        )


@pytest.mark.skip(reason="Not implemented")
def test_read_in_AligendReads_aa_ins():
    """Test the read_in_AlignedReads_aa_ins function."""
    raise NotImplementedError


def test_read_in_AlignedReads_aa_seq_and_ins(tmp_path, monkeypatch):
    """Test read_in_AlignedReads_aa_seq_and_ins including amino acid insertions."""
    # Create a dummy GeneSet with one gene "GeneX" of length 50.
    gene_name_str = "GeneX"
    gene_name = GeneName(gene_name_str)
    dummy_gene = Gene(gene_name, 50)
    dummy_gene_set = GeneSet([dummy_gene])

    # Prepare a dummy AlignedRead object; other attributes can be empty.
    dummy_aligned_read = AlignedRead(
        read_id="read1",
        unaligned_nucleotide_sequences="",
        aligned_nucleotide_sequences="",
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([gene_name_str]),
        aligned_amino_acid_sequences=AASequenceSet([gene_name_str])
    )
    aligned_reads = {"read1": dummy_aligned_read}

    # Create a temporary FASTA-like file with tab-separated fields.
    # Fields: [read_id, ignored, gene_name, pos, ignored, cigar, ignored, ignored, ignored, dummy_seq]
    # We use:
    # read_id: "read1"
    # gene_name: "GeneX"
    # pos: 10
    # cigar: "5M1I5M1D3M"
    # dummy_seq: "ABCDEFGHIJ" (ignored since we patch convert.sam_to_seq_and_indels)
    content = "read1\tX\tGeneX\t10\tY\t5M1I5M1D3M\tZ\tW\tV\tABCDEFGHIJ\n"
    dummy_file = tmp_path / "dummy_aa_alignment.txt"
    dummy_file.write_text(content)

    # Patch sam_to_seq_and_indels to return fixed values:
    # Let it return: cleartext "CLEAR", insertions [(15, "INS")], deletions [(30, 1)]
    def fake_sam_to_seq_and_indels(seq, cigar):
        return ("CLEAR", [(15, "INS")], [(30, 1)])
    monkeypatch.setattr("sr2silo.process.convert.sam_to_seq_and_indels", fake_sam_to_seq_and_indels)

    # Call the function under test.
    updated_reads = translate_align.read_in_AlignedReads_aa_seq_and_ins(aligned_reads, dummy_file, dummy_gene_set)

    # Expected padded sequence for AA alignment:
    # pad_alignment("CLEAR", pos=10, reference_length=50)
    expected_padded = convert.pad_alignment("CLEAR", 10, 50)
    # Expected AA insertion: from fake_sam_to_seq_and_indels insertion tuple (15, "INS")
    expected_insertion = AAInsertion(15, "INS")
    expected_ins_set = AAInsertionSet([gene_name])
    expected_ins_set.set_insertions_for_gene(gene_name_str, [expected_insertion])

    # Retrieve the AA insertions for gene "GeneX" from updated_reads["read1"]
    aa_insertion_set = updated_reads.get("read1").amino_acid_insertions
    assert str(gene_name_str) in aa_insertion_set.to_dict(), "GeneX not present in AA insertion set"


    # Check that insertions equal a list containing expected_insertion as string representation.
    assert str(expected_ins_set) == str(aa_insertion_set), f"Expected insertion {str(expected_ins_set)} not found; got {aa_insertion_set}"

    # Retrieve AA sequences and check they were set correctly.
    aa_seq_dict = updated_reads["read1"].aligned_amino_acid_sequences.to_dict()
    expected_padded_str = expected_padded
    assert aa_seq_dict.get(gene_name_str) == expected_padded_str, (
        f"Expected padded sequence {expected_padded_str}, got {aa_seq_dict.get(gene_name_str)}"
    )
