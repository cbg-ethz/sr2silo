"""Tests for the translation module."""

from __future__ import annotations

import copy
import logging
from pathlib import Path

import pytest

import sr2silo.process.translate_align as translate_align
from sr2silo.process.interface import (
    AAInsertion,
    AAInsertionSet,
    AASequenceSet,
    AlignedRead,
    Gene,
    GeneName,
    GeneSet,
)

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def test_read_in_AligendReads_nuc_ins(aligned_reads):
    """Test the enrich_read_with_nuc_ins function."""
    fasta_nuc_insertions_file = Path("tests/data/process/nuc_insertions.fasta")

    expected_aligned_reads = translate_align.enrich_read_with_nuc_ins(
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


def test_enrich_read_with_aa_seq_basic(tmp_path, monkeypatch):
    """Test that enrich_read_with_aa_seq correctly processes AA alignment data."""
    # Create a simple AlignedRead
    read_id = "test_read"
    gene_name = "TestGene"

    aligned_read = AlignedRead(
        read_id=read_id,
        unaligned_nucleotide_sequence="ATCG",
        aligned_nucleotide_sequence="ATCG",
        aligned_nucleotide_sequence_offset=0,
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([GeneName(gene_name)]),
        aligned_amino_acid_sequences=AASequenceSet([GeneName(gene_name)]),
    )
    aligned_reads = {read_id: aligned_read}

    # Create simple gene set
    gene_set = GeneSet([Gene(GeneName(gene_name), 100)])

    # Create mock SAM-like alignment file
    sam_line = f"{read_id}\t0\t{gene_name}\t5\t60\t3M\t*\t0\t0\tMKT\t!!!\n"
    alignment_file = tmp_path / "test_alignment.sam"
    alignment_file.write_text(sam_line)

    # Mock the convert function to return predictable results
    def mock_sam_to_seq_and_indels(seq, cigar):
        return ("MKT", [AAInsertion(2, "X")])

    monkeypatch.setattr(
        "sr2silo.process.convert.sam_to_seq_and_indels", mock_sam_to_seq_and_indels
    )

    # Test the function
    result = translate_align.enrich_read_with_aa_seq(
        aligned_reads, alignment_file, gene_set
    )

    # Verify results
    updated_read = result[read_id]

    # Check sequence was set
    sequences = updated_read.aligned_amino_acid_sequences.to_dict()
    assert sequences[gene_name] == "MKT"

    # Check insertions were set
    insertions = updated_read.get_amino_acid_insertions().to_dict()
    assert insertions[gene_name] == ["2:X"]


def test_curry_read_with_metadata(tmp_path):
    """Test that curry_read_with_metadata correctly enriches a read with metadata."""
    # Create minimal test metadata
    metadata = {"sample_id": "test123", "location": "TestLab"}

    metadata_file = tmp_path / "metadata.json"
    with open(metadata_file, "w") as f:
        import json

        json.dump(metadata, f)

    # Create minimal AlignedRead
    read = AlignedRead(
        read_id="test_read",
        unaligned_nucleotide_sequence="ACGT",
        aligned_nucleotide_sequence="ACGT",
        aligned_nucleotide_sequence_offset=0,
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([]),
        aligned_amino_acid_sequences=AASequenceSet([]),
    )

    # Test the main functionality
    enrich_func = translate_align.curry_read_with_metadata(metadata_file)
    enriched_read = enrich_func(read)

    assert enriched_read.metadata == metadata


def test_curry_read_with_metadata_file_not_found(tmp_path):
    """Test that curry_read_with_metadata raises FileNotFoundError for missing files."""
    with pytest.raises(FileNotFoundError):
        translate_align.curry_read_with_metadata(tmp_path / "nonexistent.json")


def test_curry_read_with_metadata_invalid_json(tmp_path):
    """Test that curry_read_with_metadata raises JSONDecodeError for invalid JSON."""
    import json

    invalid_file = tmp_path / "invalid.json"
    with open(invalid_file, "w") as f:
        f.write("not valid json")

    with pytest.raises(json.JSONDecodeError):
        translate_align.curry_read_with_metadata(invalid_file)


def test_curry_read_with_metadata_empty_metadata(tmp_path):
    """Test that curry_read_with_metadata raises ValueError for empty metadata."""
    empty_file = tmp_path / "empty.json"
    with open(empty_file, "w") as f:
        import json

        json.dump({}, f)

    with pytest.raises(ValueError, match="No metadata found in the file"):
        translate_align.curry_read_with_metadata(empty_file)
