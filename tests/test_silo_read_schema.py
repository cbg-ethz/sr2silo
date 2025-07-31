"""Tests for SILO read schema validation, focusing on insertion schemas."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from sr2silo.silo_read_schema import (
    AlignedReadSchema,
    AminoAcidInsertions,
    NucleotideInsertions,
)

# Test data for NucleotideInsertions
VALID_NUC_INSERTIONS = {"main": ["10:ACGT", "123:A", "9999:TACG"]}

INVALID_NUC_INSERTIONS_FORMAT = {
    "main": ["10:ACGT", "abc:ACGT"]  # non-numeric position
}

INVALID_NUC_INSERTIONS_SEQUENCE = {"main": ["10:ACGTZ"]}  # invalid nucleotide 'Z'

# Test data for AminoAcidInsertions
VALID_AA_INSERTIONS = {
    "S": ["10:AST", "123:G"],
    "ORF1a": ["45:KLM", "67:Y"],
}

INVALID_AA_INSERTIONS_FORMAT = {"S": ["10:AST", "abc:KLM"]}  # non-numeric position

INVALID_AA_INSERTIONS_SEQUENCE = {"S": ["10:ASt"]}  # lowercase 't' is invalid

# Test data for new schema format
VALID_MAIN_SEGMENT = {
    "sequence": "ACGTACGTACGTACGT",
    "insertions": ["10:ACGT", "123:A"],
    "offset": 0,
}

VALID_GENE_SEGMENTS = {
    "S": {"sequence": "KLMAST", "insertions": ["10:AST", "123:G"], "offset": 0},
    "ORF1a": {"sequence": "VWXYZ", "insertions": ["45:KLM", "67:Y"], "offset": 0},
}

# Test full AlignedReadSchema with new format
VALID_ALIGNED_READ = {
    "main": VALID_MAIN_SEGMENT,
    "read_id": "read123",
    "sample_id": "A1_05_2024_10_08",
    "batch_id": "batch001",
    "sampling_date": "2024-10-08",
    "location_name": "Lugano (TI)",
    "read_length": "250",
    "unaligned_main": "ACGTACGTACGTACGT",
    **VALID_GENE_SEGMENTS,
}


def test_valid_nucleotide_insertions():
    """Test that valid nucleotide insertions are accepted."""
    insertions = NucleotideInsertions(**VALID_NUC_INSERTIONS)
    assert insertions.main == VALID_NUC_INSERTIONS["main"]


def test_invalid_nucleotide_insertions_format():
    """Test that nucleotide insertions with invalid position format are rejected."""
    with pytest.raises(ValidationError) as exc_info:
        NucleotideInsertions(**INVALID_NUC_INSERTIONS_FORMAT)
    assert "not in the expected format" in str(exc_info.value)


def test_invalid_nucleotide_insertions_sequence():
    """Test that nucleotide insertions with invalid sequence are rejected."""
    with pytest.raises(ValidationError) as exc_info:
        NucleotideInsertions(**INVALID_NUC_INSERTIONS_SEQUENCE)
    assert "not in the expected format" in str(exc_info.value)


def test_valid_nucleotide_insertions_with_ns():
    """Test that nucleotide insertions with 'N' in the sequence are accepted."""
    valid_insertions = {"main": ["10:NNNN"]}
    insertions = NucleotideInsertions(**valid_insertions)
    assert insertions.main == valid_insertions["main"]


def test_valid_amino_acid_insertions():
    """Test that valid amino acid insertions are accepted."""
    insertions = AminoAcidInsertions(root=VALID_AA_INSERTIONS)
    assert insertions.root == VALID_AA_INSERTIONS


def test_invalid_amino_acid_insertions_format():
    """Test that amino acid insertions with invalid position format are rejected."""
    with pytest.raises(ValidationError) as exc_info:
        AminoAcidInsertions(root=INVALID_AA_INSERTIONS_FORMAT)
    assert "not in the expected format" in str(exc_info.value)


def test_invalid_amino_acid_insertions_sequence():
    """Test that amino acid insertions with invalid sequence are rejected."""
    with pytest.raises(ValidationError) as exc_info:
        AminoAcidInsertions(root=INVALID_AA_INSERTIONS_SEQUENCE)
    assert "not in the expected format" in str(exc_info.value)


def test_valid_amino_acid_insertions_with_ns():
    """Test that amino acid insertions with 'N' in the sequence are accepted."""
    valid_insertions = {"S": ["10:NNN"]}
    insertions = AminoAcidInsertions(root=valid_insertions)
    assert insertions.root == valid_insertions


def test_valid_aligned_read_schema():
    """Test that a valid complete aligned read schema is accepted."""
    aligned_read = AlignedReadSchema(**VALID_ALIGNED_READ)
    # Verify main nucleotide segment
    assert aligned_read.main.sequence == VALID_MAIN_SEGMENT["sequence"]
    assert aligned_read.main.insertions == VALID_MAIN_SEGMENT["insertions"]
    assert aligned_read.main.offset == VALID_MAIN_SEGMENT["offset"]
    # Verify gene segments are accessible via getattr or dict access
    assert hasattr(aligned_read, "S")
    assert hasattr(aligned_read, "ORF1a")
    # Verify metadata fields are accessible
    assert getattr(aligned_read, "read_id") == "read123"
    assert getattr(aligned_read, "sample_id") == "A1_05_2024_10_08"


def test_aligned_read_schema_without_metadata():
    """Test that an aligned read schema with minimal required fields is accepted."""
    from sr2silo.silo_read_schema import NucleotideSegment

    minimal_data = {"main": NucleotideSegment(**VALID_MAIN_SEGMENT)}
    aligned_read = AlignedReadSchema(**minimal_data)
    assert aligned_read.main.sequence == VALID_MAIN_SEGMENT["sequence"]


def test_invalid_aligned_read_schema_missing_required():
    """Test that an aligned read schema missing required fields is rejected."""
    data = dict(VALID_ALIGNED_READ)
    data.pop("main")  # Remove the required main field
    with pytest.raises(ValidationError):
        AlignedReadSchema(**data)


def test_invalid_aligned_read_schema_with_invalid_insertions():
    """Test that an aligned read schema with invalid insertions is rejected."""
    data = dict(VALID_ALIGNED_READ)
    # Create invalid main segment with bad insertions
    invalid_main = dict(VALID_MAIN_SEGMENT)
    invalid_main["insertions"] = ["abc:ACGT"]  # non-numeric position
    data["main"] = invalid_main
    with pytest.raises(ValidationError):
        AlignedReadSchema(**data)
