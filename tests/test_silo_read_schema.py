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

# Test data for other required components
VALID_ALIGNED_NUC_SEQ = {"main": "ACGTACGT"}

VALID_UNALIGNED_NUC_SEQ = {"main": "ACGTACGT"}

VALID_AA_SEQ = {"S": "KLMAST", "ORF1a": "VWXYZ"}

# Sample valid metadata dictionary from test_metadata_schema.py
VALID_METADATA = {
    "read_id": "read123",
    "sample_id": "A1_05_2024_10_08",
    "batch_id": "batch001",
    "sampling_date": "2024-10-08",
    "sequencing_date": "2024-10-24",
    "location_name": "Lugano (TI)",
    "read_length": "250",
    "primer_protocol": "v532",
    "location_code": "05",
    "flow_cell_serial_number": "2411515907",
    "sequencing_well_position": "A1",
    "primer_protocol_name": "SARS-CoV-2 ARTIC V5.3.2",
    "nextclade_reference": "sars-cov-2",
    "sr2silo_version": "v1.0.0",
}

# Test full AlignedReadSchema
VALID_ALIGNED_READ = {
    "metadata": VALID_METADATA,
    "nucleotideInsertions": VALID_NUC_INSERTIONS,
    "aminoAcidInsertions": VALID_AA_INSERTIONS,
    "alignedNucleotideSequences": VALID_ALIGNED_NUC_SEQ,
    "unalignedNucleotideSequences": VALID_UNALIGNED_NUC_SEQ,
    "alignedAminoAcidSequences": VALID_AA_SEQ,
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
    # Verify nucleotide insertions
    assert aligned_read.nucleotideInsertions.main == VALID_NUC_INSERTIONS["main"]
    # Verify amino acid insertions
    assert aligned_read.aminoAcidInsertions.root == VALID_AA_INSERTIONS
    # Verify metadata if present
    assert aligned_read.metadata is not None
    assert aligned_read.metadata.read_id == VALID_METADATA["read_id"]


def test_aligned_read_schema_without_metadata():
    """Test that an aligned read schema without metadata is accepted."""
    data = dict(VALID_ALIGNED_READ)
    data.pop("metadata")
    aligned_read = AlignedReadSchema(**data)
    assert aligned_read.metadata is None


def test_invalid_aligned_read_schema_missing_required():
    """Test that an aligned read schema missing required fields is rejected."""
    data = dict(VALID_ALIGNED_READ)
    data.pop("nucleotideInsertions")
    with pytest.raises(ValidationError):
        AlignedReadSchema(**data)


def test_invalid_aligned_read_schema_with_invalid_insertions():
    """Test that an aligned read schema with invalid insertions is rejected."""
    data = dict(VALID_ALIGNED_READ)
    data["nucleotideInsertions"] = INVALID_NUC_INSERTIONS_SEQUENCE
    with pytest.raises(ValidationError):
        AlignedReadSchema(**data)
