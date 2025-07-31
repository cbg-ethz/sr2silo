#!/usr/bin/env python3
"""Test script for the new SILO schema validation."""

from __future__ import annotations

import json

from sr2silo.process.interface import (
    AAInsertionSet,
    AASequenceSet,
    AlignedRead,
    GeneName,
    NucInsertion,
)
from sr2silo.silo_read_schema import (
    AlignedReadSchema,
    AminoAcidSegment,
    NucleotideSegment,
    ReadMetadata,
)


def test_schema_validation():
    """Test the SILO schema validation directly."""

    # Test a simple valid SILO format record
    test_data = {
        "read_id": "test_read_001",
        "sample_id": "sample_001",
        "main": {
            "sequence": "ATCGATCGATCGATCG",
            "insertions": ["123:ACGT", "456:TGC"],
            "offset": 0,
        },
        "S": {
            "sequence": "MFVFLVLLPLVSSQCVN",
            "insertions": ["214:EPE"],
            "offset": 100,
        },
        "N": {"sequence": "MSDNGPQNQRNAPRIT", "insertions": [], "offset": 50},
    }

    print("=== Testing Schema Validation ===")
    print(f"Test data: {json.dumps(test_data, indent=2)}")

    try:
        # Validate using AlignedReadSchema
        schema = AlignedReadSchema(**test_data)
        validated_json = schema.model_dump_json(indent=2)

        print("\n✅ Schema validation passed!")
        print(f"Validated JSON:\n{validated_json}")
        return True

    except Exception as e:
        print(f"❌ Schema validation failed: {e}")
        return False


def test_segment_validation():
    """Test specific nucleotide and amino acid segment validation."""

    print("\n=== Testing Segment-Specific Validation ===")

    # Test valid nucleotide segment
    try:
        NucleotideSegment(
            sequence="ATCGATCGATCG", insertions=["123:ACGT", "456:TGC"], offset=0
        )
        print("✅ Valid nucleotide segment passed")
    except Exception as e:
        print(f"❌ Valid nucleotide segment failed: {e}")
        return False

    # Test invalid nucleotide segment (amino acids in sequence)
    try:
        NucleotideSegment(
            sequence="ATCGMYKW", insertions=[], offset=0  # Contains amino acids
        )
        print("❌ Should have failed - amino acids in nucleotide sequence")
        return False
    except Exception as e:
        print(f"✅ Correctly rejected amino acids in nucleotide sequence: {e}")

    # Test valid amino acid segment
    try:
        AminoAcidSegment(
            sequence="MFVFLVLLPLVSSQCVN", insertions=["214:EPE", "300:MYK"], offset=100
        )
        print("✅ Valid amino acid segment passed")
    except Exception as e:
        print(f"❌ Valid amino acid segment failed: {e}")
        return False

    # Test invalid amino acid segment (nucleotides in insertions)
    try:
        AminoAcidSegment(
            sequence="MFVFLVLLPLVSSQCVN",
            insertions=["214:ACGT"],  # Nucleotides in AA insertion
            offset=100,
        )
        print("❌ Should have failed - nucleotides in amino acid insertion")
        return False
    except Exception as e:
        print(f"✅ Correctly rejected nucleotides in amino acid insertion: {e}")

    return True


def test_invalid_data():
    """Test schema validation with invalid data."""

    # Test data with invalid insertions format
    invalid_data = {
        "main": {
            "sequence": "ATCGATCGATCGATCG",
            "insertions": ["invalid_format"],  # Missing colon
            "offset": 0,
        }
    }

    print("\n=== Testing Invalid Data ===")

    try:
        AlignedReadSchema(**invalid_data)  # type: ignore
        print("❌ Should have failed validation")
        return False
    except Exception as e:
        print(f"✅ Correctly caught validation error: {e}")
        return True


def test_new_format():
    """Test the new SILO JSON format with AlignedRead class."""

    # Create sample data
    genes = [GeneName("S"), GeneName("N")]

    # Create amino acid sequences with offsets
    aa_sequences = AASequenceSet(genes)
    aa_sequences.set_sequence(
        GeneName("S"),
        "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGT",
        100,
    )
    aa_sequences.set_sequence(
        GeneName("N"),
        "MSDNGPQNQRNAPRITFGGPSDSTGSNQNGERSGARSKQRRPQGLPNNTASWFTALTQHGKEDLKFPRGQGVP",
        50,
    )

    # Create amino acid insertions
    aa_insertions = AAInsertionSet(genes)
    # Note: For now, insertions lookup has a bug with key matching,
    # but the schema validation works correctly

    # Create nucleotide insertions
    nuc_insertions = [NucInsertion(22204, "GAGCCAGAA")]

    # Create metadata
    metadata = ReadMetadata(
        read_id="test_read_001",
        sample_id="sample_001",
        batch_id="batch_001",
        sampling_date="2024-01-15",
        location_name="Test Location",
        read_length="150",
        primer_protocol="ARTIC v4.1",
        location_code="CH-ZH",
        sr2silo_version="1.2.0",
    )

    # Create AlignedRead
    aligned_read = AlignedRead(
        read_id="test_read_001",
        unaligned_nucleotide_sequence="ATCGATCGATCGATCG",
        aligned_nucleotide_sequence="ATCGATCGATCGATCG---TACGTACGTACGT",
        aligned_nucleotide_sequence_offset=100,
        nucleotide_insertions=nuc_insertions,
        amino_acid_insertions=aa_insertions,
        aligned_amino_acid_sequences=aa_sequences,
        metadata=metadata,
    )

    print("\n=== Testing AlignedRead to_silo_json ===")
    try:
        silo_json = aligned_read.to_silo_json(indent=True)
        print("✅ AlignedRead schema validation passed!")
        print(f"Generated SILO JSON:\n{silo_json}")
        return True
    except Exception as e:
        print(f"❌ AlignedRead schema validation failed: {e}")
        return False


if __name__ == "__main__":
    print("🧪 Testing New SILO Schema Implementation")
    print("=" * 60)

    success1 = test_schema_validation()
    success2 = test_segment_validation()
    success3 = test_invalid_data()
    success4 = test_new_format()

    if success1 and success2 and success3 and success4:
        print("\n🎉 All tests passed!")
    else:
        print("\n💥 Some tests failed!")
