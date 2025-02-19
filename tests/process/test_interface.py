"""A placeholder test file for the interface module."""

from __future__ import annotations

import pytest


# Placeholder tests for NucInsertion
@pytest.mark.skip(reason="Not yet implemented")
def test_nuc_insertion():
    """Test NucInsertion functionality."""
    raise NotImplementedError


# Placeholder tests for AAInsertion
@pytest.mark.skip(reason="Not yet implemented")
def test_aa_insertion():
    """Test AAInsertion functionality."""
    raise NotImplementedError


# Placeholder tests for AlignedRead
@pytest.mark.skip(reason="Not yet implemented")
def test_aligned_read():
    """Test AlignedRead functionality."""
    raise NotImplementedError


# Placeholder tests for GeneName
@pytest.mark.skip(reason="Not yet implemented")
def test_gene_name():
    """Test GeneName functionality."""
    raise NotImplementedError


# Placeholder tests for Gene
@pytest.mark.skip(reason="Not yet implemented")
def test_gene():
    """Test Gene functionality."""
    raise NotImplementedError


# Placeholder tests for GeneSet
@pytest.mark.skip(reason="Not yet implemented")
def test_gene_set():
    """Test GeneSet functionality."""
    raise NotImplementedError


# Placeholder tests for AAInsertionSet
@pytest.mark.skip(reason="Not yet implemented")
def test_aa_insertion_set():
    """Test AAInsertionSet functionality."""
    raise NotImplementedError


# Placeholder tests for AASequenceSet
@pytest.mark.skip(reason="Not yet implemented")
def test_aa_sequence_set():
    """Test AASequenceSet functionality."""
    raise NotImplementedError


# test the to_silo_json method
def test_to_silo_json(aligned_reads):
    """Test to_silo_json functionality."""

    # for all reads in aligned_reads, test the to_silo_json method
    for read in aligned_reads.values():
        print(read.to_silo_json())

    raise NotImplementedError
