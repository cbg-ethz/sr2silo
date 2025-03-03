"""A tests for the interface module."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from sr2silo.process.interface import (
    AAInsertion,
    AAInsertionSet,
    AASequenceSet,
    AlignedRead,
    Gene,
    GeneName,
    GeneSet,
    NucInsertion,
)
from sr2silo.silo_read_schema import ReadMetadata


def test_nuc_insertion():
    """Test NucInsertion functionality."""
    insertion = NucInsertion(10, "ACTG")
    assert insertion.position == 10
    assert insertion.sequence == "ACTG"
    assert str(insertion) == "10 : ACTG"


def test_aa_insertion():
    """Test AAInsertion functionality."""
    insertion = AAInsertion(5, "MKT")
    assert insertion.position == 5
    assert insertion.sequence == "MKT"
    assert str(insertion) == "5 : MKT"


def test_aligned_read():
    """Test AlignedRead functionality."""
    read = AlignedRead(
        read_id="read1",
        unaligned_nucleotide_sequences="ACTG",
        aligned_nucleotide_sequences="ACTG",
        nucleotide_insertions=[NucInsertion(10, "ACTG")],
        amino_acid_insertions=AAInsertionSet([GeneName("gene1")]),
        aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
    )
    assert read.read_id == "read1"
    assert read.unaligned_nucleotide_sequences == "ACTG"
    assert read.aligned_nucleotide_sequences == "ACTG"
    assert len(read.nucleotide_insertions) == 1
    assert isinstance(read.nucleotide_insertions[0], NucInsertion)
    assert isinstance(read.amino_acid_insertions, AAInsertionSet)
    assert isinstance(read.aligned_amino_acid_sequences, AASequenceSet)


def test_gene_name():
    """Test GeneName functionality."""
    gene_name = GeneName("gene1")
    assert gene_name.name == "gene1"
    assert str(gene_name) == "gene1"


def test_gene():
    """Test Gene functionality."""
    gene_name = GeneName("gene1")
    gene = Gene(gene_name, 1000)
    assert gene.name == gene_name
    assert gene.gene_length == 1000
    assert gene.to_dict() == {"gene_name": gene_name, "gene_length": 1000}


def test_gene_set():
    """Test GeneSet functionality."""
    gene1 = Gene(GeneName("gene1"), 1000)
    gene2 = Gene(GeneName("gene2"), 2000)
    gene_set = GeneSet([gene1, gene2])
    assert gene_set.get_gene(GeneName("gene1")) == gene1
    assert gene_set.get_gene_length(GeneName("gene2")) == 2000
    assert gene_set.to_dict() == {
        "gene1": {"gene_name": "gene1", "gene_length": 1000},
        "gene2": {"gene_name": "gene2", "gene_length": 2000},
    }


def test_aa_insertion_set():
    """Test AAInsertionSet functionality."""
    gene_name = GeneName("gene1")
    aa_insertion_set = AAInsertionSet([gene_name])
    aa_insertion_set.set_insertions_for_gene(gene_name, [AAInsertion(5, "MKT")])
    assert aa_insertion_set.to_dict() == {"gene1": ["5 : MKT"]}


def test_aa_sequence_set():
    """Test AASequenceSet functionality."""
    gene_name = GeneName("gene1")
    aa_sequence_set = AASequenceSet([gene_name])
    aa_sequence_set.set_sequence(gene_name, "MKT")
    assert aa_sequence_set.to_dict() == {"gene1": "MKT"}


def test_to_silo_json():
    """Test to_silo_json functionality."""
    aligned_reads = {
        "read1": AlignedRead(
            read_id="read1",
            unaligned_nucleotide_sequences="ACTG",
            aligned_nucleotide_sequences="ACTG",
            nucleotide_insertions=[NucInsertion(10, "ACTG")],
            amino_acid_insertions=AAInsertionSet([GeneName("gene1")]),
            aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
        )
    }
    # add mock metadata to aligned_reads
    metadata = {
        "sequencing_date": "2024-10-18",
        "location_name": "Zürich (ZH)",
        "batch_id": "20241018_AAG55WNM5",
        "read_length": "250",
        "primer_protocol": "v532",
        "location_code": "10",
        "flow_cell_serial_number": "AAG55WNM5",
        "nextclade_reference": "sars-cov-2",
        "sequencing_well_position": "A1",
        "sample_id": "A1_10_2024_09_30",
        "sampling_date": "2024-09-30",
        "primer_protocol_name": "SARS-CoV-2 ARTIC V5.3.2",
    }

    for read_id, read in aligned_reads.items():
        read_metadata = metadata.copy()
        read_metadata["read_id"] = read_id
        validated_metadata = ReadMetadata(**read_metadata)
        read.set_metadata(validated_metadata)

    # for all reads in aligned_reads, test the to_silo_json method
    for read in aligned_reads.values():
        try:
            read.to_silo_json()
        except ValidationError as e:
            pytest.fail(f"Validation error: {e}")


def test_set_nuc_insertion():
    """Test setting additional nucleotide insertion in AlignedRead."""
    from sr2silo.process.interface import (
        AAInsertionSet,
        AASequenceSet,
        AlignedRead,
        GeneName,
        NucInsertion,
    )

    read = AlignedRead(
        read_id="insertion_test",
        unaligned_nucleotide_sequences="AAAA",
        aligned_nucleotide_sequences="AAAA",
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([GeneName("gene1")]),
        aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
    )
    initial_len = len(read.nucleotide_insertions)
    new_insertion = NucInsertion(20, "TTTT")
    read.set_nuc_insertion(new_insertion)
    assert len(read.nucleotide_insertions) == initial_len + 1
    assert read.nucleotide_insertions[-1].position == 20
    assert read.nucleotide_insertions[-1].sequence == "TTTT"


def test_get_amino_acid_insertions():
    """Test get_amino_acid_insertions method."""
    from sr2silo.process.interface import (
        AAInsertion,
        AAInsertionSet,
        AASequenceSet,
        AlignedRead,
        GeneName,
    )

    aa_ins_set = AAInsertionSet([GeneName("gene1")])
    aa_ins_set.set_insertions_for_gene(GeneName("gene1"), [AAInsertion(7, "ABC")])
    read = AlignedRead(
        read_id="aa_ins_test",
        unaligned_nucleotide_sequences="CCCC",
        aligned_nucleotide_sequences="CCCC",
        nucleotide_insertions=[],
        amino_acid_insertions=aa_ins_set,
        aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
    )
    retrieved = read.get_amino_acid_insertions().to_dict()
    assert "gene1" in retrieved
    assert retrieved["gene1"] == ["7 : ABC"]


def test_get_metadata_without_setting():
    """Test get_metadata when metadata is not set."""
    from sr2silo.process.interface import (
        AAInsertionSet,
        AASequenceSet,
        AlignedRead,
        GeneName,
    )

    read = AlignedRead(
        read_id="meta_test",
        unaligned_nucleotide_sequences="GGGG",
        aligned_nucleotide_sequences="GGGG",
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([GeneName("gene1")]),
        aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
    )
    metadata = read.get_metadata()
    assert metadata is None


def test_gene_set_get_gene_name_list():
    """Test GeneSet.get_gene_name_list returns gene names as strings."""
    from sr2silo.process.interface import Gene, GeneName, GeneSet

    gene1 = Gene(GeneName("geneA"), 500)
    gene2 = Gene(GeneName("geneB"), 1500)
    gene_set = GeneSet([gene1, gene2])
    names = gene_set.get_gene_name_list()
    # The names are stored as strings
    assert "geneA" in names
    assert "geneB" in names


def test_str_methods():
    """Test __str__ methods for various classes."""

    # Test NucInsertion.__str__
    nuc = NucInsertion(30, "GGCC")
    assert str(nuc) == "30 : GGCC"

    # Test AAInsertion.__str__
    aa = AAInsertion(10, "MLK")
    assert str(aa) == "10 : MLK"

    # Test AAInsertionSet.__str__ equals its to_dict string
    aa_ins_set = AAInsertionSet([GeneName("gene1")])
    aa_ins_set.set_insertions_for_gene(GeneName("gene1"), [aa])
    assert str(aa_ins_set) == str(aa_ins_set.to_dict())

    # Test AASequenceSet.__str__ equals its to_dict string
    aa_seq_set = AASequenceSet([GeneName("gene1")])
    aa_seq_set.set_sequence(GeneName("gene1"), "MLKMLK")
    assert str(aa_seq_set) == str(aa_seq_set.to_dict())

    # Test Gene.__str__ via its to_dict conversion in GeneSet
    gene = Gene(GeneName("geneX"), 750)
    expected = "{gene_name: geneX, gene_length: 750}"
    assert str(gene) == expected
