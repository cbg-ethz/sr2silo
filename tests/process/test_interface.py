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
    assert str(insertion) == "10:ACTG"


def test_aa_insertion():
    """Test AAInsertion functionality."""
    insertion = AAInsertion(5, "MKT")
    assert insertion.position == 5
    assert insertion.sequence == "MKT"
    assert str(insertion) == "5:MKT"


def test_aligned_read():
    """Test AlignedRead functionality."""
    read = AlignedRead(
        read_id="read1",
        unaligned_nucleotide_sequence="ACTG",
        aligned_nucleotide_sequence="ACTG",
        aligned_nucleotide_sequence_offset=0,
        nucleotide_insertions=[NucInsertion(10, "ACTG")],
        amino_acid_insertions=AAInsertionSet([GeneName("gene1")]),
        aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
    )
    assert read.read_id == "read1"
    assert read.unaligned_nucleotide_sequence == "ACTG"
    assert read.aligned_nucleotide_sequence == "ACTG"
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
    assert aa_insertion_set.to_dict() == {"gene1": ["5:MKT"]}


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
            unaligned_nucleotide_sequence="ACTG",
            aligned_nucleotide_sequence="ACTG",
            aligned_nucleotide_sequence_offset=0,
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
        "sr2silo_version": "v1.0.0",
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
        unaligned_nucleotide_sequence="AAAA",
        aligned_nucleotide_sequence="AAAA",
        aligned_nucleotide_sequence_offset=0,
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
        unaligned_nucleotide_sequence="CCCC",
        aligned_nucleotide_sequence="CCCC",
        aligned_nucleotide_sequence_offset=0,
        nucleotide_insertions=[],
        amino_acid_insertions=aa_ins_set,
        aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
    )
    retrieved = read.get_amino_acid_insertions().to_dict()
    assert "gene1" in retrieved
    assert retrieved["gene1"] == ["7:ABC"]


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
        unaligned_nucleotide_sequence="GGGG",
        aligned_nucleotide_sequence="GGGG",
        aligned_nucleotide_sequence_offset=0,
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
    assert str(nuc) == "30:GGCC"

    # Test AAInsertion.__str__
    aa = AAInsertion(10, "MLK")
    assert str(aa) == "10:MLK"

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


def test_insertion_equality():
    """Test the __eq__ method for Insertion class."""
    # Test base Insertion class equality
    ins1 = NucInsertion(10, "ACTG")
    ins2 = NucInsertion(10, "ACTG")
    ins3 = NucInsertion(10, "ACTA")
    ins4 = NucInsertion(11, "ACTG")

    # Same position and sequence should be equal
    assert ins1 == ins2

    # Different sequence should not be equal
    assert ins1 != ins3

    # Different position should not be equal
    assert ins1 != ins4

    # Comparing with a non-Insertion object should return False
    assert ins1 != "not_an_insertion"

    # Test with AAInsertion
    aa_ins1 = AAInsertion(5, "MET")
    aa_ins2 = AAInsertion(5, "MET")
    aa_ins3 = AAInsertion(5, "ALA")

    assert aa_ins1 == aa_ins2
    assert aa_ins1 != aa_ins3

    # Test cross-subclass comparison (should be False even with same position/sequence)
    mixed_ins = AAInsertion(10, "ACTG")
    assert ins1 != mixed_ins


def test_aa_insertion_set_equality():
    """Test the __eq__ method for AAInsertionSet class."""
    # Create two identical sets
    gene_name1 = GeneName("gene1")
    gene_name2 = GeneName("gene2")

    aa_set1 = AAInsertionSet([gene_name1, gene_name2])
    aa_set1.set_insertions_for_gene(
        gene_name1, [AAInsertion(5, "MKT"), AAInsertion(10, "LVD")]
    )
    aa_set1.set_insertions_for_gene(gene_name2, [AAInsertion(15, "RST")])

    aa_set2 = AAInsertionSet([gene_name1, gene_name2])
    aa_set2.set_insertions_for_gene(
        gene_name1, [AAInsertion(5, "MKT"), AAInsertion(10, "LVD")]
    )
    aa_set2.set_insertions_for_gene(gene_name2, [AAInsertion(15, "RST")])

    # Sets with identical content should be equal
    assert aa_set1 == aa_set2

    # Sets with different insertion order should still be equal
    aa_set3 = AAInsertionSet([gene_name1, gene_name2])
    aa_set3.set_insertions_for_gene(
        gene_name1, [AAInsertion(10, "LVD"), AAInsertion(5, "MKT")]
    )
    aa_set3.set_insertions_for_gene(gene_name2, [AAInsertion(15, "RST")])

    assert aa_set1 == aa_set3

    # Sets with different content should not be equal
    aa_set4 = AAInsertionSet([gene_name1, gene_name2])
    aa_set4.set_insertions_for_gene(
        gene_name1, [AAInsertion(5, "MKT"), AAInsertion(10, "LVX")]
    )
    aa_set4.set_insertions_for_gene(gene_name2, [AAInsertion(15, "RST")])

    assert aa_set1 != aa_set4

    # Sets with different keys should not be equal
    aa_set5 = AAInsertionSet([gene_name1])
    aa_set5.set_insertions_for_gene(
        gene_name1, [AAInsertion(5, "MKT"), AAInsertion(10, "LVD")]
    )

    assert aa_set1 != aa_set5

    # Comparing with a non-AAInsertionSet object should return False
    assert aa_set1 != "not_an_aa_insertion_set"


def test_read_id_in_json_output():
    """Test that read_id appears as the first field in JSON output."""
    from sr2silo.process.interface import (
        AAInsertionSet,
        AASequenceSet,
        AlignedRead,
        GeneName,
    )

    read_id = "test_read_123"
    read = AlignedRead(
        read_id=read_id,
        unaligned_nucleotide_sequence="ACTG",
        aligned_nucleotide_sequence="ACTG",
        aligned_nucleotide_sequence_offset=0,
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([GeneName("gene1")]),
        aligned_amino_acid_sequences=AASequenceSet([GeneName("gene1")]),
    )

    # Test to_dict includes read_id
    result_dict = read.to_dict()

    # Check that read_id is present
    assert "read_id" in result_dict, "read_id missing from to_dict() output"
    assert result_dict["read_id"] == read_id, (
        f"Expected {read_id}, got {result_dict['read_id']}"
    )

    # Check that read_id is the first key
    first_key = list(result_dict.keys())[0]
    assert first_key == "read_id", (
        f"read_id should be first key, but {first_key} was first"
    )

    # Test that the SILO JSON validation works with read_id
    try:
        json_output = read.to_silo_json()
        assert json_output is not None

        # Parse the JSON to verify structure
        import json

        parsed = json.loads(json_output)
        assert "read_id" in parsed
        assert parsed["read_id"] == read_id
    except Exception as e:
        pytest.fail(f"SILO JSON validation failed with read_id: {e}")


def test_empty_gene_representation():
    """Test that genes with no sequence are represented as null."""
    from sr2silo.process.interface import (
        AAInsertionSet,
        AASequenceSet,
        AlignedRead,
        GeneName,
    )

    # Create a read with genes - one empty, one with content
    gene1 = GeneName("gene1")
    gene2 = GeneName("gene2")

    # Set up amino acid sequences - gene1 has sequence, gene2 is empty
    aa_seq_set = AASequenceSet([gene1, gene2])
    aa_seq_set.set_sequence(gene1, "MKTSF", 0)
    aa_seq_set.set_sequence(gene2, "", 0)  # Empty sequence

    # Set up amino acid insertions - note that insertions without sequence
    # context are non-sense
    aa_ins_set = AAInsertionSet([gene1, gene2])
    aa_ins_set.set_insertions_for_gene(gene1, [AAInsertion(5, "VLK")])
    aa_ins_set.set_insertions_for_gene(
        gene2, [AAInsertion(10, "ABC")]
    )  # Nonsensical without sequence

    read = AlignedRead(
        read_id="empty_gene_test",
        unaligned_nucleotide_sequence="ATGAAGACCTCGTTC",
        aligned_nucleotide_sequence="ATGAAGACCTCGTTC",
        aligned_nucleotide_sequence_offset=0,
        nucleotide_insertions=[],
        amino_acid_insertions=aa_ins_set,
        aligned_amino_acid_sequences=aa_seq_set,
    )

    # Get the dictionary representation
    result_dict = read.to_dict()

    # Check that gene1 is represented normally
    assert "gene1" in result_dict
    assert result_dict["gene1"] == {
        "sequence": "MKTSF",
        "offset": 0,
        "insertions": ["5:VLK"],
    }

    # Check that gene2 is represented as null (only sequence emptiness matters)
    assert "gene2" in result_dict
    assert result_dict["gene2"] is None

    # Test that the JSON schema validation still works
    try:
        json_output = read.to_silo_json()
        assert json_output is not None
        # Parse the JSON to verify structure
        import json

        parsed = json.loads(json_output)
        assert parsed["gene2"] is None
        assert "gene1" in parsed
    except Exception as e:
        pytest.fail(f"Schema validation failed for empty gene: {e}")
