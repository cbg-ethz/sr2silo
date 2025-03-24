"""Tests for the translation module."""

from __future__ import annotations

import copy
import logging
import tempfile
from pathlib import Path

import pytest

import sr2silo.process.translate_align as translate_align
from sr2silo.process import convert
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


def test_translate():
    """Test the translation function."""

    with tempfile.TemporaryDirectory() as tmpdirname:
        output_dir = Path(tmpdirname) / "output"
        output_dir.mkdir(parents=True, exist_ok=True)

        translate_align.translate_nextclade(
            [Path("tests/data/merged_expected.fasta")],
            output_dir,
            "nextstrain/sars-cov-2/XBB",
        )
    assert True


# Nota bene: The output tested against is not validated, yet a failure in test
# notes a change of output here.
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
            print(f"act_attr: {act_attr}")
            print(f"exp_attr: {exp_attr}")
            assert str(act_attr) == str(exp_attr), (
                f"Mismatch at read {i} for attribute {attr}: "
                f"Expected {exp_attr} but got {act_attr}"
            )

def test_parse_translate_align_synth(micro_bam_fp, micro_reference_fp, micro_aa_reference_fp):
    """Test the parse_translate_align function with synthetic data."""

    aligned_reads = translate_align.parse_translate_align(
        nuc_reference_fp = micro_reference_fp,
        aa_reference_fp = micro_aa_reference_fp,
        nuc_alignment_fp = micro_bam_fp,
    )




@pytest.mark.skip(reason="Not implemented")
def test_read_in_AligendReads_nuc_seq():
    """Test the make_read_with_nuc_seq function."""
    raise NotImplementedError


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


@pytest.mark.skip(reason="Not implemented")
def test_read_in_AligendReads_aa_ins():
    """Test the read_in_AlignedReads_aa_ins function."""
    raise NotImplementedError


@pytest.mark.parametrize(
    "read_id, gene_name_str, pos, cigar, expected_seq, expected_insertions",
    [
        (
            "read1",
            "GeneX",
            10,
            "5M1I5M1D3M2I2M",
            "CLEARINS",
            [AAInsertion(15, "INS"), AAInsertion(18, "IN")],
        ),
        ("read2", "GeneY", 5, "3M2I4M", "SEQ", [AAInsertion(7, "IN")]),
    ],
)
def test_enrich_read_with_aa_seq_multiple(
    tmp_path,
    monkeypatch,
    read_id,
    gene_name_str,
    pos,
    cigar,
    expected_seq,
    expected_insertions,
):
    """Test enrich_read_with_aa_seq with multiple insertions."""
    # Create a dummy GeneSet with one gene of variable length.
    gene_name = GeneName(gene_name_str)
    dummy_gene = Gene(gene_name, 50)
    dummy_gene_set = GeneSet([dummy_gene])

    # Prepare a dummy AlignedRead object; other attributes can be empty.
    dummy_aligned_read = AlignedRead(
        read_id=read_id,
        unaligned_nucleotide_sequences="",
        aligned_nucleotide_sequences="",
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([gene_name_str]),
        aligned_amino_acid_sequences=AASequenceSet([gene_name_str]),
    )
    aligned_reads = {read_id: dummy_aligned_read}

    # Create a temporary FASTA-like file with tab-separated fields.
    content = f"{read_id}\tX\t{gene_name_str}\t{pos}\tY\t{cigar}\tZ\tW\tV\tABCDEFGHIJ\n"
    dummy_file = tmp_path / "dummy_aa_alignment.txt"
    dummy_file.write_text(content)

    # Patch sam_to_seq_and_indels to return fixed values:
    def fake_sam_to_seq_and_indels(seq, cigar):
        """Return fixed values for testing."""
        return (expected_seq, expected_insertions, [(30, 1)])

    monkeypatch.setattr(
        "sr2silo.process.convert.sam_to_seq_and_indels", fake_sam_to_seq_and_indels
    )

    # Call the function under test.
    updated_reads = translate_align.enrich_read_with_aa_seq(
        aligned_reads, dummy_file, dummy_gene_set
    )

    # Expected padded sequence for AA alignment:
    expected_padded = convert.pad_alignment(expected_seq, pos, 50)
    # Expected AA insertions
    expected_ins_set = AAInsertionSet([gene_name])
    expected_ins_set.set_insertions_for_gene(gene_name_str, expected_insertions)

    # Retrieve the AA insertions for gene from updated_reads
    read = updated_reads.get(read_id)
    if read is not None:
        aa_insertion_set = read.get_amino_acid_insertions()
    else:
        raise ValueError(f"Read {read_id} not found in updated_reads")
    assert (
        str(gene_name_str) in aa_insertion_set.to_dict()
    ), f"{gene_name_str} not present in AA insertion set"

    # Check that insertions equal a list containing expected_insertions
    # as string representation.
    assert str(expected_ins_set) == str(
        aa_insertion_set
    ), f"Expected insertion {str(expected_ins_set)} not found; got {aa_insertion_set}"

    # Retrieve AA sequences and check they were set correctly.
    aa_seq_dict = updated_reads[read_id].aligned_amino_acid_sequences.to_dict()
    expected_padded_str = expected_padded
    assert aa_seq_dict.get(gene_name_str) == expected_padded_str, (
        f"Expected padded sequence {expected_padded_str}, "
        f"got {aa_seq_dict.get(gene_name_str)}"
    )


def test_curry_read_with_metadata(tmp_path):
    """Test the curry_read_with_metadata function."""
    # Create a temporary metadata JSON file
    metadata = {
        "read_id": "test_read",
        "sample_id": "sample123",
        "batch_id": "batch456",
        "sampling_date": "2023-01-01",
        "sequencing_date": "2023-01-15",
        "location_name": "Test Location",
        "read_length": "150",
        "primer_protocol": "v1",
        "location_code": "TL",
        "flow_cell_serial_number": "FC123",
        "sequencing_well_position": "A1",
        "primer_protocol_name": "Test Protocol",
        "nextclade_reference": "ref123",
    }

    metadata_file = tmp_path / "test_metadata.json"
    with open(metadata_file, "w") as f:
        import json

        json.dump(metadata, f)

    # Get the enrich function using curry_read_with_metadata
    enrich_read = translate_align.curry_read_with_metadata(metadata_file)

    # Create a dummy AlignedRead
    dummy_read = AlignedRead(
        read_id="read1",
        unaligned_nucleotide_sequences="ACGT",
        aligned_nucleotide_sequences="ACGT",
        nucleotide_insertions=[],
        amino_acid_insertions=AAInsertionSet([]),
        aligned_amino_acid_sequences=AASequenceSet([]),
    )

    # Enrich the read with metadata
    enriched_read = enrich_read(dummy_read)

    # Verify metadata was attached correctly
    assert enriched_read.metadata == metadata

    # Test error cases
    with pytest.raises(FileNotFoundError):
        translate_align.curry_read_with_metadata(tmp_path / "nonexistent_file.json")

    # Test invalid JSON
    invalid_json_file = tmp_path / "invalid.json"
    with open(invalid_json_file, "w") as f:
        f.write("This is not valid JSON")

    with pytest.raises(json.JSONDecodeError):
        translate_align.curry_read_with_metadata(invalid_json_file)

    # Test empty metadata
    empty_metadata_file = tmp_path / "empty.json"
    with open(empty_metadata_file, "w") as f:
        json.dump({}, f)

    with pytest.raises(ValueError, match="No metadata found in the file"):
        translate_align.curry_read_with_metadata(empty_metadata_file)
