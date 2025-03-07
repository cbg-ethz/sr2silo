"""Tests for the translation module."""

from __future__ import annotations

import copy
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import List

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


def translate_align_nextclade(
    input_files: List[Path], result_dir: Path, nextclade_reference: str
) -> None:
    """Nextclades alignment of reads in fasta format and translation and
    alignment of the reads in amino acid space.

    Note: This function is a wrapper around the nextclade command line tool.

    This implementation is meant only of orthogonal testing.

    Args:
        input_file (str): The path to the input file.
                          the nucleotide sequences in fasta format.
        result_dir (str): The path to the directory to save the results.
        nextclade_reference (str): The path to the nextclade reference.
                                    e.g. nextstrain/sars-cov-2/XBB
                                    see `nextclade dataset list`
    """

    with tempfile.TemporaryDirectory() as temp_dir:
        logging.debug(f"temp_dir: {temp_dir}")
        # first get the test dataset from the gff3 file
        command = [
            "nextclade",
            "dataset",
            "get",
            "--name",
            f"{nextclade_reference}",
            "--output-dir",
            temp_dir,
        ]
        logging.debug(f"Running command: {command}")
        subprocess.run(command, check=True)

        for input_file in input_files:
            logging.info(f"Translating {input_file}")

            # then replace the sequences.fasta in the temp_dir
            #  with the sequences.fasta from the input file
            command = ["cp", input_file, f"{temp_dir}/sequences.fasta"]
            logging.debug(f"Running command: {command}")
            subprocess.run(command, check=True)

            # then run the nextclade run command
            command = [
                "nextclade",
                "run",
                "--input-dataset",
                temp_dir,
                f"--output-all={result_dir}/",
                f"{temp_dir}/sequences.fasta",
            ]

            logging.debug(f"Running nextclade: {command}")

            try:
                result = subprocess.run(
                    command, check=True, capture_output=True, text=True
                )
                logging.debug(result.stdout)
                logging.debug(result.stderr)
            except subprocess.CalledProcessError as e:
                logging.error(f"nextclade failed with exit code {e.returncode}")
                logging.error(e.stderr)
                raise

            # move the results to the result_dir
            result_path = result_dir / input_file.stem
            command = ["mv", f"{temp_dir}/results", str(result_path)]


def test_translate_align_nextclade():
    """Test that the translate_align using nextclade executes
    without error.
    """

    with tempfile.TemporaryDirectory() as tmpdir:

        translate_align_nextclade(
            [Path("tests/data/merged_expected.fasta")],
            Path(tmpdir),
            "nextstrain/sars-cov-2/XBB",
        )


def test_parse_translate_align_orth_nextclade():
    """Test the translate_align() orthogonally using nextclade."""

    # with tempfile.TemporaryDirectory() as tmpdirname:

    tempdirname = Path("nextclade_output")
    tempdirname.mkdir(parents=True, exist_ok=True)

    output_dir = Path(tempdirname) / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Test Data
    # TODO: get the test data with a lot of insertions

    translate_align_nextclade(
        [Path("tests/data/merged_expected.fasta")],
        output_dir,
        "nextstrain/sars-cov-2/XBB",
    )

    # Parse the Nextclade file to AlignedReads

    # Run the parse_translate_align function

    # The do attribute comparison


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


# TODO: Implement the following tests
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


# TODO: Implement the following tests
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
