"""Tests for the translation module."""

from __future__ import annotations

import copy
import json
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
    NucInsertion,
)

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def make_reference_nextclade(reference_fasta_fp: Path) -> None:
    """This function makes a reference file out of sr2silos reference file,
    to be used with nextclade.
    """

    # get reference of sr2silo - resources/sars-cov-2/nuc_reference_genomes.fasta
    sr2silo_reference = (
        Path("resources/sars-cov-2/nuc_reference_genomes.fasta").read_text().split("\n")
    )
    sr2silo_reference = "".join(
        [line for line in sr2silo_reference if not line.startswith(">")]
    )
    # write to file with one header > and the seqeunce in one line
    with open(reference_fasta_fp, "w") as f:
        f.write(">reference\n")
        f.write(sr2silo_reference)


def translate_align_nextclade_ref(
    input_files: List[Path], result_dir: Path, reference_fasta_fp: Path
) -> None:
    """Nextclades alignment of reads in fasta format and translation and
    alignment of the reads in amino acid space.

    Note: This function is a wrapper around the nextclade command line tool.

    This implementation is meant only of orthogonal testing.

    Args:
        input_file (str): The path to the input file.
                          the nucleotide sequences in fasta format.
        result_dir (str): The path to the directory to save the results.

        reference_fasta_fp (str): The path to the reference fasta file.
    """

    with tempfile.TemporaryDirectory() as temp_dir:
        logging.debug(f"temp_dir: {temp_dir}")
        # first get the test dataset from the gff3 file
        command = [
            "nextclade",
            "dataset",
            "get",
            "--name",
            "nextstrain/sars-cov-2/XBB",
            "--output-dir",
            temp_dir,
        ]
        logging.debug(f"Running command: {command}")
        subprocess.run(command, check=True)

        # then replace the reference.fasta in the temp_dir
        #  with the reference.fasta from the input file
        command = ["cp", reference_fasta_fp, f"{temp_dir}/reference.fasta"]
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


def test_translate_align_nextclade(temp_dir):
    """Test that the translate_align using nextclade executes
    without error.
    """
    translate_align_nextclade(
        [Path("tests/data/merged_expected.fasta")],
        Path(temp_dir),
        "nextstrain/sars-cov-2/XBB",
    )


def test_validate_nextclade_reference():

    nextclade_reference = "nextstrain/sars-cov-2/XBB"

    temp_dir = tempfile.mkdtemp()
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

    # read in the reference.fasta file to str, skip lines with > as they are headers
    reference_fasta = Path(temp_dir) / "reference.fasta"
    reference_fasta_str = reference_fasta.read_text().split("\n")
    reference_fasta_str = "".join(
        [line for line in reference_fasta_str if not line.startswith(">")]
    )

    # get reference of sr2silo - resources/sars-cov-2/nuc_reference_genomes.fasta
    sr2silo_reference = (
        Path("resources/sars-cov-2/nuc_reference_genomes.fasta").read_text().split("\n")
    )
    sr2silo_reference = "".join(
        [line for line in sr2silo_reference if not line.startswith(">")]
    )
    # write to file with one header > and the seqeunce in one line
    with open("reference.fasta", "w") as f:
        f.write(">reference\n")
        f.write(reference_fasta_str)

    # check that the reference.fasta file is the same as the reference.fasta file in sr2silo
    # these are long strings please show me the difference explicitly
    assert (
        reference_fasta_str == sr2silo_reference
    ), "The reference.fasta file is not the same as the reference.fasta file in sr2silo"


def test_parse_translate_align_orth_nextclade(bam_and_fasta_raw_data):
    """Test the translate_align() orthogonally using nextclade."""

    # get the test data
    fasta_raw_data = bam_and_fasta_raw_data[1]
    bam_path = bam_and_fasta_raw_data[0]

    # with tempfile.TemporaryDirectory() as tmpdirname:

    tempdirname = Path("nextclade_output")
    tempdirname.mkdir(parents=True, exist_ok=True)

    output_dir = Path(tempdirname) / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    nuc_reference_genomes_fp = output_dir / "nuc_reference_genomes.fasta"
    make_reference_nextclade(nuc_reference_genomes_fp)

    translate_align_nextclade_ref(
        [fasta_raw_data], output_dir, nuc_reference_genomes_fp
    )

    ### Parse the Nextclade file to AlignedReads
    ## Read in Insertions from nextclade.ndjson, "insertions" for NucInsertions and "aaInsertions" for AAInsertions
    # for each line in the json file, load the JSON and get insertions

    nuc_insertions_store = {}
    for line in (output_dir / "nextclade.ndjson").read_text().split("\n"):
        if not line:
            continue
        json_data = json.loads(line)
        read_id = json_data["seqName"]
        nuc_insertions = []
        # now json_data["insertions"] is a list as  [{'pos': 1929, 'ins': 'CTA'}]  read in these objects as NucInsertions
        for nuc_ins in json_data["insertions"]:
            nuc_insertions.append(NucInsertion(nuc_ins["pos"], nuc_ins["ins"]))
            nuc_insertions_store[read_id] = nuc_insertions
            # TODO: add ammino acid insertions test from "aaInsertions" if test data contains them

    ## from nextclade.aligned.fasta --> alignedNucSequence  (need to adjust Padding Char from - to N) - also adjust deletion char to from XXX to -
    alignedSequenceStore = {}
    with open(output_dir / "nextclade.aligned.fasta") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                read_id = line[1:]
                alignedSequenceStore[read_id] = ""
            else:
                alignedSequenceStore[read_id] += line  # type: ignore

    # replace padding char from - to N
    for read_id, aligned_sequence in alignedSequenceStore.items():
        alignedSequenceStore[read_id] = aligned_sequence.replace("-", "N")

    # TODO: replace deletion char from XXX to - (deletion char nextclade unknown)

    ## from nextclade.cds_translation.<<GeneName>>.fasta --> alignedAASequence (need to adjust Padding Char from - to X) - also adjust deletion char to from XXX to -
    alignedAASequenceStore = {}

    # get all gene names from the nextclade.cds_translation files
    gene_names = [
        file.stem.split(".")[2]
        for file in (output_dir).glob("nextclade.cds_translation.*.fasta")
    ]
    # make List[GeneName]
    gene_names = [GeneName(gene_name) for gene_name in gene_names]

    # get all file names in dir that fit the pattern nextclade.cds_translation.<<GeneName>>.fasta
    for file in (output_dir).glob("nextclade.cds_translation.*.fasta"):
        gene_name = file.stem.split(".")[2]
        with open(file) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    read_id = line[1:]
                    # if read id not in store make new AASequenceSet
                    if read_id not in alignedAASequenceStore:
                        alignedAASequenceStore[read_id] = AASequenceSet(gene_names)
                else:
                    # convert padding char from - to N and also the undertermined char X to N
                    alignedAASequenceStore[read_id].set_sequence(gene_name, line.replace("-", "N").replace("X", "N"))  # type: ignore

                    #  TODO: replace deletion char from XXX to -   (deletion char nextclade unknown)

    # make AlignedReads from the data
    aligned_reads = {}
    for read_id in alignedSequenceStore.keys():
        aligned_read = AlignedRead(
            read_id=read_id,
            unaligned_nucleotide_sequences="",
            aligned_nucleotide_sequences=alignedSequenceStore[read_id],
            nucleotide_insertions=nuc_insertions_store.get(read_id, []),
            amino_acid_insertions=AAInsertionSet(gene_names),
            aligned_amino_acid_sequences=alignedAASequenceStore.get(
                read_id, AASequenceSet(gene_names)
            ),
        )
        aligned_reads[read_id] = aligned_read

    ### Run the parse_translate_align function
    aligned_read_actual = translate_align.parse_translate_align(
        nuc_reference_fp=Path(
            "resources/sars-cov-2/nuc_reference_genomes.fasta"
        ),  # TODO: replace with fixture
        aa_reference_fp=Path(
            "resources/sars-cov-2/aa_reference_genomes.fasta"
        ),  # TODO: replace with fixture
        nuc_alignment_fp=bam_path,
    )

    # check that the read_ids are the same
    assert aligned_reads.keys() == aligned_read_actual.keys()

    # write all reads expected and actual to a file
    with open("aligned_reads.ndjson", "w") as f:
        for read_id in aligned_reads.keys():
            f.write(str(aligned_reads[read_id]) + "\n")

    with open("aligned_reads_actual.ndjson", "w") as f:
        for read_id in aligned_reads.keys():
            f.write(str(aligned_read_actual[read_id]) + "\n")

    # for a given read_id check that the aligned sequences are the same
    for read_id in aligned_reads.keys():
        for field in [
            # "aligned_nucleotide_sequences", # skip ad V-Pipe does the alignment
            # "aligned_amino_acid_sequences",
            "amino_acid_insertions",
        ]:
            expected_value = getattr(aligned_reads[read_id], field)
            actual_value = getattr(aligned_read_actual[read_id], field)

            # Use direct comparison
            assert expected_value == actual_value, (
                f"Mismatch for read_id '{read_id}' in field '{field}': "
                f"expected {expected_value}, got {actual_value}"
            )

        # Special handling for nucleotide insertions - compare by string representation
        expected_nuc_insertions = aligned_reads[read_id].nucleotide_insertions
        actual_nuc_insertions = getattr(
            aligned_read_actual[read_id], "nucleotide_insertions"
        )

        # Convert both lists to sorted lists of string representations for stable comparison
        expected_nuc_insertions_str = sorted(
            str(ins) for ins in expected_nuc_insertions
        )
        actual_nuc_insertions_str = sorted(str(ins) for ins in actual_nuc_insertions)

        assert expected_nuc_insertions_str == actual_nuc_insertions_str, (
            f"Mismatch for read_id '{read_id}' in nucleotide_insertions: "
            f"expected {expected_nuc_insertions_str}, got {actual_nuc_insertions_str}"
        )


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
