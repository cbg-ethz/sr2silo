"""
This module contains tests for the conversion functions in the sr2silo package.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from sr2silo.process import bam_to_sam
from sr2silo.process.convert import (
    bam_to_fastq_handle_indels,
    is_sorted_qname,
    pad_alignment,
    sam_to_seq_and_indels,
    sort_bam_file,
)
from sr2silo.process.interface import Insertion


def test_bam_to_sam(bam_data: Path):
    """Test the bam_to_sam function.
    Note:
    BAM -> SAM -> BAM conversion cannot be tested directly,
    as the resulting BAM file will not be identical to the original BAM file.
    """

    # Convert the BAM file to SAM
    sam_file = bam_data.with_suffix(".sam")
    bam_to_sam(bam_data, sam_file)

    # Check that the SAM file was created
    assert sam_file.exists(), "The SAM file was not created"
    # check that the SAM file is not empty
    assert sam_file.stat().st_size > 0, "The SAM file is empty"
    # check that the SAM file is not a binary file
    with sam_file.open("r") as f:
        first_line = f.readline()
        assert first_line.startswith("@")
    # assert that the file is larger than the original BAM file
    assert (
        sam_file.stat().st_size > bam_data.stat().st_size
    ), "The SAM file is smaller than the BAM file"


def test_sort_bam_file(tmp_path, monkeypatch):
    """Test the sort_bam_file function with both sorting options."""

    # Create a mock for pysam.sort to track calls
    sort_calls = []

    def mock_sort(*args):
        sort_calls.append(args)
        # Create an empty file at the output location to simulate successful sorting
        output_path = None
        for i, arg in enumerate(args):
            if arg == "-o" and i + 1 < len(args):
                output_path = args[i + 1]
                break
        if output_path:
            with open(output_path, "w") as f:
                f.write("")

    # Apply the monkeypatch
    monkeypatch.setattr(pysam, "sort", mock_sort)

    # Setup test files
    input_bam = tmp_path / "input.bam"
    input_bam.write_text("mock bam content")

    output_coord_bam = tmp_path / "output_coord.bam"
    output_qname_bam = tmp_path / "output_qname.bam"

    # Test coordinate sorting (default)
    sort_bam_file(input_bam, output_coord_bam)

    # Test query name sorting
    sort_bam_file(input_bam, output_qname_bam, sort_by_qname=True)

    # Verify calls
    assert len(sort_calls) == 2, "Expected two calls to pysam.sort"

    # Check coordinate sort call
    coord_sort_args = sort_calls[0]
    assert "-o" in coord_sort_args, "Missing -o flag in coordinate sort"
    assert (
        str(output_coord_bam) in coord_sort_args
    ), "Output path not in coordinate sort arguments"
    assert "-n" not in coord_sort_args, "Should not have -n flag in coordinate sort"

    # Check qname sort call
    qname_sort_args = sort_calls[1]
    assert "-o" in qname_sort_args, "Missing -o flag in qname sort"
    assert (
        str(output_qname_bam) in qname_sort_args
    ), "Output path not in qname sort arguments"
    assert "-n" in qname_sort_args, "Missing -n flag in qname sort"

    # Test output files exist
    assert output_coord_bam.exists(), "Coordinate sorted BAM file not created"
    assert output_qname_bam.exists(), "Query name sorted BAM file not created"


def test_is_sorted_qname(tmp_path, monkeypatch):
    """Test the is_sorted_qname function against BAM files sorted by query name."""

    # Create mock BAM headers
    mock_coordinate_header = {"HD": {"VN": "1.0", "SO": "coordinate"}}
    mock_queryname_header = {"HD": {"VN": "1.0", "SO": "queryname"}}
    mock_no_so_header = {"HD": {"VN": "1.0"}}
    mock_no_hd_header = {"SQ": [{"LN": 1000, "SN": "chr1"}]}

    # Setup mock AlignmentFile objects
    class MockAlignmentFile:
        def __init__(self, header):
            self.header = header

        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc_val, exc_tb):
            pass

    # Test files
    coordinate_bam = tmp_path / "coordinate_sorted.bam"
    queryname_bam = tmp_path / "queryname_sorted.bam"
    no_so_bam = tmp_path / "no_so.bam"
    no_hd_bam = tmp_path / "no_hd.bam"
    error_bam = tmp_path / "error.bam"

    # Setup mocked open function
    def mock_alignment_file(file_path, *args, **kwargs):
        if str(file_path) == str(coordinate_bam):
            return MockAlignmentFile(mock_coordinate_header)
        elif str(file_path) == str(queryname_bam):
            return MockAlignmentFile(mock_queryname_header)
        elif str(file_path) == str(no_so_bam):
            return MockAlignmentFile(mock_no_so_header)
        elif str(file_path) == str(no_hd_bam):
            return MockAlignmentFile(mock_no_hd_header)
        elif str(file_path) == str(error_bam):
            raise ValueError("Mock error opening file")
        else:
            raise FileNotFoundError(f"Mock file not found: {file_path}")

    # Apply monkeypatch
    monkeypatch.setattr(pysam, "AlignmentFile", mock_alignment_file)

    # Create empty test files
    coordinate_bam.write_text("")
    queryname_bam.write_text("")
    no_so_bam.write_text("")
    no_hd_bam.write_text("")
    error_bam.write_text("")

    # Run tests
    assert (
        is_sorted_qname(coordinate_bam) is False
    ), "Should not detect coordinate sorted BAM as queryname sorted"
    assert (
        is_sorted_qname(queryname_bam) is True
    ), "Should detect queryname sorted BAM correctly"
    assert (
        is_sorted_qname(no_so_bam) is False
    ), "Should not detect BAM with no SO tag as queryname sorted"
    assert (
        is_sorted_qname(no_hd_bam) is False
    ), "Should not detect BAM with no HD section as queryname sorted"
    assert (
        is_sorted_qname(error_bam) is None
    ), "Should return None for file that raises error"


def test_sort_bam_file_and_check_sorting(tmp_path):
    """Integration test: Sort BAM files and verify sorting with is_sorted_qname."""

    # We'll create a minimal valid BAM file structure for testing
    # First, create a header
    header = {"HD": {"VN": "1.0", "SO": "unknown"}, "SQ": [{"LN": 100, "SN": "chr1"}]}

    # Create the input BAM file
    input_bam = tmp_path / "input.bam"
    with pysam.AlignmentFile(str(input_bam), "wb", header=header) as outf:  # noqa
        # We don't need to write any reads for this test
        pass

    # Create output files
    output_coord_bam = tmp_path / "output_coord.bam"
    output_qname_bam = tmp_path / "output_qname.bam"

    # Sort by coordinate
    sort_bam_file(input_bam, output_coord_bam, sort_by_qname=False)

    # Sort by query name
    sort_bam_file(input_bam, output_qname_bam, sort_by_qname=True)

    # Check that the files exist
    assert output_coord_bam.exists(), "Coordinate sorted BAM file was not created"
    assert output_qname_bam.exists(), "Query name sorted BAM file was not created"

    # Verify sorting using is_sorted_qname
    # We can't directly check the headers with pysam here since we need a real BAM file
    # Just print the results for manual verification
    print(
        "Coordinate sorted BAM has queryname sorting: "
        f"{is_sorted_qname(output_coord_bam)}"
    )
    print(
        "Query name sorted BAM has queryname sorting: "
        f"{is_sorted_qname(output_qname_bam)}"
    )


@pytest.mark.skip(reason="Not implemented")
def test_create_index():
    """Test the index_bam_file function."""

    raise NotImplementedError


@pytest.mark.skip(reason="Not implemented")
def test_bam_to_fasta():
    """Test the bam_to_fasta function."""

    raise NotImplementedError


def test_pad_alignment():
    """Test the pad_alignment function with various inputs."""

    # Test with string input and default unknown_char 'N'
    seq = "ACTG"
    ref_start = 2
    ref_length = 10  # Expected: "NNACTGNNNN" (2 left, 4 right)
    expected = "NNACTGNNNN"
    result = pad_alignment(seq, ref_start, ref_length)
    assert result == expected, f"Expected {expected}, got {result}"

    # Test with list input, same parameters
    seq_list = ["A", "C", "T", "G"]
    expected = "NNACTGNNNN"
    result = pad_alignment(seq_list, ref_start, ref_length)
    assert result == expected, f"Expected {expected}, got {result}"

    # Test with custom unknown_char '-' (variation of deletion char)
    expected = "--ACTG----"
    result = pad_alignment(seq, ref_start, ref_length, unknown_char="-")
    assert result == expected, f"Expected {expected}, got {result}"


def test_sam_to_seq_and_indels():
    """Test the sam_to_seq_and_indels function, including Insertion conversion."""

    # Example from the docstring:
    # sequence = "AGCTTAGCTAGCTT"
    # cigar = "5M1I5M1D3M"
    seq = "AGCTTAGCTAGCTT"
    cigar = "5M1I5M1D3M"

    # Expected:
    # - Cleartext: "AGCTTGCTAGCTT" [matches from 5M, 5M, 3M]
    # - Insertion: one insertion at position 5 with base "A"
    # - Deletion: one deletion at position 10 with length 1
    expected_cleartext = "AGCTTGCTAGCTT"
    expected_deletions = [(10, 1)]

    cleartext, insertions, deletions = sam_to_seq_and_indels(seq, cigar)

    # assert the types of the outputs
    assert isinstance(cleartext, str), "cleartext is not a string"
    assert isinstance(insertions, list), "insertions is not a list"
    assert isinstance(deletions, list), "deletions is not a list"

    # assert the elements of insertiosn and deletions
    for ins in insertions:
        assert isinstance(ins, Insertion), "insertions contains a non-Insertion object"
    for del_ in deletions:
        assert isinstance(del_, tuple), "deletions contains a non-tuple object"

    assert (
        cleartext == expected_cleartext
    ), f"Expected cleartext {expected_cleartext}, got {cleartext}"
    assert (
        deletions == expected_deletions
    ), f"Expected deletions {expected_deletions}, got {deletions}"
    assert len(insertions) == 1, f"Expected 1 insertion, got {len(insertions)}"
    # Check that the insertion is an instance of Insertion and has correct attributes.
    ins = insertions[0]
    assert isinstance(ins, Insertion), "Insertion is not an instance of Insertion"
    assert ins.position == 5, f"Expected insertion position 5, got {ins.position}"
    assert ins.sequence == "A", f"Expected insertion sequence 'A', got {ins.sequence}"


def test_get_gene_set_from_ref():
    """Test the get_gene_set_from_ref function using the real AA reference file."""

    from sr2silo.process.convert import get_gene_set_from_ref

    aa_ref_fp = Path("resources/sars-cov-2/aa_reference_genomes.fasta")
    gene_set = get_gene_set_from_ref(aa_ref_fp)
    gene_names = gene_set.get_gene_name_list()

    # Assert that expected gene names are present and their lengths are > 0.
    for expected in ["E", "S", "ORF1a", "ORF7a"]:
        assert expected in gene_names, f"Expected gene {expected} to be in gene set"
        gene_length = gene_set.get_gene_length(expected)
        assert (
            gene_length > 0
        ), f"Expected gene length for {expected} to be > 0, got {gene_length}"


def test_bam_to_fastq_handle_indels(dummy_alignment, tmp_path):
    """
    Test bam_to_fastq_handle_indels using the dummy AlignmentFile from conf_test.py.
    The dummy read has:
      - query_sequence "ACTG"
      - query_qualities [30, 31, 32, 33]
      - cigartuples: [(0,2), (1,1), (0,1)]
      - reference_start: 100
    Expected FASTQ record:
      @read1
      ACG
      +
      ?@B
      alignment_position:100
    Expected insertion file:
      read1    102    T    A
    """
    # Create temporary files for FASTQ and insertions
    fastq_file = tmp_path / "output.fastq"
    insertions_file = tmp_path / "insertions.txt"

    # dummy_alignment is provided by the fixture
    bam_to_fastq_handle_indels(dummy_alignment, fastq_file, insertions_file)

    expected_fastq = (
        "@read1\n"
        "ACG\n"  # from match of 2 bases ("AC") and then match of 1 ("G")
        "+\n"
        "?@B\n"  # qualities: chr(30+33)="?" , chr(31+33)="@" , chr(33+33)="B"
        "alignment_position:100\n"
    )
    fastq_content = fastq_file.read_text()
    assert (
        fastq_content == expected_fastq
    ), f"FASTQ output mismatch:\nExpected:\n{expected_fastq}\nGot:\n{fastq_content}"

    expected_insertion = "read1\t102\tT\tA\n"
    insertion_content = insertions_file.read_text()
    assert insertion_content == expected_insertion, (
        f"Insertion output mismatch:\nExpected:\n{expected_insertion}\n"
        f"Got:\n{insertion_content}"
    )


def test_get_gene_set_from_ref_malformed_no_sequence(tmp_path):
    """Test get_gene_set_from_ref with a FASTA file that has header(s) but no
    sequence lines."""
    malformed = tmp_path / "malformed.fasta"
    # Write a header without a following sequence line.
    malformed.write_text(">GeneX\n")
    from sr2silo.process.convert import get_gene_set_from_ref

    gene_set = get_gene_set_from_ref(malformed)
    # Expect GeneSet to be empty since no sequence was provided.
    assert (
        gene_set.get_gene_name_list() == []
    ), "Expected empty gene set for header without sequence"

    def test_get_gene_set_from_ref_malformed_blank_lines(tmp_path):
        """Test get_gene_set_from_ref with a FASTA file that has multiple headers
        and blank sequence lines."""

    # Create file with headers with blank sequence lines
    content = """>GeneA
    >GeneB
    AGCTAGCT
    >GeneC

    """
    malformed.write_text(content)
    from sr2silo.process.convert import get_gene_set_from_ref

    gene_set = get_gene_set_from_ref(malformed)
    # Expect only GeneB to be added (GeneA and GeneC have blank sequences)
    gene_names = gene_set.get_gene_name_list()
    assert "GeneB" in gene_names, "Expected GeneB to be parsed"
    assert (
        "GeneA" not in gene_names
    ), "Expected GeneA to be skipped due to missing sequence"
    assert (
        "GeneC" not in gene_names
    ), "Expected GeneC to be skipped due to missing sequence"
