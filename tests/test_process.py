"""This is a sample python file for testing functions from the source code."""

from __future__ import annotations

from sr2silo.process import parse_cigar
from sr2silo.process import pair_normalize_reads


def test_parse_cigar(sam_data):
    """Test the parse_cigar function."""

    cigars = []
    for line in sam_data.splitlines():
        if line.startswith("@"):
            continue
        fields = line.strip().split("\t")
        CIGAR = parse_cigar(fields[5])  # CIGAR string
        cigars.append(CIGAR)

    # validate that the CIGAR strings are parsed correctly
    # the output is a list of lists of tuples of first s capital letter
    # and second s integer
    expected_cigars = [
        [("S", 31), ("M", 220)],
        [("S", 22), ("M", 220)],
        [("M", 225), ("S", 26)],
        [("M", 223), ("S", 28)],
        [("S", 23), ("M", 224)],
        [("M", 220), ("S", 29)],
        [("S", 23), ("M", 228)],
        [("M", 226), ("S", 25)],
        [("S", 33), ("M", 218)],
        [("S", 33), ("M", 218)],
        [("M", 218), ("S", 33)],
        [("M", 218), ("S", 33)],
        [("S", 31), ("M", 217)],
        [("M", 221), ("S", 30)],
        [("S", 27), ("M", 224)],
        [("M", 221), ("S", 30)],
        [("S", 26), ("M", 225)],
        [("S", 25), ("M", 225)],
        [("M", 227), ("S", 24)],
        [("M", 213), ("S", 24)],
        [("S", 24), ("M", 227)],
        [("M", 229), ("S", 22)],
        [("S", 29), ("M", 222)],
        [("S", 29), ("M", 222)],
        [("M", 227), ("S", 24)],
        [("M", 220), ("S", 24)],
        [("S", 27), ("M", 224)],
        [("M", 209), ("S", 30)],
        [("S", 29), ("M", 213)],
        [("S", 29), ("M", 222)],
        [("S", 26), ("M", 225)],
        [("M", 224), ("S", 27)],
        [("M", 224), ("S", 27)],
        [("M", 224), ("S", 27)],
        [("S", 26), ("M", 225)],
        [("M", 220), ("S", 27)],
        [("S", 24), ("M", 219)],
        [("S", 24), ("M", 227)],
        [("S", 11), ("M", 219)],
        [("M", 226), ("S", 25)],
        [("M", 190), ("S", 24)],
        [("S", 34), ("M", 217)],
        [("M", 225), ("S", 26)],
        [("S", 26), ("M", 225)],
        [("M", 218), ("S", 33)],
        [("S", 27), ("M", 220)],
        [("M", 199), ("S", 30)],
        [("S", 33), ("M", 4), ("D", 3), ("M", 214)],
        [("M", 246), ("S", 5)],
        [("S", 26), ("M", 225)],
        [("S", 26), ("M", 191)],
        [("M", 221), ("S", 30)],
        [("M", 215), ("S", 30)],
        [("S", 30), ("M", 205)],
        [("M", 211), ("S", 33)],
        [("S", 20), ("M", 231)],
        [("S", 20), ("M", 231)],
        [("M", 227), ("S", 24)],
        [("M", 217), ("S", 24)],
        [("S", 24), ("M", 227)],
        [("M", 185), ("S", 21)],
        [("S", 24), ("M", 189)],
        [("M", 134), ("D", 26), ("M", 81), ("S", 33)],
    ]
    assert cigars == expected_cigars, f"Expected {expected_cigars}, but got {cigars}"


def test_pair_normalize_reads(sam_data, temp_dir):
    """Test the pair_normalize_reads function."""

    out_fasta = temp_dir / "out_merged.fasta"
    out_insertions = temp_dir / "out_insertions.txt"

    # Pipe the SAM data to the function as stdin
    pair_normalize_reads(sam_data, out_fasta, out_insertions)

    # Check the output files for expected content
    with open(out_fasta, "r") as fasta_file:
        fasta_content = fasta_file.read()
    with open(out_insertions, "r") as insertions_file:
        insertions_content = insertions_file.read()

    # Load the actual expected content from the files
    with open("tests/data/merged_expected.fasta", "r") as expected_fasta_file:
        expected_fasta_content = expected_fasta_file.read()
    with open("tests/data/nuc_insertions_expected.txt", "r") as expected_insertions_file:
        expected_insertions_content = expected_insertions_file.read()

    assert fasta_content == expected_fasta_content, f"Expected {expected_fasta_content}, but got {fasta_content}"
    assert (
        insertions_content == expected_insertions_content
    ), f"Expected {expected_insertions_content}, but got {insertions_content}"
