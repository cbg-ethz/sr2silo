"""Script to convert SAM to json gene, sequence insertions
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple

from sr2silo.process import pad_alignment

def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """Parse the CIGAR string into a list of tuples."""
    import re

    return [
        (int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


def process_sequence(
    seq: str, cigar: str
) -> Tuple[str, List[Tuple[int, str]], List[Tuple[int, int]]]:
    """
    Processes a SAM file-style sequence and a CIGAR string to return the
    cleartext sequence, along with detailed information about insertions
    and deletions.

    Args:
        seq (str): The sequence string from the SAM file, representing the read.
        cigar (str): The CIGAR string that describes how the sequence aligns
                    to a reference.

    Returns:
        tuple: A tuple containing:
            - cleartext_sequence (str): The sequence aligned to the reference,
                                         excluding insertions and deletions.
            - insertions (list of tuples): A list of tuples, each containing:
                - position (int): The position in the reference where the
                                  insertion occurs.
                - inserted_sequence (str): The sequence that is inserted at
                                           the given position.
            - deletions (list of tuples): A list of tuples, each containing:
                - position (int): The position in the reference where the
                                  deletion starts.
                - length (int): The length of the deletion.

    Example:
        sequence = "AGCTTAGCTAGCTT"
        cigar = "5M1I5M1D3M"
        cleartext, insertions, deletions = process_sequence(sequence, cigar)

        # Output:
        # Cleartext Sequence: AGCTTAGCTAGC
        # Insertions: [(5, 'A')]
        # Deletions: [(11, 1)]

    Notes:
        - The function assumes that the input sequence and CIGAR string are
          valid and correctly formatted.
        - The CIGAR operations handled include:
            - 'M', '=', 'X': Match or mismatch (aligned to the reference).
            - 'I': Insertion to the reference.
            - 'D': Deletion from the reference.
            - 'N': Skipped region from the reference.
            - 'S': Soft clipping (clipped sequences present in SEQ).
            - 'H': Hard clipping (clipped sequences NOT present in SEQ).
            - 'P': Padding (silent deletion from padded reference).
    """
    parsed_cigar = parse_cigar(cigar)
    cleartext_sequence = []
    insertions = []
    deletions = []

    seq_index = 0
    ref_position = 0

    for length, op in parsed_cigar:
        if op == "M" or op == "=" or op == "X":  # Match or mismatch
            cleartext_sequence.append(seq[seq_index : seq_index + length])
            seq_index += length
            ref_position += length
        elif op == "I":  # Insertion to the reference
            insertions.append((ref_position, seq[seq_index : seq_index + length]))
            seq_index += length
        elif op == "D":  # Deletion from the reference
            deletions.append((ref_position, length))
            ref_position += length
        elif op == "N":  # Skipped region from the reference
            ref_position += length
        elif op == "S":  # Soft clipping (clipped sequences present in SEQ)
            seq_index += length
        elif op == "H":  # Hard clipping (clipped sequences NOT present in SEQ)
            pass
        elif op == "P":  # Padding (silent deletion from padded reference)
            pass

    return "".join(cleartext_sequence), insertions, deletions


class AlignedRead:
    """Class to represent an aligned read."""

    def __init__(
        self,
        read_id: str,
        unaligned_nucleotide_sequences: str,
        aligned_nucleotide_sequences: str,
        nucleotide_insertions: List[Tuple[int, str]],
        amino_acid_insertions: List[Tuple[int, str]],
        aligned_amino_acid_sequences: Dict[str, str],
    ):
        self.read_id = read_id
        self.unaligned_nucleotide_sequences = unaligned_nucleotide_sequences
        self.aligned_nucleotide_sequences = aligned_nucleotide_sequences
        self.nucleotide_insertions = nucleotide_insertions
        self.amino_acid_insertions = amino_acid_insertions
        self.aligned_amino_acid_sequences = aligned_amino_acid_sequences

    def to_dict(self) -> Dict[str, any]: # noqa
        return {
            "read_id": self.read_id,
            "unaligned_nucleotide_sequences": self.unaligned_nucleotide_sequences,
            "aligned_nucleotide_sequences": self.aligned_nucleotide_sequences,
            "nucleotide_insertions": self.nucleotide_insertions,
            "amino_acid_insertions": self.amino_acid_insertions,
            "aligned_amino_acid_sequences": self.aligned_amino_acid_sequences,
        }


class Gene:
    def __init__(self, gene_name: str, gene_length: int):
        self.gene_name = gene_name
        self.gene_length = gene_length

    def to_dict(self) -> Dict[str, int]:
        return {
            "gene_name": self.gene_name,
            "gene_length": self.gene_length,
        }


def get_genes_and_lengths_from_ref(reference_fp: Path) -> Dict[str, Gene]:
    """Load the gene ref fasta and get all the gene names."""
    genes = dict()

    with open(reference_fp, "r") as f:
        awaiting_next_line = False
        for line in f:
            if line.startswith(">"):
                gene = line[1:].strip()
                awaiting_next_line = True
            elif awaiting_next_line:
                reference_length = len(line.strip())
                genes[gene] = Gene(gene, reference_length)
                awaiting_next_line = False
            else:
                continue

    return genes


# Define the input and reference file paths
INPUT_FILE = "diamond_blastx.sam"
REFERENCE_FILE = "../resources/sars-cov-2/reference_genomes.fasta"



gene_dict = get_genes_and_lengths_from_ref(REFERENCE_FILE)


with open(INPUT_FILE, "r") as f:
    for line in f:
        # Skip header lines
        if line.startswith("@"):
            continue
        # Split the line into fields
        fields = line.strip().split("\t")
        read_id = fields[0]
        gene_name = fields[2]
        pos = int(fields[3])
        cigar = fields[5]
        seq = fields[9]

        aa_aligned, aa_insertions, aa_deletions = process_sequence(seq, cigar)

        # Make a dict to hold the aligned amino acid sequences
        aligned_amino_acid_sequences = {}
        # Write a null for all gene names
        for gene in get_genes_and_lengths_from_ref(REFERENCE_FILE).keys():
            aligned_amino_acid_sequences[gene] = None


        # pad the alignment
        padded_alignment = pad_alignment(
            seq, pos, gene_dict[gene_name].gene_length
        )

        print(padded_alignment)

        # Print the results
        print(f"Read ID: {read_id}")
        print(f"Gene Name: {gene_name}")
        print(f"AA Aligned: {aa_aligned}")
        print(f"AA Insertions: {aa_insertions}")
        print(f"AA Deletions: {aa_deletions}")

