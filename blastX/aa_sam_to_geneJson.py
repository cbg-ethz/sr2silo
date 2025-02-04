"""Script to convert SAM to json gene, sequence insertions
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

from sr2silo.process import pad_alignment


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """Parse the CIGAR string into a list of tuples."""
    import re

    return [
        (int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


class AAInsertion:
    def __init__(self, position: int, sequence: str):
        self.position = position
        self.sequence = sequence

    def __str__(self) -> str:
        return f"{self.position} : {self.sequence}"


def process_sequence(
    seq: str, cigar: str
) -> Tuple[str, List[AAInsertion], List[Tuple[int, int]]]:
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

    # convert insertions to AAInsertion objects
    insertions = [AAInsertion(position, sequence) for position, sequence in insertions]

    return "".join(cleartext_sequence), insertions, deletions


class AlignedRead:
    """Class to represent an aligned read."""

    def __init__(
        self,
        read_id: str,
        unaligned_nucleotide_sequences: str,
        aligned_nucleotide_sequences: str,
        nucleotide_insertions: any,
        amino_acid_insertions: AAInsertionSet,
        aligned_amino_acid_sequences: Dict[str, str],
    ):
        self.read_id = read_id
        self.unaligned_nucleotide_sequences = unaligned_nucleotide_sequences
        self.aligned_nucleotide_sequences = aligned_nucleotide_sequences
        self.nucleotide_insertions = nucleotide_insertions
        self.amino_acid_insertions = amino_acid_insertions
        self.aligned_amino_acid_sequences = aligned_amino_acid_sequences

    def to_dict(self) -> Dict[str, any]:  # noqa
        return {
            "read_id": self.read_id,
            "unaligned_nucleotide_sequences": self.unaligned_nucleotide_sequences,
            "aligned_nucleotide_sequences": self.aligned_nucleotide_sequences,
            "nucleotide_insertions": self.nucleotide_insertions,
            "amino_acid_insertions": self.amino_acid_insertions,
            "aligned_amino_acid_sequences": self.aligned_amino_acid_sequences,
        }

    def get_amino_acid_insertions(self) -> Dict[str, List[AAInsertion]]:
        print(self.amino_acid_insertions)
        return self.amino_acid_insertions

    def to_json(self) -> str:
        json_representation = {
            "nucleotideInsertions": {
                "main": self.nucleotide_insertions,
            },
            "aminoAcidInsertions": {
                k: [i.__str__() for i in v]
                for k, v in self.get_amino_acid_insertions().items()
            },
            "alignedNucleotideSequences": {
                "main": self.aligned_nucleotide_sequences,
            },
            "unalignedNucleotideSequences": {
                "main": self.unaligned_nucleotide_sequences,
            },
            "alignedAminoAcidSequences": {
                "main": self.aligned_amino_acid_sequences,
            },
        }
        return json.dumps(json_representation, indent=4)


class Gene:
    def __init__(self, gene_name: str, gene_length: int):
        self.gene_name = gene_name
        self.gene_length = gene_length

    def to_dict(self) -> Dict[str, int]:
        return {
            "gene_name": self.gene_name,
            "gene_length": self.gene_length,
        }


class AAInsertionSet:
    """Class to represent the set of amino acid insertions for a full set of genes."""

    def __init__(
        self, gene_dict: Dict[Gene], aa_insertions: List[AAInsertion], gene_name: str
    ):
        aa_insertions_set = {}

        for gene in gene_dict.keys():
            aa_insertions_set[gene] = []

        aa_insertions_set[gene_name] = aa_insertions

        self.aa_insertions_set = aa_insertions_set

    def to_dict(self) -> Dict[str, List[AAInsertion]]:
        return self.aa_insertions_set

    def items(self):
        return self.aa_insertions_set.items()


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


reads: List[AlignedRead] = []

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

        aa_insertions = [
            AAInsertion(34, "A"),
            AAInsertion(56, "B"),
            AAInsertion(78, "C"),
        ]
        # convert aa_insertions to dict of all gene names and add insertions to the correct gene

        aa_insertions_fmt = {}
        for gene in gene_dict.keys():
            aa_insertions_fmt[gene] = []
        aa_insertions_fmt[gene_name] = aa_insertions

        aa_insertions = AAInsertionSet(gene_dict, aa_insertions, gene_name)

        # Make a dict to hold the aligned amino acid sequences
        aligned_amino_acid_sequences = {}
        # Write a null for all gene names
        for gene in get_genes_and_lengths_from_ref(REFERENCE_FILE).keys():
            aligned_amino_acid_sequences[gene] = None

        # pad the alignment
        padded_alignment = pad_alignment(seq, pos, gene_dict[gene_name].gene_length)

        reads.append(
            AlignedRead(
                read_id=read_id,
                unaligned_nucleotide_sequences="null",
                aligned_nucleotide_sequences="null",
                nucleotide_insertions=[],
                amino_acid_insertions=aa_insertions,
                aligned_amino_acid_sequences=padded_alignment,
            )
        )

print(reads[0].to_json())
