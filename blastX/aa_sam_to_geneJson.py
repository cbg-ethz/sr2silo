"""Script to convert SAM to json gene, sequence insertions
"""
from __future__ import annotations


def parse_cigar(cigar):
    """Parse the CIGAR string into a list of tuples."""
    import re

    return [
        (int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


def process_sequence(seq, cigar):
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


#### PSUEDOCODE ####

## read in each line in the sam file

### if the read_id is not in the dictionary, add it

###### make a dictionary for AA_sequences for the gene

###### make a dictionary for insertions for the gene

### get the gene name from the reference name

### get the cleartext sequence from the aligned query sequence and the cigar string

### get the insertions from the aligned query sequence and the cigar string

### add the gene name: sequence, and insertions to the dictionaries

## write the dictionaries to a json file


INPUT_FILE = "diamond_blastx.sam"

# make a dictionary for the AA_sequences
aa_sequences = {}

# make a dictionary for the insertions
aa_insertions = {}

with open(INPUT_FILE, "r") as f:
    for line in f:
        # Skip header lines
        if line.startswith("@"):
            continue
        # Split the line into fields
        fields = line.strip().split("\t")
        read_id = fields[0]
        gene_name = fields[2]
        cigar = fields[5]
        seq = fields[9]

        cleartext, insertions, deletions = process_sequence(seq, cigar)

        if (
            gene_name not in aa_sequences.keys()
            and gene_name not in aa_insertions.keys()
        ):
            aa_sequences[gene_name] = ""
            aa_insertions[gene_name] = ""
