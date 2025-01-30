# Script to convert SAM to json gene, sequence insertions
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


# TODO: validate the long names
protein_names = {
    "YP_009742608.1": "leader protein",
    "YP_009742609.1": "nsp2",
    "YP_009742610.1": "nsp3",
    "YP_009742611.1": "nsp4",
    "YP_009742612.1": "3C-like proteinase",
    "YP_009742613.1": "nsp6",
    "YP_009742614.1": "nsp7",
    "YP_009742615.1": "nsp8",
    "YP_009742616.1": "nsp9",
    "YP_009742617.1": "nsp10",
    "YP_009725318.1": "ORF7b",
    "YP_009725295.1": "ORF1a polyprotein",
    "YP_009725297.1": "leader protein",
    "YP_009725298.1": "nsp2",
    "YP_009725299.1": "nsp3",
    "YP_009725300.1": "nsp4",
    "YP_009725301.1": "3C-like proteinase",
    "YP_009725302.1": "nsp6",
    "YP_009725303.1": "nsp7",
    "YP_009725304.1": "nsp8",
    "YP_009725305.1": "nsp9",
    "YP_009725306.1": "nsp10",
    "YP_009725307.1": "RNA-dependent RNA polymerase",
    "YP_009725308.1": "helicase",
    "YP_009725309.1": "3'-to-5' exonuclease",
    "YP_009725310.1": "endoRNAse",
    "YP_009725311.1": "2'-O-ribose methyltransferase",
    "YP_009725312.1": "nsp11",
    "YP_009725255.1": "ORF10 protein",
    "YP_009724389.1": "ORF1ab polyprotein",
    "YP_009724390.1": "surface glycoprotein",
    "YP_009724391.1": "ORF3a protein",
    "YP_009724392.1": "envelope protein",
    "YP_009724393.1": "membrane glycoprotein",
    "YP_009724394.1": "ORF6 protein",
    "YP_009724395.1": "ORF7a protein",
    "YP_009724396.1": "ORF8 protein",
    "YP_009724397.2": "nucleocapsid phosphoprotein",
}

# TODO: validate the short names
short_names = {
    "leader protein": "leader",
    "nsp2": "nsp2",
    "nsp3": "nsp3",
    "nsp4": "nsp4",
    "3C-like proteinase": "3CLpro",
    "nsp6": "nsp6",
    "nsp7": "nsp7",
    "nsp8": "nsp8",
    "nsp9": "nsp9",
    "nsp10": "nsp10",
    "ORF7b": "ORF7b",
    "ORF1a polyprotein": "ORF1a",
    "RNA-dependent RNA polymerase": "RdRp",
    "helicase": "Helicase",
    "3'-to-5' exonuclease": "ExoN",
    "endoRNAse": "EndoU",
    "2'-O-ribose methyltransferase": "NMT",
    "nsp11": "nsp11",
    "ORF10 protein": "ORF10",
    "ORF1ab polyprotein": "ORF1ab",
    "surface glycoprotein": "S",
    "ORF3a protein": "ORF3a",
    "envelope protein": "E",
    "membrane glycoprotein": "M",
    "ORF6 protein": "ORF6",
    "ORF7a protein": "ORF7a",
    "ORF8 protein": "ORF8",
    "nucleocapsid phosphoprotein": "N",
}

protein_short_names = {}

for key, value in protein_names.items():
    if value in short_names:
        protein_short_names[key] = short_names[value]
    else:
        protein_short_names[
            key
        ] = value  # Use the full name if no short name is available

print(protein_short_names)

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


for ref_name in aa_sequences.keys():
    print(protein_short_names[ref_name])
