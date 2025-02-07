"""Script to convert SAM to json gene, sequence insertions
"""
from __future__ import annotations

import json
import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple

import bam_to_fasta
import pysam

from sr2silo.process import pad_alignment

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """Parse the CIGAR string into a list of tuples."""
    return [
        (int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)
    ]


class NucInsertion:
    def __init__(self, position: int, sequence: str):
        self.position = position
        self.sequence = sequence

    def __str__(self) -> str:
        return f"{self.position} : {self.sequence}"


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


INPUT_NUC_ALIGMENT_FILE = "input/combined.bam"
# INPUT_NUC_ALIGMENT_FILE = "input/Ref_aln.bam"
FASTQ_NUC_ALIGMENT_FILE = "output.fastq"
FASTA_NUC_INSERTIONS_FILE = "output_ins.fasta"
AA_ALIGMENT_FILE = "diamond_blastx.sam"
AA_REFERENCE_FILE = "../resources/sars-cov-2/aa_reference_genomes.fasta"
NUC_REFERENCE_FILE = "../resources/sars-cov-2/nuc_reference_genomes.fasta"

bam_file = Path("input/sorted.bam")
bam_to_fasta.sort_bam_file(INPUT_NUC_ALIGMENT_FILE, bam_file)

bai_file = bam_file.with_suffix(".bai")
if not bai_file.exists() or bam_file.stat().st_mtime > bai_file.stat().st_mtime:
    print("Creating index for BAM file")
    bam_to_fasta.create_index(str(bam_file))

print("Converting BAM to FASTQ with INDELS to show NUC alignment")

# make a mofidifed path for indels
FASTQ_NUC_ALIGMENT_FILE_WI = "output_with_indels.fastq"
FASTA_NUC_INSERTIONS_FILE_WI = "output_ins_with_indels.fasta"
bam_to_fasta.bam_to_fastq_handle_indels(
    "input/sorted.bam", FASTQ_NUC_ALIGMENT_FILE_WI, FASTA_NUC_INSERTIONS_FILE_WI
)

print("Converting BAM to FASTQ for AA alignment")
bam_to_fasta.bam_to_fastq("input/sorted.bam", FASTQ_NUC_ALIGMENT_FILE)

try:
    # translate and align to AA
    # ==== Make Sequence DB ====
    print("== Making Sequence DB ==")
    result = os.system(f"diamond makedb --in {AA_REFERENCE_FILE} -d ref/hxb_pol_db")
    if result != 0:
        raise RuntimeError(
            "Error occurred while making sequence DB with diamond makedb"
        )

    # ==== Alignment ====
    print("== Aligning to AA ==")
    result = os.system(
        f"""
    diamond blastx -d ref/hxb_pol_db \
        -q {FASTQ_NUC_ALIGMENT_FILE} \
        -o {AA_ALIGMENT_FILE} \
        --evalue 1 \
        --gapopen 6 \
        --gapextend 2 \
        --outfmt 101 \
        --matrix BLOSUM62 \
        --unal 0 \
        --max-hsps 1 \
        --block-size 0.5
    """
    )
    if result != 0:
        raise RuntimeError("Error occurred while aligning to AA with diamond blastx")
except Exception as e:
    print(f"An error occurred: {e}")
    raise


with open(NUC_REFERENCE_FILE, "r") as f:
    nuc_reference = f.read()
nuc_reference_length = len(nuc_reference)

gene_dict = get_genes_and_lengths_from_ref(AA_REFERENCE_FILE)

reads: List[AlignedRead] = []

# load in the nuc insertions file - NucInsertion(position, sequence)
nuc_insertions: Dict[str, List[NucInsertion]] = dict()

with open(FASTA_NUC_INSERTIONS_FILE, "r") as f:
    for line in f:
        if line.startswith(">"):
            continue
        fields = line.strip().split("\t")
        read_id = fields[0]
        position = int(fields[1])
        sequence = fields[2]
        nuc_ins = NucInsertion(position, sequence)
        if read_id not in nuc_insertions:
            nuc_insertions[read_id] = []
        nuc_insertions[read_id].append(nuc_ins)


with pysam.AlignmentFile("input/sorted.bam", "rb") as bam:
    for entry in bam:
        for read in bam.fetch():
            read_id = read.query_name
            seq = read.query_sequence
            qual = "".join(chr(ord("!") + q) for q in read.query_qualities)
            pos = read.qstart

            aligned_nucleotide_sequences = pad_alignment(seq, pos, nuc_reference_length)

            if read_id in nuc_insertions:
                nucleotide_insertions = nuc_insertions[read_id]
            else:
                nucleotide_insertions = []

            reads.append(
                AlignedRead(
                    read_id=read_id,
                    unaligned_nucleotide_sequences=seq,
                    aligned_nucleotide_sequences=aligned_nucleotide_sequences,
                    nucleotide_insertions=nucleotide_insertions,
                    amino_acid_insertions="null",
                    aligned_amino_acid_sequences="null",
                )
            )

with open(AA_ALIGMENT_FILE, "r") as f:
    count = 0
    for line in f:
        count += 1
        if count % 1000 == 0:
            print(f"AA Alignment of read {count} ")
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
        # convert aa_insertions to dict of all gene names and add insertions to the correct gene

        aa_insertions_fmt = {}
        for gene in gene_dict.keys():
            aa_insertions_fmt[gene] = []
        aa_insertions_fmt[gene_name] = aa_insertions

        aa_insertions = AAInsertionSet(gene_dict, aa_insertions, gene_name)

        # Make a dict to hold the aligned amino acid sequences
        aligned_amino_acid_sequences = {}
        # Write a null for all gene names
        for gene in get_genes_and_lengths_from_ref(AA_REFERENCE_FILE).keys():
            aligned_amino_acid_sequences[gene] = None

        # pad the alignment
        padded_aa_alignment = pad_alignment(seq, pos, gene_dict[gene_name].gene_length)

        # find the correct read by read_id in reads
        for read in reads:
            if read.read_id == read_id:
                read.aligned_amino_acid_sequences = padded_aa_alignment
                read.amino_acid_insertions = aa_insertions
                break

print(reads[0].to_json())
print(reads[-1].to_json())
print(f"Reads: {len(reads)}")
print("Done!")
