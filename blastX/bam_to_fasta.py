# Script to convert BAM to FASTA with quality scores

from __future__ import annotations

import os
from pathlib import Path

import pysam


def sort_bam_file(input_bam_path: Path, output_bam_path: Path):
    """
    Sorts a BAM file using pysam.sort to avoid loading all alignments into memory.
    """
    try:
        # Use os.fspath() to convert Path to str for pysam.sort.
        pysam.sort("-o", os.fspath(output_bam_path), os.fspath(input_bam_path))
        print(f"BAM file has been sorted and saved to {output_bam_path}")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise Exception(f"An error occurred: {e}")


def create_index(bam_file):
    """
    Create an index for a BAM file using pysam.

    Args:
        bam_file (str): Path to the input BAM file.
    """
    try:
        # Convert bam_file to str if it is a Path object.
        pysam.index(os.fspath(bam_file))
        print(f"Index created for {bam_file}")
    except Exception as e:
        print(f"An error occurred: {e}")


def bam_to_fastq(bam_file, fastq_file):
    """
    Convert a BAM file to a FASTQ file. Bluntly resolved the sam to fastq.

    Args:
        bam_file (str): Path to the input BAM file.
        fastq_file (str): Path to the output FASTQ file.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(fastq_file, "w") as fq:
            for read in bam.fetch():
                if not read.is_unmapped:
                    name = read.query_name
                    seq = read.query_sequence
                    qual = "".join(chr(q + 33) for q in read.query_qualities)
                    fq.write(f"@{name}\n{seq}\n+\n{qual}\n")


def bam_to_fastq_handle_indels(
    bam_file, fastq_file, insertions_file, deletion_char="-"
):
    """
    Convert a BAM file to a FASTQ file, removing insertions and adding a special character for deletions.
    Save the insertions to a separate file. Include alignment positions in the FASTQ file.

    Used to look at the cleartext nucleotide sequence of the reads.

    :param bam_file: Path to the input BAM file
    :param fastq_file: Path to the output FASTQ file
    :param insertions_file: Path to the output file containing insertions
    :param deletion_char: Special character to use for deletions
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam, open(
        fastq_file, "w"
    ) as fastq, open(insertions_file, "w") as insertions:
        for read in bam.fetch():
            if not read.is_unmapped:
                query_sequence = read.query_sequence
                query_qualities = read.query_qualities
                new_sequence = []
                new_qualities = []
                insertion_positions = []

                query_pos = 0
                ref_align_start = read.reference_start
                ref_pos = read.reference_start

                for cigar in read.cigartuples:
                    if cigar[0] == 0:  # Match or mismatch
                        new_sequence.extend(
                            query_sequence[query_pos : query_pos + cigar[1]]
                        )
                        new_qualities.extend(
                            query_qualities[query_pos : query_pos + cigar[1]]
                        )
                        query_pos += cigar[1]
                        ref_pos += cigar[1]
                    elif cigar[0] == 1:  # Insertion
                        insertion_seq = query_sequence[query_pos : query_pos + cigar[1]]
                        insertion_qual = [
                            chr(q + 33)
                            for q in query_qualities[query_pos : query_pos + cigar[1]]
                        ]
                        insertion_positions.append(
                            (ref_pos, insertion_seq, insertion_qual)
                        )
                        query_pos += cigar[1]
                    elif cigar[0] == 2:  # Deletion
                        new_sequence.extend([deletion_char] * cigar[1])
                        new_qualities.extend(
                            [0] * cigar[1]
                        )  # Assigning a low-quality score for deletions
                        ref_pos += cigar[1]

                # Write the modified read to the FASTQ file
                fastq.write(f"@{read.query_name}\n")
                fastq.write(f"{''.join(new_sequence)}\n")
                fastq.write("+\n")
                fastq.write(
                    f"{''.join(chr(q + 33) for q in new_qualities)}\n"
                )  # Phred33 encoding

                # Write the alignment positions to the FASTQ file
                fastq.write(f"aligment_position:{ref_align_start}\n")

                # Write the insertions to the insertions file
                for insertion_pos, insertion_seq, insertion_qual in insertion_positions:
                    insertions.write(
                        f"{read.query_name}\t{insertion_pos}\t{''.join(insertion_seq)}\t{''.join(insertion_qual)}\n"
                    )
