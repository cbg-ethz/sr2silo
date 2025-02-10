# Script to convert BAM to FASTA with quality scores

from __future__ import annotations

from pathlib import Path

import pysam


import logging


logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def sort_bam_file(input_bam_path: Path, output_bam_path: Path):
    """
    Sorts a BAM file using pysam.sort to avoid loading all alignments into memory.

    Args:
        input_bam_path (Path): Path to the input BAM file.
        output_bam_path (Path): Path to the output sorted BAM file.
    """
    try:
        # Convert Path objects to strings for pysam compatibility
        input_bam_str = str(input_bam_path)
        output_bam_str = str(output_bam_path)

        # Using pysam.sort command to sort the BAM file and write to disk incrementally.
        pysam.sort("-o", output_bam_str, input_bam_str)
        print(f"BAM file has been sorted and saved to {output_bam_str}")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise Exception(f"An error occurred: {e}")


def create_index(bam_file: Path):
    """
    Create an index for a BAM file using pysam.

    Args:
        bam_file (str): Path to the input BAM file.
    """
    try:
        # Convert Path object to string for pysam compatibility
        bam_file_str = str(bam_file)

        # Open BamFile and save with 'bai' extension
        pysam.index(bam_file_str)
        print(f"Index created for {bam_file}")
    except Exception as e:
        print(f"An error occurred: {e}")


def bam_to_fasta(bam_file, fasta_file):
    """
    Convert a BAM file to a FASTA file. Bluntly resolved the sam to fasta.

    Args:
        bam_file (str): Path to the input BAM file.
        fasta_file (str): Path to the output FASTQ file.
    """
    # check for proper format
    if not bam_file.endswith(".bam"):
        raise ValueError("Input file is not a BAM file")
    if not fasta_file.endswith(".fasta"):
        raise ValueError("Output file is not a FASTA file")

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(fasta_file, "w") as fq:
            for read in bam.fetch():
                if not read.is_unmapped:
                    name = read.query_name
                    seq = read.query_sequence
                    fq.write(f">{name}\n{seq}\n")


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
                fastq.write(f"alignment_position:{ref_align_start}\n")

                # Write the insertions to the insertions file
                for insertion_pos, insertion_seq, insertion_qual in insertion_positions:
                    insertions.write(
                        f"{read.query_name}\t{insertion_pos}\t{''.join(insertion_seq)}\t{''.join(insertion_qual)}\n"
                    )


def check_fastq_format(fastq_path):
    """Utility to check that the FASTQ file has 4 lines per record."""
    with open(fastq_path, "r") as f:
        line_count = 0
        # if line starts with @, it is a header, so parse the read id
        for line in f:
            if line.startswith("@"):
                # get the text after tjhe @ symbol
                read_id = line.strip()[1:]
                # get next line which is the sequence
                seq = next(f).strip()
                # skip the next line (which is a +)
                next(f)
                # get the quality scores
                qual = next(f).strip()
                # check if the sequence and quality scores are the same length
                if len(seq) != len(qual):
                    raise ValueError(f"Sequence and quality lines have different lengths for {read_id}")
                if not seq or not qual:
                    raise ValueError(f"Empty sequence or quality line for {read_id}")
                line_count += 4
        if line_count % 4 != 0:
            raise ValueError(
                f"FASTQ format error: {fastq_path} has {line_count} lines (not a multiple of 4)"
            )
        logging.info(f"FASTQ file {fastq_path} has correct format")