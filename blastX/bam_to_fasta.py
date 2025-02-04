# Script to convert BAM to FASTA with quality scores

from __future__ import annotations

import pysam


def sort_bam_file(input_bam_path, output_bam_path):
    """
    Sorts a BAM file by coordinate.

    :param input_bam_path: Path to the input BAM file.
    :param output_bam_path: Path where the sorted BAM file will be saved.
    """
    try:
        # Open the input BAM file
        with pysam.AlignmentFile(input_bam_path, "rb") as bam_file:
            # Sort the alignments by coordinate and write them to a new file
            sorted_alignments = bam_file.fetch(until_eof=True)
            with pysam.AlignmentFile(
                output_bam_path, "wb", template=bam_file
            ) as output_file:
                for alignment in sorted(
                    sorted_alignments,
                    key=lambda x: (x.reference_name, x.reference_start),
                ):
                    output_file.write(alignment)

        print(f"BAM file has been sorted and saved to {output_bam_path}")

    except Exception as e:
        print(f"An error occurred: {e}")


def create_index(bam_file):
    """
    Create an index for a BAM file using pysam.

    Args:
        bam_file (str): Path to the input BAM file.
    """
    try:
        # Open BamFile and save with 'bai' extension
        pysam.index(bam_file)
        print(f"Index created for {bam_file}")
    except Exception as e:
        print(f"An error occurred: {e}")


def bam_to_fastq(bam_file, fastq_file):
    """
    Convert a BAM file to a FASTQ file.

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
                    qual = "".join(chr(ord("!") + q) for q in read.query_qualities)
                    fq.write(f"@{name}\n{seq}\n+\n{qual}\n")


def bam_to_fastq_handle_indels(
    bam_file, fastq_file, insertions_file, deletion_char="-"
):
    """
    Convert a BAM file to a FASTQ file, removing insertions and adding a special character for deletions.
    Save the insertions to a separate file.

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
                        insertion_qual = query_qualities[
                            query_pos : query_pos + cigar[1]
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
                fastq.write(f"{''.join(chr(q + 33) for q in new_qualities)}\n")

                # Write the insertions to the insertions file
                for insertion_pos, insertion_seq, insertion_qual in insertion_positions:
                    insertions.write(
                        f"{read.query_name}\t{insertion_pos}\t{''.join(insertion_seq)}\t{''.join(chr(q + 33) for q in insertion_qual)}\n"
                    )


# Example usage:

sort_bam_file("input/combined.bam", "input/sorted.bam")
create_index("input/sorted.bam")
bam_to_fastq("input/sorted.bam", "output.fastq")
