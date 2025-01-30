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


# Example usage
input_bam = "path/to/your/input.bam"
output_bam = "path/to/your/output.sorted.bam"
sort_bam_file(input_bam, output_bam)


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


# Example usage:

sort_bam_file("input/combined.bam", "input/sorted.bam")
create_index("input/sorted.bam")
bam_to_fastq("input/sorted.bam", "output.fastq")
