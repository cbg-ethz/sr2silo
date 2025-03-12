"""Implements the read pairing and normalization of the reads."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path


logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# TODO: check that the input file is sorted by QNAME and needs @SQ header
def paired_end_read_merger(
    nuc_align_sam_fp: Path, ref_genome_fasta_fp: Path, output_merged_sam_fp: Path
) -> None:
    """Merge paired-end reads using `smallgenomutilities`tool `paired_end_read_merger`.

    Args:
        nuc_align_sam_fp: Path to the nucleotide alignment SAM file.
        ref_genome_fasta_fp: Path to the reference genome FASTA file.
        output_merged_sam_fp: Path to the output merged SAM file.

    Returns:
        None
    """

    command = [
        "paired_end_read_merger",
        str(nuc_align_sam_fp),
        "-f",
        str(ref_genome_fasta_fp),
        "-o",
        str(output_merged_sam_fp),
    ]

    logger.info(f"Pairing reads based on the alignment file: {nuc_align_sam_fp}")
    logger.debug(f"Running command: {' '.join(command)}")

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error: {e}")
        raise e
    except Exception as e:
        logger.error(f"Error: {e}")
        raise e


def is_bam_sorted_qname(bam_file):
    """Checks if a BAM file is sorted by genomic coordinates using pysam.

    Args:
        bam_file (str): Path to the BAM file.

    Returns:
        bool: True if the BAM file is sorted, False otherwise.
        Returns None if there's an issue opening the file.
    """
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")  # Open in read-binary mode
        is_sorted = (
            bam.header.get("HD", {}).get("SO") == "coordinate"  # pyright: ignore
        )  # pyright: ignore
        bam.close()  # Important: Close the file!
        return is_sorted
    except ValueError as e:
        print(f"Error opening BAM file {bam_file}: {e}")  # Handle file errors
        return None  # Indicate an issue
    except Exception as e:  # Catch other potential errors (like missing index)
        print(f"An unexpected error occurred: {e}")
        return None


def had_SQ_header():
    """Check if the input file has a @SQ header."""

    return NotImplementedError
