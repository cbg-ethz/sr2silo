"""Implements the read pairing and normalization of the reads."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path

import pysam

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


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

    if not nuc_align_sam_fp.exists():
        raise FileNotFoundError(f"File not found: {nuc_align_sam_fp}")

    if not ref_genome_fasta_fp.exists():
        raise FileNotFoundError(f"File not found: {ref_genome_fasta_fp}")

    if output_merged_sam_fp.exists():
        raise FileExistsError(f"File already exists: {output_merged_sam_fp}")

    if not is_bam_sorted_qname(nuc_align_sam_fp):
        raise ValueError(f"Input file {nuc_align_sam_fp} is not sorted by QNAME.")

    if not had_SQ_header(nuc_align_sam_fp):
        logging.error(f"Input file {nuc_align_sam_fp} does not have @SQ headers.")
        raise ValueError(f"Input file {nuc_align_sam_fp} does not have @SQ headers.")

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
    """Checks if a BAM/SAM file is sorted by QNAME using pysam.

    Args:
        bam_file (str): Path to the BAM/SAM file.

    Returns:
        bool: True if the file is sorted by query name, False otherwise.
        None if there's an issue opening the file.
    """
    try:
        with pysam.AlignmentFile(bam_file, "rb") as af:
            # Header HD should define SO as "queryname" for a QNAME sorted file.
            so = af.header.get("HD", {}).get("SO")  # type: ignore
            return so == "queryname"
    except ValueError as e:
        print(f"Error opening file {bam_file}: {e}")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None


def had_SQ_header(sam_file):
    """Checks if a SAM/BAM file has @SQ headers using pysam.

    Args:
        sam_file (str): Path to the SAM/BAM file.

    Returns:
        bool: True if the file contains @SQ headers, False otherwise.
    """
    try:
        with pysam.AlignmentFile(sam_file, "r") as af:
            # Check if there is an 'SQ' entry in the header
            return "SQ" in af.header and len(af.header["SQ"]) > 0  # type: ignore
    except Exception as e:
        print(f"Error reading SAM file {sam_file}: {e}")
        return False


def sort_sam_by_qname(input_sam_path: Path, output_sam_path: Path):
    """
    Sorts a sam file using pysam.sort command by query name.

    Args:
        input_bam_path (Path): Path to the input BAM file.
        output_bam_path (Path): Path to the output sorted BAM file.
    """
    try:
        # Convert Path objects to strings for pysam compatibility
        input_sam_str = str(input_sam_path)
        output_sam_str = str(output_sam_path)

        # Using pysam.sort command to sort the BAM file and write to disk incrementally.
        pysam.sort("-n", "-o", output_sam_str, input_sam_str)
        logging.info(f"SAM file has been sorted and saved to {output_sam_str}")
    except Exception as e:
        print(f"An error occurred: {e}")
        raise Exception(f"An error occurred: {e}")
