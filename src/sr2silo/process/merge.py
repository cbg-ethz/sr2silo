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
        paired_end_read_merger,
        str(nuc_align_sam_fp),
        str(ref_genome_fasta_fp),
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
