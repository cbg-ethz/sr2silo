"""Implements the read pairing and normalization of the reads."""

from __future__ import annotations

import logging
import subprocess
from pathlib import Path


from sr2silo.process.convert import is_bam_sorted_qname, had_SQ_header

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

    if not had_SQ_header(nuc_align_sam_fp):  # pragma: no cover
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
    except subprocess.CalledProcessError as e:  # pragma: no cover
        logger.error(f"Error: {e}")
        raise e
    except Exception as e:  # pragma: no cover
        logger.error(f"Error: {e}")
        raise e

