"""Process a single BAM file."""

from __future__ import annotations

import logging
from pathlib import Path

import click

from sr2silo.convert import bam_to_sam
from sr2silo.process import pair_normalize_reads
from sr2silo.translation import translate

logging.basicConfig(level=logging.INFO)


@click.command()
@click.option(
    "--file_path",
    type=click.Path(exists=True),
    required=True,
    help="Path to the input BAM file.",
)
@click.option(
    "--output_dir",
    type=click.Path(),
    required=True,
    help="Directory to save the output files.",
)
def process_file(file_path, output_dir):
    """Process a single BAM file.

    Args:
        file_path: Path to the input BAM file.
        output_dir: Directory to save the output files.

    Returns:
        Full output of nextclade translation.
        and currently the insertions and normalised reads.
    """
    directory = Path(file_path).parent
    result_dir = Path(output_dir)
    result_dir.mkdir(parents=True, exist_ok=True)
    logging.info(f"Processing file: {file_path}")

    # Convert BAM to SAM
    logging.info(f"Converting BAM to SAM: {file_path}")
    bam_file = Path(file_path)
    sam_data = bam_to_sam(bam_file)

    # Process SAM to FASTA
    logging.info(f"Processing SAM to FASTA: {file_path}")
    fasta_file = result_dir / "normalized_reads.fasta"
    insertions_file = result_dir / "normalized_insertions.txt"
    pair_normalize_reads(sam_data, fasta_file, insertions_file)

    # Translate nucleotides to amino acids
    logging.info(f"Translating nucleotides to amino acids: {file_path}")
    nextclade_reference = "nextstrain/sars-cov-2/wuhan-hu-1/orfs"
    translate([fasta_file], result_dir, nextclade_reference)


if __name__ == "__main__":
    process_file()
