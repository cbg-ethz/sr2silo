"""Converts V-PIPE's outputs of nucleotide read sequencs output to
   ready to import nuclotides and amino acid sequences to the SILO database."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

from sr2silo.convert import bam_to_sam
from sr2silo.process import pair_normalize_reads
from sr2silo.translation import translate

logging.basicConfig(level=logging.DEBUG)


def load_config(config_file: Path) -> dict:
    """Load configuration from a JSON file."""
    with config_file.open() as f:
        config = json.load(f)
    return config


def process_directory(
    directory: Path, result_dir: Path, nextclade_reference: str
) -> None:
    """Process all files in a given directory."""
    logging.info(f"Processing directory: {directory}")

    # Read metadata from JSON file
    metadata_file = directory / "metadata.json"
    with metadata_file.open() as f:
        metadata = json.load(f)

    # Convert BAM to SAM
    bam_file = directory / "reads.bam"
    sam_data = bam_to_sam(bam_file)

    # Process SAM to FASTA
    fasta_file = directory / "reads.fasta"
    insertions_file = directory / "insertions.txt"
    pair_normalize_reads(sam_data, fasta_file, insertions_file)

    # Translate nucleotides to amino acids
    translate([fasta_file], result_dir, nextclade_reference)

    # Save metadata and results
    result_metadata_file = result_dir / f"{directory.name}_metadata.json"
    with result_metadata_file.open("w") as f:
        json.dump(metadata, f, indent=4)


def read_timeline(timeline_file: Path) -> list[Path]:
    """Read the timeline.tsv file and return a list of directory paths to process."""
    directories = []
    with timeline_file.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Assuming the directory path is constructed from metadata in the row
            directory = Path(row["base_dir"]) / row["sub_dir"]
            directories.append(directory)
    return directories


def main(config_file: Path) -> None:
    """Main function to process all subdirectories."""
    config = load_config(config_file)
    base_dir = Path(config["base_dir"])
    result_dir = Path(config["result_dir"])
    timeline_file = Path(config["timeline_file"])

    directories = read_timeline(timeline_file)
    for subdir in directories:
        if subdir.is_dir():
            logging.debug(f"Processing directory: {subdir}")
            # process_directory(subdir, result_dir)


if __name__ == "__main__":
    config_file = Path("vp_transformer_config.json")
    main(config_file)
