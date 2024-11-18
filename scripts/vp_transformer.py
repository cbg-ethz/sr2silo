"""Converts V-PIPE's outputs of nucleotide read sequencs output to
   ready to import nuclotides and amino acid sequences to the SILO database."""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path

import click

from sr2silo.convert import bam_to_sam
from sr2silo.process import pair_normalize_reads
from sr2silo.translation import translate

logging.basicConfig(level=logging.INFO)


def load_config(config_file: Path) -> dict:
    """Load configuration from a JSON file."""
    return NotImplementedError
    try:
        with config_file.open() as f:
            config = json.load(f)
        logging.debug(f"Loaded config: {config}")
        return config
    except FileNotFoundError:
        logging.error(f"Config file {config_file} not found.")
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from config file {config_file}: {e}")
        raise


def sample_id_decoder(sample_id: str) -> dict:
    """Decode the sample ID into individual components.

    Args:
        sample_id (str): The sample ID to decode.

    Returns:
        dict: A dictionary containing the decoded components.
              containing the following keys:
                - sequencing_well_position (str : sequencing well position)
                - location_code (int : code of the location)
                - sampling_date (str : date of the sampling)
    """
    components = sample_id.split("_")
    # Assign components to meaningful variable names
    well_position = components[0]  # A1
    location_code = components[1]  # 10
    sampling_date = f"{components[2]}-{components[3]}-{components[4]}"  # 2024-09-30
    return {
        "sequencing_well_position": well_position,
        "location_code": location_code,
        "sampling_date": sampling_date,
    }


def batch_id_decoder(batch_id: str) -> dict:
    """Decode the batch ID into individual components.

    Args:
        batch_id (str): The batch ID to decode.

    Returns:
        dict: A dictionary containing the decoded components.
              containing the following keys:
                - sequencing_date (str : date of the sequencing)
                - flow_cell_serial_number (str : serial number of the flow cell)
    """
    components = batch_id.split("_")
    # Assign components to meaningful variable names
    sequencing_date = (
        f"{components[0][:4]}-{components[0][4:6]}-{components[0][6:]}"  # 2024-10-18
    )
    flow_cell_serial_number = components[1]  # AAG55WNM5
    return {
        "sequencing_date": sequencing_date,
        "flow_cell_serial_number": flow_cell_serial_number,
    }


def get_metadata(directory: Path, timeline: Path) -> dict:
    """
    Get metadata for a given sample and batch directory.
    Cross-references the directory with the timeline file to get the metadata.

    Args:
        directory (Path): The directory to extract metadata from.
        timeline (Path): The timeline file to cross-reference the metadata.

    Returns:
        dict: A dictionary containing the metadata.
    """

    # Extract sample and batch IDs from the directory name
    # samples/{sample_id}/{batch_id}/alignments/REF_aln_trim.bam
    metadata = {}
    metadata["sample_id"] = directory.parent.name
    metadata["batch_id"] = directory.name

    # Decompose the ids into individual components
    logging.info(f"Decoding sample_id: {metadata['sample_id']}")
    sample_id = metadata["sample_id"]
    metadata.update(sample_id_decoder(sample_id))
    logging.info(f"Decoding batch_id: {metadata['batch_id']}")
    batch_id = metadata["batch_id"]
    metadata.update(batch_id_decoder(batch_id))

    # Read the timeline file to get additional metadata
    # find row with matching sample_id and batch_id
    # timline has headers:
    #  sample_id	batch_id	read_length	primer_protocol	location_code	sampling_date	location_name
    # get read length, primer protocol, location name
    # double checl if location code and location code are the same
    with timeline.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row[0] == metadata["sample_id"] and row[1] == metadata["batch_id"]:
                logging.info(
                    f"Enriching metadata with timeline data e.g. read_length, primer_protocol, location_name"
                )
                metadata["read_length"] = row[2]
                metadata["primer_protocol"] = row[3]
                metadata["location_name"] = row[6]
                # Convert sampling_date to ISO format for comparison
                timeline_sampling_date = f"{row[5][:4]}-{row[5][4:6]}-{row[5][6:]}"
                if (
                    metadata["location_code"] != row[4]
                    or metadata["sampling_date"] != timeline_sampling_date
                ):
                    logging.warning(
                        f"Mismatch in location code or sampling date for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
                    )
                break
        else:
            raise ValueError(
                f"No matching entry found in timeline for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
            )
    return metadata


def process_directory(
    input_dir: Path,
    result_dir: Path,
    nextclade_reference: str,
    timeline_file: Path,
    file_name: str = "REF_aln_trim.bam",
) -> None:
    """Process all files in a given directory.

    Args:
        input_dir (Path): The directory to process.
        result_dir (Path): The directory to save the results.
        nextclade_reference (str): The reference to use for nextclade.
        timeline_file (Path): The timeline file to cross-reference the metadata.
        file_name (str): The name of the file to process

    Returns:
        None (writes results to the result_dir)
    """

    # check that one was given a directory and not a file and it exists
    if not input_dir.is_dir():
        logging.error(f"Input directory not found, is it a directory?: {input_dir}")
        return

    logging.info(f"Processing directory: {input_dir}")
    logging.info(f"Assuming the input file is: {file_name}")

    # Get Sample and Batch metadata and write to a file
    metadata = get_metadata(input_dir, timeline_file)
    # add nextclade reference to metadata
    metadata["nextclade_reference"] = nextclade_reference
    metadata_file = result_dir / "metadata.json"
    result_dir.mkdir(parents=True, exist_ok=True)
    with metadata_file.open("w") as f:
        json.dump(metadata, f, indent=4)
    logging.info(f"Metadata saved to: {metadata_file}")

    # Convert BAM to SAM
    logging.info(f"Converting BAM to SAM")
    bam_file = input_dir / file_name
    sam_data = bam_to_sam(bam_file)

    # Process SAM to FASTA
    logging.info(f"Processing SAM to FASTA (pair, merge, and normalize reads)")
    fasta_file = result_dir / "reads.fasta"
    insertions_file = result_dir / "insertions.txt"
    pair_normalize_reads(sam_data, fasta_file, insertions_file)

    # Translate nucleotides to amino acids
    logging.info(f"Aliging and translating sequences")
    translate([fasta_file], result_dir, nextclade_reference)

    logging.info(f"Results saved to: {result_dir}")


# TODO: Implement the read_timeline function
def read_timeline(timeline_file: Path) -> list[Path]:
    """Read the timeline.tsv file and return a list of directory paths to process."""
    return NotImplementedError
    directories = []
    with timeline_file.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Assuming the directory path is constructed from metadata in the row
            directory = Path(row["base_dir"]) / row["sub_dir"]
            directories.append(directory)
    return directories


# TODO: Implement the main function
@click.command()
@click.option(
    "--config", default="vp_transformer_config.json", help="Path to the config file."
)
def main(config_file: Path) -> None:
    """Main function to process all subdirectories."""
    return NotImplementedError
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
    # Load the configuration
    config = load_config(Path("scripts/vp_config.json"))

    timeline_file = Path(config["timeline_file"])
    result_dir = Path(config["result_dir"])
    nextclade_reference = config["nextclade_reference"]

    # main(config_file)

    # process a directory: batch / sample
    process_directory(
        Path("../../../data/sr2silo/samples/A1_05_2024_10_08/20241024_2411515907"),
        Path("results"),
        "nextstrain/sars-cov-2/wuhan-hu-1/orfs",
    )
