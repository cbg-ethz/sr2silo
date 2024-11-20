"""Converts V-PIPE's outputs of nucleotide read sequencs output to
   ready to import nuclotides and amino acid sequences to the SILO database."""

from __future__ import annotations

import csv
import datetime
import json
import logging
from pathlib import Path

import click

from sr2silo.convert import bam_to_sam
from sr2silo.process import pair_normalize_reads
from sr2silo.translation import translate

logging.basicConfig(
    level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s"
)


def load_config(config_file: Path) -> dict:
    """Load a JSON configuration file."""
    try:
        with config_file.open() as f:
            return json.load(f)
    except FileNotFoundError:
        logging.error(f"Config file not found: {config_file}")
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from config file: {config_file} - {e}")
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


def convert_to_iso_date(date: str) -> str:
    """Convert a date string to ISO 8601 format (date only)."""
    # Parse the date string
    date_obj = datetime.datetime.strptime(date, "%Y-%m-%d")
    # Format the date as ISO 8601 (date only)
    return date_obj.date().isoformat()


def get_metadata(sample_id: str, batch_id: str, timeline: Path) -> dict:
    """
    Get metadata for a given sample and batch directory.
    Cross-references the directory with the timeline file to get the metadata.

    Args:
        sample_id (str): The sample ID to use for metadata.
        batch_id (str): The batch ID to use for metadata.
        timeline (Path): The timeline file to cross-reference the metadata.

    Returns:
        dict: A dictionary containing the metadata.

    """

    metadata = {}
    metadata["sample_id"] = sample_id
    metadata["batch_id"] = batch_id

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
                timeline_sampling_date = convert_to_iso_date(row[5])
                if int(metadata["location_code"]) != int(row[4]):
                    # output both location codes for comparison and their types for debugging
                    logging.warning(
                        f"Mismatch in location code for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
                    )
                    logging.debug(
                        f"Location code mismatch: {metadata['location_code']} (sample_id) vs {row[4]} (timeline)"
                    )
                    logging.debug(
                        f"Location code types: {type(metadata['location_code'])} (sample_id) vs {type(row[4])} (timeline)"
                    )
                if metadata["sampling_date"] != timeline_sampling_date:
                    # output both sampling dates for comparison and their types for debugging
                    logging.warning(
                        f"Mismatch in sampling date for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
                    )
                    logging.debug(
                        f"Sampling date mismatch: {metadata['sampling_date']} (sample_id) vs {timeline_sampling_date} (timeline)"
                    )
                    logging.debug(
                        f"Sampling date types: {type(metadata['sampling_date'])} (sample_id) vs {type(timeline_sampling_date)} (timeline)"
                    )
                break
        else:
            raise ValueError(
                f"No matching entry found in timeline for sample_id {metadata['sample_id']} and batch_id {metadata['batch_id']}"
            )
    return metadata


def process_directory(
    input_dir: Path,
    sample_id: str,
    batch_id: str,
    result_dir: Path,
    nextclade_reference: str,
    timeline_file: Path,
    file_name: str = "REF_aln_trim.bam",
) -> None:
    """Process all files in a given directory.

    Args:
        input_dir (Path): The directory to process. i.e. the directory containing the BAM file.
                          to reach samples/A1_05_2024_10_08/20241024_2411515907/alignments/
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
        raise FileNotFoundError(f"Directory not found: {input_dir}")

    logging.info(f"Processing directory: {input_dir}")
    logging.info(f"Assuming the input file is: {file_name}")
    # check that the file exists and also it's .bai file
    sample_fp = input_dir / file_name
    if not sample_fp.exists():
        logging.error(f"Input file not found: {sample_fp}")
        raise FileNotFoundError(f"Input file not found: {sample_fp}")

    # Get Sample and Batch metadata and write to a file
    metadata = get_metadata(sample_id, batch_id, timeline_file)
    # add nextclade reference to metadata
    metadata["nextclade_reference"] = nextclade_reference
    metadata_file = result_dir / "metadata.json"
    result_dir.mkdir(parents=True, exist_ok=True)
    with metadata_file.open("w") as f:
        json.dump(metadata, f, indent=4)
    logging.info(f"Metadata saved to: {metadata_file}")

    # Convert BAM to SAM
    logging.info(f"Converting BAM to SAM")
    bam_file = sample_fp
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


@click.command()
@click.option("--sample_dir", envvar="SAMPLE_DIR", help="Path to the sample directory.")
@click.option("--sample_id", envvar="SAMPLE_ID", help="sample_id to use for metadata.")
@click.option("--batch_id", envvar="BATCH_ID", help="batch_id to use for metadata.")
@click.option(
    "--result_dir", envvar="RESULTS_DIR", help="Path to the results directory."
)
@click.option(
    "--timeline_file", envvar="TIMELINE_FILE", help="Path to the timeline file."
)
@click.option(
    "--nextclade_reference",
    envvar="NEXTCLADE_REFERENCE",
    default="sars-cov-2",
    help="Nextclade reference.",
)
def main(
    sample_dir, sample_id, batch_id, result_dir, timeline_file, nextclade_reference
):
    """Process a sample directory."""
    logging.info(f"Processing sample directory: {sample_dir}")
    logging.info(f"Saving results to: {result_dir}")
    logging.info(f"Using timeline file: {timeline_file}")
    logging.info(f"Using Nextclade reference: {nextclade_reference}")
    logging.info(f"Using sample_id: {sample_id}")
    logging.info(f"Using batch_id: {batch_id}")

    process_directory(
        input_dir=Path("sample"),
        sample_id=sample_id,
        batch_id=batch_id,
        result_dir=Path("results"),
        timeline_file=Path("timeline.tsv"),
        nextclade_reference=nextclade_reference,
    )


if __name__ == "__main__":
    main()
