"""Extract metadata from V-Pipe Filenaming Conventions."""

from __future__ import annotations

import csv
import datetime
import logging
from pathlib import Path

import yaml


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
        batch_id (str): The batch ID to decode. Can be empty string.

    Returns:
        dict: A dictionary contains the decoded components.
              dict: A dictionary contains the decoded components.
                - sequencing_date (str : date of the sequencing or empty)
                - flow_cell_serial_number (str : serial number of the flow cell or empty)
    """
    if not batch_id or batch_id.strip() == "":
        return {
            "sequencing_date": "",
            "flow_cell_serial_number": "",
        }

    components = batch_id.split("_")
    if len(components) < 2:
        # Handle malformed batch_id
        return {
            "sequencing_date": "",
            "flow_cell_serial_number": "",
        }

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


def enrich_metadata_from_timeline(metadata: dict[str, str], timeline: Path) -> None:
    """Enrich metadata from the timeline file."""
    if not timeline.is_file():
        logging.error(f"Timeline file not found or is not a file: {timeline}")
        raise FileNotFoundError(f"Timeline file not found or is not a file: {timeline}")

    with timeline.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            sample_id_match = row[0] == metadata["sample_id"]

            if sample_id_match:
                logging.info(
                    "Enriching metadata with timeline data e.g. read_length, "
                    "primer_protocol, location_name"
                )
                metadata["read_length"] = row[2]
                metadata["primer_protocol"] = row[3]
                metadata["location_name"] = row[6]

                # Convert sampling_date to ISO format for comparison
                timeline_sampling_date = convert_to_iso_date(row[5])

                location_code_mismatch = int(metadata["location_code"]) != int(row[4])
                sampling_date_mismatch = (
                    metadata["sampling_date"] != timeline_sampling_date
                )

                if location_code_mismatch:
                    logging.warning(
                        f"Mismatch in location code for sample_id "
                        f"{metadata['sample_id']}"
                    )
                    logging.debug(
                        f"Location code mismatch: {metadata['location_code']} "
                        f"(sample_id) vs {row[4]} (timeline)"
                    )
                    logging.debug(
                        f"Location code types: {type(metadata['location_code'])} "
                        f"(sample_id) vs {type(row[4])} (timeline)"
                    )

                if sampling_date_mismatch:
                    logging.warning(
                        f"Mismatch in sampling date for sample_id "
                        f"{metadata['sample_id']}"
                    )
                    logging.debug(
                        f"Sampling date mismatch: {metadata['sampling_date']} "
                        f"(sample_id) vs {timeline_sampling_date} (timeline)"
                    )
                    logging.debug(
                        f"Sampling date types: {type(metadata['sampling_date'])} "
                        f"(sample_id) vs {type(timeline_sampling_date)} (timeline)"
                    )
                break
        else:
            raise ValueError(
                f"No matching entry found in timeline for sample_id "
                f"{metadata['sample_id']}"
            )


def get_primer_protocol_name(primer_protocol: str, primers: Path) -> str:
    """Get the name of the primer protocol from the primers file.

    Args:
        primer_protocol (str): The primer protocol short name.
        primers (Path): The primers file to with the long, canonical name

    Returns:
        str: The long name of the primer protocol.
    """
    if not primers.is_file():
        logging.error(f"Primers file not found or is not a file: {primers}")
        raise FileNotFoundError(f"Primers file not found or is not a file: {primers}")

    # Load YAML file
    with primers.open() as f:
        primers_conf = yaml.safe_load(f)

    for primer in primers_conf.keys():
        if primer == primer_protocol:
            return primers_conf[primer]["name"]

    raise ValueError(
        f"No matching entry found in primers for primer_protocol {primer_protocol}"
    )


def get_metadata(
    sample_id: str, batch_id: str | None, timeline: Path, primers: Path
) -> dict[str, str]:
    """
    Get metadata for a given sample and batch directory.
    Cross-references the directory with the timeline file to get the metadata.

    Args:
        sample_id (str): The sample ID to use for metadata.
        batch_id (str | None): The batch ID to use for metadata. Can be None or empty.
        timeline (Path): The timeline file to cross-reference the metadata.
        primers (Path): The primers file to cross-reference the metadata.

    Returns:
        dict: A dictionary containing the metadata.

    """

    metadata = {}
    metadata["sample_id"] = sample_id

    # Handle None or empty batch_id
    if batch_id is None:
        batch_id = ""
    metadata["batch_id"] = batch_id

    # Decompose the ids into individual components
    logging.info(f"Decoding sample_id: {metadata['sample_id']}")
    sample_id = metadata["sample_id"]
    metadata.update(sample_id_decoder(sample_id))
    logging.info(f"Decoding batch_id: '{metadata['batch_id']}'")
    batch_id = metadata["batch_id"]
    metadata.update(batch_id_decoder(batch_id))

    enrich_metadata_from_timeline(metadata, timeline)

    metadata["primer_protocol_name"] = get_primer_protocol_name(
        metadata["primer_protocol"], primers
    )

    return metadata
