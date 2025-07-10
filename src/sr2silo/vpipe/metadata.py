"""Extract metadata from V-Pipe Timeline Files."""

from __future__ import annotations

import csv
import datetime
import logging
from pathlib import Path


def convert_to_iso_date(date: str) -> str:
    """Convert a date string to ISO 8601 format (date only)."""
    # Parse the date string
    date_obj = datetime.datetime.strptime(date, "%Y-%m-%d")
    # Format the date as ISO 8601 (date only)
    return date_obj.date().isoformat()


def get_metadata_from_timeline(
    sample_id: str, batch_id: str, timeline: Path
) -> dict[str, str]:
    """Get metadata from the timeline file for a given sample and batch.

    Args:
        sample_id (str): The sample ID to use for metadata.
        batch_id (str): The batch ID to use for metadata.  
        timeline (Path): The timeline file to get the metadata from.

    Returns:
        dict: A dictionary containing the metadata.

    """
    if not timeline.is_file():
        logging.error(f"Timeline file not found or is not a file: {timeline}")
        raise FileNotFoundError(f"Timeline file not found or is not a file: {timeline}")

    metadata = {
        "sample_id": sample_id,
        "batch_id": batch_id,
    }

    with timeline.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            # Timeline format: sample, batch, reads, proto, location_code, date, location
            if len(row) >= 7 and row[0] == sample_id and row[1] == batch_id:
                logging.info(
                    f"Found timeline entry for sample_id {sample_id} and batch_id {batch_id}"
                )
                
                # Extract metadata from timeline with graceful handling of missing data
                metadata["read_length"] = row[2] if row[2] else None
                metadata["primer_protocol"] = row[3] if row[3] else None
                metadata["location_code"] = row[4] if row[4] else None
                
                # Handle sampling date
                if row[5]:
                    try:
                        metadata["sampling_date"] = convert_to_iso_date(row[5])
                    except ValueError as e:
                        logging.warning(f"Invalid date format in timeline: {row[5]}, error: {e}")
                        metadata["sampling_date"] = None
                else:
                    metadata["sampling_date"] = None
                    
                metadata["location_name"] = row[6] if row[6] else None

                return metadata
                
        # If we reach here, no matching entry was found
        logging.warning(
            f"No matching entry found in timeline for sample_id "
            f"{sample_id} and batch_id {batch_id}"
        )
        # Return metadata with None values for missing fields
        metadata.update({
            "read_length": None,
            "primer_protocol": None,
            "location_code": None,
            "sampling_date": None,
            "location_name": None,
        })
        return metadata


def get_metadata(
    sample_id: str, batch_id: str, timeline: Path
) -> dict[str, str]:
    """
    Get metadata for a given sample and batch from timeline file only.

    Args:
        sample_id (str): The sample ID to use for metadata.
        batch_id (str): The batch ID to use for metadata.
        timeline (Path): The timeline file to get the metadata from.

    Returns:
        dict: A dictionary containing the metadata.

    """
    return get_metadata_from_timeline(sample_id, batch_id, timeline)
