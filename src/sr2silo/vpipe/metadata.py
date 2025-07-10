"""Extract metadata from Timeline file only."""

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
) -> dict[str, str] | None:
    """Get metadata from the timeline file.
    
    Args:
        sample_id (str): The sample ID to search for.
        batch_id (str): The batch ID to search for.
        timeline (Path): The timeline file to search in.
        
    Returns:
        dict[str, str] | None: The metadata if found, None otherwise.
    """
    if not timeline.is_file():
        logging.error(f"Timeline file not found or is not a file: {timeline}")
        raise FileNotFoundError(f"Timeline file not found or is not a file: {timeline}")

    with timeline.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            sample_id_match = row[0] == sample_id
            batch_id_match = row[1] == batch_id

            if sample_id_match:
                logging.info(
                    "Found metadata in timeline for sample_id %s and batch_id %s",
                    sample_id, batch_id
                )
                
                # Extract metadata from timeline row
                # Timeline format: sample	batch	reads	proto	location_code	date	location
                metadata = {
                    "sample_id": sample_id,
                    "batch_id": batch_id,
                    "read_length": row[2],
                    "primer_protocol": row[3],
                    "location_code": row[4],
                    "sampling_date": convert_to_iso_date(row[5]),
                    "location_name": row[6],
                }
                
                logging.info("Extracted metadata from timeline: %s", metadata)
                return metadata
        
        # No matching entry found
        logging.warning(
            "No matching entry found in timeline for sample_id %s and batch_id %s",
            sample_id, batch_id
        )
        return None


def get_metadata(
    sample_id: str, batch_id: str, timeline: Path
) -> dict[str, str]:
    """
    Get metadata for a given sample and batch from timeline file only.

    Args:
        sample_id (str): The sample ID to use for metadata.
        batch_id (str | None): The batch ID to use for metadata. Can be None or empty.
        timeline (Path): The timeline file to cross-reference the metadata.

    Returns:
        dict: A dictionary containing the metadata, or empty dict if not found.

    """
    metadata = get_metadata_from_timeline(sample_id, batch_id, timeline)
    
    if metadata is None:
        # Return basic metadata structure if not found in timeline
        logging.warning(
            "Timeline entry not found for sample_id %s and batch_id %s. "
            "Returning basic metadata structure with None values.",
            sample_id, batch_id
        )
        metadata = {
            "sample_id": sample_id,
            "batch_id": batch_id,
            "read_length": None,
            "primer_protocol": None,
            "location_code": None,
            "sampling_date": None,
            "location_name": None,
        }
    
    return metadata
