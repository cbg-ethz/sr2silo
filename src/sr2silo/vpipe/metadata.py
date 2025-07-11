"""Extract metadata from Timeline file only."""

from __future__ import annotations

import csv
import datetime
import logging
from pathlib import Path


def convert_to_iso_date(date: str) -> str:
    """Convert a date string to ISO 8601 format (date only).

    Args:
        date (str): Date string in YYYY-MM-DD format.

    Returns:
        str: ISO 8601 formatted date string.

    Raises:
        ValueError: If date string cannot be parsed or is empty.
    """
    if not date or not date.strip():
        raise ValueError("Date string cannot be empty")

    try:
        # Parse the date string
        date_obj = datetime.datetime.strptime(date.strip(), "%Y-%m-%d")
        # Format the date as ISO 8601 (date only)
        return date_obj.date().isoformat()
    except ValueError as e:
        raise ValueError(f"Invalid date format '{date}': {e}") from e


def get_metadata_from_timeline(sample_id: str, timeline: Path) -> dict[str, str] | None:
    """Get metadata from the timeline file.

    Args:
        sample_id (str): The sample ID to search for.
        timeline (Path): The timeline file to search in.

    Returns:
        dict[str, str] | None: The metadata if found, None otherwise.

    Raises:
        ValueError: If critical fields (sample_id, location_name, sampling_date) are
                   missing or invalid.
        FileNotFoundError: If timeline file is not found.
    """
    if not timeline.is_file():
        logging.error(f"Timeline file not found or is not a file: {timeline}")
        raise FileNotFoundError(f"Timeline file not found or is not a file: {timeline}")

    with timeline.open() as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            # Check if row has enough columns
            if len(row) < 7:
                logging.warning(
                    f"Skipping malformed row with {len(row)} columns: {row}"
                )
                continue

            if row[0] == sample_id:
                logging.info(
                    "Found metadata in timeline for sample_id %s",
                    sample_id,
                )

                # Validate critical fields before processing
                # Critical fields: sample_id, location_name (row[6]), sampling_date
                if not sample_id or not sample_id.strip():
                    raise ValueError("Sample ID cannot be empty")

                if not row[6] or not row[6].strip():
                    raise ValueError(f"Location name is missing for sample {sample_id}")

                if not row[5] or not row[5].strip():
                    raise ValueError(f"Sampling date is missing for sample {sample_id}")

                # Extract metadata from timeline row
                # Timeline format: sample	batch	reads	proto	location_code	date	location
                try:
                    sampling_date = convert_to_iso_date(row[5])
                except ValueError as e:
                    raise ValueError(
                        f"Invalid sampling date for sample {sample_id}: {e}"
                    ) from e

                metadata = {
                    "sample_id": sample_id,
                    "batch_id": row[1] if row[1] and row[1].strip() else None,
                    "read_length": row[2] if row[2] and row[2].strip() else None,
                    "primer_protocol": row[3] if row[3] and row[3].strip() else None,
                    "location_code": row[4] if row[4] and row[4].strip() else None,
                    "sampling_date": sampling_date,
                    "location_name": row[6].strip(),
                }

                logging.info("Extracted metadata from timeline: %s", metadata)
                return metadata

        # No matching entry found
        logging.warning(
            "No matching entry found in timeline for sample_id %s",
            sample_id,
        )
        return None


def get_metadata(sample_id: str, timeline: Path) -> dict[str, str]:
    """
    Get metadata for a given sample from timeline file only.

    Args:
        sample_id (str): The sample ID to use for metadata.
        timeline (Path): The timeline file to cross-reference the metadata.

    Returns:
        dict: A dictionary containing the metadata, or empty dict if not found.

    Raises:
        ValueError: If critical fields (sample_id, location_name, sampling_date) are
                   missing or invalid when found in timeline.
    """
    # Validate critical input
    if not sample_id or not sample_id.strip():
        raise ValueError("Sample ID cannot be empty")

    metadata = get_metadata_from_timeline(sample_id, timeline)

    if metadata is None:
        # Return basic metadata structure if not found in timeline
        logging.warning(
            "Timeline entry not found for sample_id %s. "
            "Returning basic metadata structure with None values.",
            sample_id,
        )
        metadata = {
            "sample_id": sample_id,
            "batch_id": None,
            "read_length": None,
            "primer_protocol": None,
            "location_code": None,
            "sampling_date": None,
            "location_name": None,
        }

    return metadata
