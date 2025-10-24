"""Organism-specific field mapping configuration.

This module provides simple field mappings from SILO ReadMetadata fields
(snake_case) to Loculus submission fields (camelCase) for different organisms.

The mappings are based on the ReadMetadata schema and can be easily extended
for new organisms.
"""

from __future__ import annotations

from typing import Dict

# Mapping from SILO ReadMetadata fields to Loculus submission fields
# Based on src/sr2silo/silo_read_schema.ReadMetadata
ORGANISM_FIELD_MAPPINGS: Dict[str, Dict[str, str]] = {
    "covid": {
        "sample_id": "sampleId",
        "batch_id": "batchId",
        "location_code": "locationCode",
        "sampling_date": "samplingDate",
        "location_name": "locationName",
        "sr2silo_version": "sr2siloVersion",
    },
    # Add more organisms here as needed
    # "mpox": {
    #     "sample_id": "sampleId",
    #     ...
    # },
}


def get_field_mapping(organism: str) -> Dict[str, str]:
    """Get field mapping for a specific organism.

    Args:
        organism: Organism identifier (e.g., 'covid', 'mpox')

    Returns:
        Dictionary mapping SILO field names (snake_case) to Loculus field names (camelCase)

    Raises:
        ValueError: If organism is not supported
    """
    if organism not in ORGANISM_FIELD_MAPPINGS:
        raise ValueError(
            f"Organism '{organism}' not supported. "
            f"Available organisms: {list(ORGANISM_FIELD_MAPPINGS.keys())}"
        )
    return ORGANISM_FIELD_MAPPINGS[organism].copy()
