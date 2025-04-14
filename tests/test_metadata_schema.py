"""Tests for metadata schema validation."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from sr2silo.silo_read_schema import ReadMetadata

# Sample valid metadata dictionary simulating a database record
VALID_METADATA = {
    "read_id": "read123",
    "sample_id": "A1_05_2024_10_08",
    "batch_id": "batch001",
    "sampling_date": "2024-10-08",
    "sequencing_date": "2024-10-24",
    "location_name": "Lugano (TI)",
    "read_length": "250",
    "primer_protocol": "v532",
    "location_code": "05",
    "flow_cell_serial_number": "2411515907",
    "sequencing_well_position": "A1",
    "primer_protocol_name": "SARS-CoV-2 ARTIC V5.3.2",
    "nextclade_reference": "sars-cov-2",
    "sr2silo_version": "v1.0.0",
}


def test_valid_metadata():
    # Validate that a valid metadata dictionary creates a ReadMetadata instance
    metadata = ReadMetadata(**VALID_METADATA)
    for key, value in VALID_METADATA.items():
        assert getattr(metadata, key) == value


def test_invalid_metadata_missing_field():
    # Remove a required field to simulate an invalid record
    invalid_metadata = dict(VALID_METADATA)
    invalid_metadata.pop("read_id")
    with pytest.raises(ValidationError):
        ReadMetadata(**invalid_metadata)
