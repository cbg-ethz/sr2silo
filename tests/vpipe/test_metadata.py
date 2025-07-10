"""Implement tests for the metadata extraction functions."""

from __future__ import annotations

from pathlib import Path

from sr2silo.vpipe.metadata import get_metadata


def test_get_metadata(timeline: Path):
    """Test the get_metadata function."""

    metadata = get_metadata(
        sample_id="A1_05_2024_10_08",
        batch_id="20241024_2411515907",
        timeline=timeline,
    )

    expected_metadata = {
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "20241024_2411515907",
        "read_length": "250",
        "primer_protocol": "v532",
        "location_code": "5",
        "sampling_date": "2024-10-08",
        "location_name": "Lugano (TI)",
    }

    assert metadata == expected_metadata


def test_get_metadata_not_found(timeline: Path):
    """Test the get_metadata function when no matching entry is found."""

    metadata = get_metadata(
        sample_id="NONEXISTENT",
        batch_id="NONEXISTENT",
        timeline=timeline,
    )

    expected_metadata = {
        "sample_id": "NONEXISTENT",
        "batch_id": "NONEXISTENT",
        "read_length": None,
        "primer_protocol": None,
        "location_code": None,
        "sampling_date": None,
        "location_name": None,
    }

    assert metadata == expected_metadata
