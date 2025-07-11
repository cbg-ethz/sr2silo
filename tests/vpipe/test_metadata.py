"""Implement tests for the metadata extraction functions."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest

from sr2silo.vpipe.metadata import get_metadata


def test_get_metadata(timeline: Path):
    """Test the get_metadata function."""

    metadata = get_metadata(
        sample_id="A1_05_2024_10_08",
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
        timeline=timeline,
    )

    expected_metadata = {
        "sample_id": "NONEXISTENT",
        "batch_id": None,
        "read_length": None,
        "primer_protocol": None,
        "location_code": None,
        "sampling_date": None,
        "location_name": None,
    }

    assert metadata == expected_metadata


def test_get_metadata_finds_sample(timeline: Path):
    """Test the get_metadata function finds sample correctly."""

    metadata = get_metadata(
        sample_id="A1_05_2024_10_08",
        timeline=timeline,
    )

    expected_metadata = {
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "20241024_2411515907",  # Should use the one from timeline
        "read_length": "250",
        "primer_protocol": "v532",
        "location_code": "5",
        "sampling_date": "2024-10-08",
        "location_name": "Lugano (TI)",
    }

    assert metadata == expected_metadata


def test_get_metadata_single_lookup(timeline: Path):
    """Test the get_metadata function works with sample-only lookup."""

    metadata = get_metadata(
        sample_id="A1_05_2024_10_08",
        timeline=timeline,
    )

    expected_metadata = {
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "20241024_2411515907",  # Should use the one from timeline
        "read_length": "250",
        "primer_protocol": "v532",
        "location_code": "5",
        "sampling_date": "2024-10-08",
        "location_name": "Lugano (TI)",
    }

    assert metadata == expected_metadata


def test_get_metadata_empty_sample_id_error(timeline: Path):
    """Test that empty sample_id raises ValueError."""

    with pytest.raises(ValueError, match="Sample ID cannot be empty"):
        get_metadata(
            sample_id="",
            timeline=timeline,
        )

    with pytest.raises(ValueError, match="Sample ID cannot be empty"):
        get_metadata(
            sample_id="  ",  # whitespace only
            timeline=timeline,
        )


def test_get_metadata_missing_location_error():
    """Test that missing location_name raises ValueError."""

    # Create a temporary timeline file with missing location
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation\n")
        f.write(
            "A1_05_2024_10_08\t20241024_2411515907\t250\tv532\t5\t2024-10-08\t\n"
        )  # empty location
        temp_timeline = Path(f.name)

    try:
        with pytest.raises(
            ValueError, match="Location name is missing for sample A1_05_2024_10_08"
        ):
            get_metadata(
                sample_id="A1_05_2024_10_08",
                timeline=temp_timeline,
            )
    finally:
        temp_timeline.unlink()


def test_get_metadata_missing_date_error():
    """Test that missing sampling_date raises ValueError."""

    # Create a temporary timeline file with missing date
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation\n")
        f.write(
            "A1_05_2024_10_08\t20241024_2411515907\t250\tv532\t5\t\tLugano (TI)\n"
        )  # empty date
        temp_timeline = Path(f.name)

    try:
        with pytest.raises(
            ValueError, match="Sampling date is missing for sample A1_05_2024_10_08"
        ):
            get_metadata(
                sample_id="A1_05_2024_10_08",
                timeline=temp_timeline,
            )
    finally:
        temp_timeline.unlink()


def test_get_metadata_invalid_date_error():
    """Test that invalid sampling_date raises ValueError."""

    # Create a temporary timeline file with invalid date
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation\n")
        f.write(
            "A1_05_2024_10_08\t20241024_2411515907\t250\tv532\t5\tinvalid-date\tLugano (TI)\n"
        )
        temp_timeline = Path(f.name)

    try:
        with pytest.raises(
            ValueError, match="Invalid sampling date for sample A1_05_2024_10_08"
        ):
            get_metadata(
                sample_id="A1_05_2024_10_08",
                timeline=temp_timeline,
            )
    finally:
        temp_timeline.unlink()


def test_get_metadata_graceful_empty_fields():
    """Test that non-critical empty fields are handled gracefully."""

    # Create a temporary timeline file with some empty non-critical fields
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation\n")
        f.write(
            "A1_05_2024_10_08\t\t\t\t\t2024-10-08\tLugano (TI)\n"
        )  # empty batch, reads, proto, location_code
        temp_timeline = Path(f.name)

    try:
        metadata = get_metadata(
            sample_id="A1_05_2024_10_08",
            timeline=temp_timeline,
        )

        expected_metadata = {
            "sample_id": "A1_05_2024_10_08",
            "batch_id": None,  # empty in timeline
            "read_length": None,  # empty in timeline
            "primer_protocol": None,  # empty in timeline
            "location_code": None,  # empty in timeline
            "sampling_date": "2024-10-08",
            "location_name": "Lugano (TI)",
        }

        assert metadata == expected_metadata
    finally:
        temp_timeline.unlink()


def test_get_metadata_malformed_row_skipped():
    """Test that malformed rows (too few columns) are skipped."""

    # Create a temporary timeline file with malformed row
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation\n")
        f.write(
            "A1_05_2024_10_08\t20241024_2411515907\t250\n"
        )  # only 3 columns - malformed
        f.write(
            "A1_05_2024_10_08\t20241024_2411515907\t250\tv532\t5\t2024-10-08\tLugano (TI)\n"
        )  # correct row
        temp_timeline = Path(f.name)

    try:
        metadata = get_metadata(
            sample_id="A1_05_2024_10_08",
            timeline=temp_timeline,
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
    finally:
        temp_timeline.unlink()
