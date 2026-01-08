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
        error_msg = "Location name is missing for sample A1_05_2024_10_08"
        with pytest.raises(ValueError, match=error_msg):
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
            "A1_05_2024_10_08\t20241024_2411515907\t250\tv532\t5\t"
            "invalid-date\tLugano (TI)\n"
        )
        temp_timeline = Path(f.name)

    try:
        error_msg = "Invalid sampling date for sample A1_05_2024_10_08"
        with pytest.raises(ValueError, match=error_msg):
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


def test_get_metadata_rsva_organism(rsva_timeline):
    """Test the get_metadata function with RSV-A organism-specific columns.

    RSV uses different column names: submissionId instead of sample,
    primerProtocol instead of proto.
    """
    metadata = get_metadata(
        sample_id="A1_05_2025_11_05",
        timeline=rsva_timeline,
        organism="rsva",
    )

    # Check that metadata was extracted correctly using RSV column names
    assert metadata["sample_id"] == "A1_05_2025_11_05"
    assert metadata["batch_id"] == "20251128_2511665243"
    assert metadata["read_length"] == "250"
    assert metadata["primer_protocol"] == "Eawag-2024-v532_pooled"
    assert metadata["location_code"] == "05"
    assert metadata["sampling_date"] == "2025-11-05"
    assert metadata["location_name"] == "Lugano"


def test_get_metadata_rsva_vs_covid_columns():
    """Test that RSV and COVID parse correctly despite different column names.

    This ensures the column mapping system works correctly for both organisms.
    """
    # Create a COVID-style timeline (sample, proto columns)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation\n")
        f.write(
            "COVID_SAMPLE\t20241024_2411515907\t250\tCOVID_PROTO\t5\t"
            "2024-10-08\tLugano (TI)\n"
        )
        covid_timeline = Path(f.name)

    # Create an RSV-style timeline (submissionId, primerProtocol columns)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write(
            "submissionId\tbatch\treads\tprimerProtocol\tlocation_code\tdate\tlocation\n"
        )
        f.write(
            "RSV_SAMPLE\t20251128_2511665243\t250\tRSV_PROTO\t05\t2025-11-05\tLugano\n"
        )
        rsv_timeline = Path(f.name)

    try:
        # Parse COVID timeline with COVID organism (default)
        covid_metadata = get_metadata(
            sample_id="COVID_SAMPLE",
            timeline=covid_timeline,
            organism="covid",
        )

        # Parse RSV timeline with RSV organism
        rsv_metadata = get_metadata(
            sample_id="RSV_SAMPLE",
            timeline=rsv_timeline,
            organism="rsva",
        )

        # Both should parse successfully with correct values
        assert covid_metadata["sample_id"] == "COVID_SAMPLE"
        assert covid_metadata["primer_protocol"] == "COVID_PROTO"
        assert covid_metadata["batch_id"] == "20241024_2411515907"

        assert rsv_metadata["sample_id"] == "RSV_SAMPLE"
        assert rsv_metadata["primer_protocol"] == "RSV_PROTO"
        assert rsv_metadata["batch_id"] == "20251128_2511665243"

    finally:
        covid_timeline.unlink()
        rsv_timeline.unlink()


def test_get_metadata_malformed_row_skipped():
    """Test that malformed rows (too few columns) are skipped."""

    # Create a temporary timeline file with malformed row
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
        f.write("sample\tbatch\treads\tproto\tlocation_code\tdate\tlocation\n")
        f.write(
            "A1_05_2024_10_08\t20241024_2411515907\t250\n"
        )  # only 3 columns - malformed
        # correct row with all required fields
        f.write(
            "A1_05_2024_10_08\t20241024_2411515907\t250\tv532\t5\t"
            "2024-10-08\tLugano (TI)\n"
        )
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
