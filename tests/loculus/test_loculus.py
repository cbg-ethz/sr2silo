"""Tests for loculus module functions."""

from __future__ import annotations

from pathlib import Path

import pytest

from sr2silo.loculus.loculus import Submission


@pytest.fixture
def test_silo_input_uncompressed():
    """Fixture providing path to uncompressed test SILO input file."""
    return Path(__file__).parent / "test_silo_input.ndjson"


def test_parse_metadata_uncompressed(test_silo_input_uncompressed):
    """Test parsing metadata from uncompressed SILO input file."""
    metadata = Submission.parse_metadata(test_silo_input_uncompressed)

    # Check that metadata is returned as dictionary
    assert isinstance(metadata, dict)

    # Check expected metadata fields in camelCase (excluding readId)
    expected_fields = {
        "sampleId": "A1_05_2025_06_18",
        "batchId": "20250711_2443602573",
        "locationCode": "5",
        "locationName": "Lugano (TI)",
        "samplingDate": "2025-06-18",
        "sr2siloVersion": "1.3.0 (v1.1.0-5-g31b0623)",
    }

    for key, expected_value in expected_fields.items():
        assert key in metadata
        assert metadata[key] == expected_value

    # Ensure readId is excluded (function should exclude it)
    assert "readId" not in metadata


def test_count_reads_uncompressed(test_silo_input_uncompressed):
    """Test counting reads in uncompressed SILO input file."""
    count = Submission.count_reads(test_silo_input_uncompressed)

    # The test file has 2 lines (no final newline), so expect 2 reads
    assert count == 3
