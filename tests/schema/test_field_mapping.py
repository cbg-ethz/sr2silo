"""Tests for field mapping."""

from __future__ import annotations

import pytest

from sr2silo.schema import get_field_mapping


def test_get_field_mapping_covid():
    """Test field mapping for COVID organism."""
    mapping = get_field_mapping("covid")

    # Check expected mappings
    assert mapping["sample_id"] == "sampleId"
    assert mapping["batch_id"] == "batchId"
    assert mapping["location_code"] == "locationCode"
    assert mapping["sampling_date"] == "samplingDate"
    assert mapping["location_name"] == "locationName"
    assert mapping["sr2silo_version"] == "sr2siloVersion"

    # Ensure it's the complete expected mapping
    assert len(mapping) == 6


def test_get_field_mapping_unknown_organism():
    """Test that unknown organism raises ValueError."""
    with pytest.raises(ValueError) as exc_info:
        get_field_mapping("unknown_organism")

    assert "not supported" in str(exc_info.value)
    assert "unknown_organism" in str(exc_info.value)


def test_get_field_mapping_returns_copy():
    """Test that get_field_mapping returns a copy, not the original."""
    mapping1 = get_field_mapping("covid")
    mapping2 = get_field_mapping("covid")

    # Modify one mapping
    mapping1["new_field"] = "newField"

    # Ensure the other is not affected
    assert "new_field" not in mapping2
