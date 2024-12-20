"""Implement tests for the metadata extraction functions."""


from __future__ import annotations

from sr2silo.vpipe.metadata import sample_id_decoder


def test_sample_id_decoder():
    """Test the sample_id_decoder function."""
    sample_id = "A1_10_2024_09_30"
    result = sample_id_decoder(sample_id)
    expected = {
        "sequencing_well_position": "A1",
        "location_code": "10",
        "sampling_date": "2024-09-30",
    }
    assert result == expected
