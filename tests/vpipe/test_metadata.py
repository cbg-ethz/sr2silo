"""Implement tests for the metadata extraction functions."""


from __future__ import annotations

from sr2silo.vpipe.metadata import batch_id_decoder, sample_id_decoder


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


def test_batch_id_decoder():
    """Test the batch_id_decoder function."""
    batch_id = "20241018_AAG55WNM5"
    result = batch_id_decoder(batch_id)
    expected = {
        "sequencing_date": "2024-10-18",
        "flow_cell_serial_number": "AAG55WNM5",
    }
    assert result == expected
