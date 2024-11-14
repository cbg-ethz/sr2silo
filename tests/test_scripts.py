"""Tests for the scripts in the scripts directory."""

from __future__ import annotations

import sys
from pathlib import Path

# Add the scripts directory to the Python path
scripts_dir = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(scripts_dir))

# Import the process_directory function from vp_transformer.py
from vp_transformer import get_metadata  # noqa: E402 # pyright: ignore
from vp_transformer import process_directory  # noqa: E402 # pyright: ignore


def test_get_metadata():
    """Test the get_metadata function."""
    metadata = get_metadata(
        Path("tests/data/samples/A1_05_2024_10_08/20241024_2411515907")
    )

    expected_metadata = {
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "20241024_2411515907",
        "sequencing_well_position": "A1",
        "location_code": "05",
        "sampling_date": "2024_10_08",
        "sequencing_date": "20241024",
        "flow_cell_serial_number": "2411515907",
    }

    assert metadata == expected_metadata


def test_process_directory():
    """Test the process_directory function."""
    process_directory(
        Path("tests/data/samples/A1_05_2024_10_08/20241024_2411515907"),
        Path("tests/output"),
        "nextstrain/sars-cov-2/wuhan-hu-1/orfs",
    )
    assert True
