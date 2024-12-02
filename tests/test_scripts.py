"""Tests for the scripts in the scripts directory."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# Add the scripts directory to the Python path
scripts_dir = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(scripts_dir))


# TODO: refactor the metadata formatter out of the scripts
@pytest.mark.skip(
    reason="Won't run on GitHub Actions, as it requires the rust extension"
)
def test_get_metadata():
    """Test the get_metadata function."""

    # Import the process_directory function from vp_transformer.py
    from vp_transformer import get_metadata  # noqa: E402 # pyright: ignore

    metadata = get_metadata(
        sample_id="A1_05_2024_10_08",
        batch_id="20241024_2411515907",
        timeline=Path("tests/data/samples/timeline_A1_05_2024_10_08.tsv"),
        primers=Path("tests/data/samples/primers.yaml"),
    )

    print(metadata)

    expected_metadata = {
        "sample_id": "A1_05_2024_10_08",
        "batch_id": "20241024_2411515907",
        "sequencing_well_position": "A1",
        "location_code": "05",
        "sampling_date": "2024-10-08",
        "sequencing_date": "2024-10-24",
        "flow_cell_serial_number": "2411515907",
        "read_length": "250",
        "primer_protocol": "v532",
        "location_name": "Lugano (TI)",
        "primer_protocol_name": "SARS-CoV-2 ARTIC V5.3.2",
    }

    assert metadata == expected_metadata


# TODO: refactor the function outside of the scripts
@pytest.mark.skip(
    reason="Won't run on GitHub Actions, as it requires the rust extension"
)
def test_process_directory():
    """Test the process_directory function."""

    # Import the process_directory function from vp_transformer.py
    from vp_transformer import (
        process_directory,
    )  # noqa: E402 # pyright: ignore # pylint: disable=all, isort:skip, black:skip

    process_directory(
        Path("tests/data/samples/A1_05_2024_10_08/20241024_2411515907/alignments"),
        "A1_05_2024_10_08",
        "20241024_2411515907",
        Path("tests/output"),
        "nextstrain/sars-cov-2/wuhan-hu-1/orfs",
        Path("tests/data/samples/timeline_A1_05_2024_10_08.tsv"),
        Path("tests/data/samples/primers.yaml"),
    )
    assert True
