"""Tests for loculus module functions."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path

import pytest

from sr2silo.loculus.loculus import Submission


@pytest.fixture
def test_silo_input_uncompressed():
    """Fixture providing path to uncompressed test SILO input file."""
    return Path(__file__).parent / "test_silo_input.ndjson"


@pytest.fixture
def test_silo_input_compressed():
    """Fixture providing path to compressed test SILO input file."""
    return Path(__file__).parent / "test_silo_input.ndjson.zst"


def test_parse_metadata_uncompressed(test_silo_input_uncompressed):
    """Test parsing metadata from uncompressed SILO input file."""
    metadata = Submission.parse_metadata(test_silo_input_uncompressed)

    # Check that metadata is returned as dictionary
    assert isinstance(metadata, dict)

    # Check expected metadata fields in snake_case (excluding read_id)
    expected_fields = {
        "sample_id": "A1_05_2025_06_18",
        "batch_id": "20250711_2443602573",
        "location_code": "5",
        "location_name": "Lugano (TI)",
        "sampling_date": "2025-06-18",
        "sr2silo_version": "1.3.0 (v1.1.0-5-g31b0623)",
    }

    for key, expected_value in expected_fields.items():
        assert key in metadata
        assert metadata[key] == expected_value

    # Ensure read_id is excluded (function should exclude it)
    assert "read_id" not in metadata


def test_parse_metadata_compressed(test_silo_input_compressed):
    """Test parsing metadata from compressed SILO input file."""
    metadata = Submission.parse_metadata(test_silo_input_compressed)

    # Check that metadata is returned as dictionary
    assert isinstance(metadata, dict)

    # Check expected metadata fields in snake_case (excluding read_id)
    expected_fields = {
        "sample_id": "A1_05_2025_06_18",
        "batch_id": "20250711_2443602573",
        "location_code": "5",
        "location_name": "Lugano (TI)",
        "sampling_date": "2025-06-18",
        "sr2silo_version": "1.3.0 (v1.1.0-5-g31b0623)",
    }

    for key, expected_value in expected_fields.items():
        assert key in metadata
        assert metadata[key] == expected_value

    # Ensure read_id is excluded (function should exclude it)
    assert "read_id" not in metadata


def test_count_reads_uncompressed(test_silo_input_uncompressed):
    """Test counting reads in uncompressed SILO input file."""
    count = Submission.count_reads(test_silo_input_uncompressed)

    # The test file has 2 lines (no final newline), so expect 2 reads
    assert count == 3


def test_count_reads_compressed(test_silo_input_compressed):
    """Test counting reads in compressed SILO input file."""
    count = Submission.count_reads(test_silo_input_compressed)

    # The test file has 3 lines, so expect 3 reads
    assert count == 3


def test_create_metadata_file(test_silo_input_uncompressed):
    """Test creating metadata TSV file."""
    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy test file to temp directory to simulate processed file
        temp_processed_file = Path(temp_dir) / "test_processed.ndjson"
        temp_processed_file.write_text(test_silo_input_uncompressed.read_text())

        # Test without count_reads
        metadata_file, submission_id = Submission.create_metadata_file(
            temp_processed_file
        )

        # Check that file was created
        assert metadata_file.exists()
        assert metadata_file.name.startswith("metadata_")
        assert metadata_file.suffix == ".tsv"

        # Check that submission directory was created
        assert metadata_file.parent.name == "submission"

        # Check file contents
        with open(metadata_file, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            fieldnames = reader.fieldnames or []

            # Should have exactly 1 row of data (wide format)
            assert len(rows) == 1
            assert "submissionId" in fieldnames
            assert "date" in fieldnames

            # Check submission ID matches
            first_row = rows[0]
            assert first_row["submissionId"] == submission_id
            assert first_row["date"]  # Should have a date

            # Check that metadata fields are present as columns
            expected_metadata_fields = {
                "sampleId",
                "batchId",
                "locationCode",
                "locationName",
                "samplingDate",
                "sr2siloVersion",
            }
            for field in expected_metadata_fields:
                assert field in fieldnames, (
                    f"Expected metadata field '{field}' not found in columns"
                )
                assert first_row[field]  # Should have a value

        # Test with count_reads=True
        metadata_file2, submission_id2 = Submission.create_metadata_file(
            temp_processed_file, count_reads=True
        )

        # Should create a different file with different submission ID
        assert metadata_file2 != metadata_file
        assert submission_id2 != submission_id
        assert metadata_file2.exists()

        # Check that countReads column is included when count_reads=True
        with open(metadata_file2, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = list(reader)
            fieldnames = reader.fieldnames or []

            assert "countSiloReads" in fieldnames
            first_row = rows[0]
            assert first_row["countSiloReads"]  # Should have a count value
            assert int(first_row["countSiloReads"]) > 0  # Should be a positive number
