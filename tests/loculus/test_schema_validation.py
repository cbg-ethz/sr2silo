"""Tests for schema validation in submission."""

from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
import yaml

from sr2silo.loculus.loculus import Submission


@pytest.fixture
def test_schema():
    """Fixture providing a test organism schema."""
    schema = {
        "schema": {
            "metadata": [
                {"name": "sample_id", "type": "string"},
                {"name": "sampleId", "type": "string"},
                {"name": "batch_id", "type": "string"},
                {"name": "batchId", "type": "string"},
                {"name": "location_code", "type": "string"},
                {"name": "locationCode", "type": "string"},
                {"name": "sampling_date", "type": "date"},
                {"name": "samplingDate", "type": "date"},
                {"name": "location_name", "type": "string"},
                {"name": "locationName", "type": "string"},
                {"name": "sr2silo_version", "type": "string"},
                {"name": "sr2siloVersion", "type": "string"},
                {"name": "read_id", "type": "string"},
                {"name": "read_length", "type": "string"},
                {"name": "primer_protocol", "type": "string"},
            ],
            "primaryKey": "read_id",
        }
    }
    return schema


@pytest.fixture
def test_schema_file(test_schema):
    """Fixture providing a test schema file."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False
    ) as schema_file:
        yaml.dump(test_schema, schema_file)
        schema_path = Path(schema_file.name)

    yield schema_path

    # Cleanup
    schema_path.unlink()


@pytest.fixture
def test_metadata():
    """Fixture providing test metadata."""
    return {
        "sample_id": "A1_05_2025_06_18",
        "batch_id": "20250711_2443602573",
        "location_code": "5",
        "location_name": "Lugano (TI)",
        "sampling_date": "2025-06-18",
        "sr2silo_version": "1.3.0",
    }


def test_validate_metadata_with_valid_schema(test_metadata, test_schema_file):
    """Test metadata validation with a valid schema."""
    is_valid, error_msg = Submission.validate_metadata_against_schema(
        test_metadata, test_schema_file
    )

    assert is_valid is True
    assert error_msg == ""


def test_validate_metadata_with_missing_schema_file(test_metadata):
    """Test metadata validation with missing schema file."""
    nonexistent_schema = Path("/tmp/nonexistent_schema.yaml")
    is_valid, error_msg = Submission.validate_metadata_against_schema(
        test_metadata, nonexistent_schema
    )

    assert is_valid is False
    assert "not found" in error_msg


def test_validate_metadata_with_invalid_yaml(test_metadata):
    """Test metadata validation with invalid YAML file."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False
    ) as schema_file:
        schema_file.write("invalid: yaml: content: [")
        schema_path = Path(schema_file.name)

    try:
        is_valid, error_msg = Submission.validate_metadata_against_schema(
            test_metadata, schema_path
        )

        assert is_valid is False
        assert "parsing" in error_msg.lower() or "yaml" in error_msg.lower()
    finally:
        schema_path.unlink()


def test_validate_metadata_with_invalid_schema_format(test_metadata):
    """Test metadata validation with invalid schema format (missing 'schema' key)."""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".yaml", delete=False
    ) as schema_file:
        yaml.dump({"invalid_key": "value"}, schema_file)
        schema_path = Path(schema_file.name)

    try:
        is_valid, error_msg = Submission.validate_metadata_against_schema(
            test_metadata, schema_path
        )

        assert is_valid is False
        assert "Invalid schema" in error_msg or "schema" in error_msg.lower()
    finally:
        schema_path.unlink()


def test_create_metadata_file_with_schema(test_schema_file):
    """Test creating metadata file with schema validation."""
    test_silo_input = Path(__file__).parent / "test_silo_input.ndjson"

    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy test file to temp directory
        temp_processed_file = Path(temp_dir) / "test_processed.ndjson"
        temp_processed_file.write_text(test_silo_input.read_text())

        # Create metadata file with schema validation
        metadata_file, submission_id = Submission.create_metadata_file(
            temp_processed_file, count_reads=False, schema_path=test_schema_file
        )

        # Check that file was created successfully
        assert metadata_file.exists()
        assert metadata_file.name.startswith("metadata_")
        assert submission_id is not None


def test_create_metadata_file_without_schema():
    """Test creating metadata file without schema validation."""
    test_silo_input = Path(__file__).parent / "test_silo_input.ndjson"

    with tempfile.TemporaryDirectory() as temp_dir:
        # Copy test file to temp directory
        temp_processed_file = Path(temp_dir) / "test_processed.ndjson"
        temp_processed_file.write_text(test_silo_input.read_text())

        # Create metadata file without schema validation
        metadata_file, submission_id = Submission.create_metadata_file(
            temp_processed_file, count_reads=False, schema_path=None
        )

        # Check that file was created successfully
        assert metadata_file.exists()
        assert metadata_file.name.startswith("metadata_")
        assert submission_id is not None
