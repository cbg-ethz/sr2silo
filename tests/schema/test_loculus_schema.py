"""Tests for Loculus schema loading."""

from __future__ import annotations

import pytest

from sr2silo.schema.loculus_schema import LoculusSchema

# Sample Loculus configuration for testing (based on WisePulse config)
MOCK_LOCULUS_CONFIG = {
    "schema": {
        "organismName": "SARS-CoV-2",
        "metadata": [
            {
                "name": "date",
                "type": "date",
                "displayName": "Submission Date",
            },
            {
                "name": "sampleId",
                "type": "string",
                "displayName": "Sample ID",
                "required": True,
            },
            {
                "name": "batchId",
                "type": "string",
                "displayName": "Batch ID",
            },
            {
                "name": "locationCode",
                "type": "string",
                "displayName": "Location Code",
            },
            {
                "name": "locationName",
                "type": "string",
                "displayName": "Location",
            },
            {
                "name": "samplingDate",
                "type": "date",
                "displayName": "Sampling Date",
            },
            {
                "name": "sr2siloVersion",
                "type": "string",
                "displayName": "sr2silo Version",
            },
            {
                "name": "countSiloReads",
                "type": "string",
                "displayName": "srSILO Read Count",
            },
        ],
        "submissionDataTypes": {
            "consensusSequences": False,
            "files": {
                "enabled": True,
                "categories": [
                    {"name": "nucleotideAlignment"},
                    {"name": "siloReads"},
                ],
            },
        },
    }
}


@pytest.fixture
def loculus_schema():
    """Create a LoculusSchema instance."""
    return LoculusSchema.from_dict("covid", MOCK_LOCULUS_CONFIG)


def test_loculus_schema_from_dict():
    """Test creating LoculusSchema from dict."""
    schema = LoculusSchema.from_dict("covid", MOCK_LOCULUS_CONFIG)
    assert schema.organism == "covid"
    assert schema._config == MOCK_LOCULUS_CONFIG


def test_get_metadata_fields(loculus_schema):
    """Test extracting metadata field names."""
    fields = loculus_schema.get_metadata_fields()

    expected = [
        "date", "sampleId", "batchId", "locationCode",
        "locationName", "samplingDate", "sr2siloVersion", "countSiloReads"
    ]
    assert fields == expected


def test_get_metadata_field_info(loculus_schema):
    """Test getting detailed metadata field information."""
    field_info = loculus_schema.get_metadata_field_info()

    assert len(field_info) == 8
    # Check first field
    assert field_info[0]["name"] == "date"
    assert field_info[0]["type"] == "date"
    assert field_info[0]["displayName"] == "Submission Date"


def test_get_field_type(loculus_schema):
    """Test getting field type."""
    assert loculus_schema.get_field_type("sampleId") == "string"
    assert loculus_schema.get_field_type("date") == "date"
    assert loculus_schema.get_field_type("samplingDate") == "date"
    assert loculus_schema.get_field_type("nonexistent") is None


def test_get_required_fields(loculus_schema):
    """Test getting required fields."""
    required = loculus_schema.get_required_fields()

    # In our mock, only sampleId is explicitly marked as required
    # Others default to required=True
    assert "sampleId" in required
    assert len(required) > 0


def test_get_file_categories(loculus_schema):
    """Test getting file categories."""
    categories = loculus_schema.get_file_categories()

    assert categories == ["nucleotideAlignment", "siloReads"]


def test_get_file_categories_disabled():
    """Test getting file categories when files are disabled."""
    config = {
        "schema": {
            "submissionDataTypes": {
                "files": {
                    "enabled": False,
                }
            }
        }
    }
    schema = LoculusSchema.from_dict("covid", config)
    categories = schema.get_file_categories()

    assert categories == []


def test_get_metadata_fields_no_config():
    """Test error when config not loaded."""
    schema = LoculusSchema("covid", config=None)

    with pytest.raises(ValueError) as exc_info:
        schema.get_metadata_fields()
    assert "Configuration not loaded" in str(exc_info.value)


def test_from_file_missing_organism(tmp_path):
    """Test error when organism not in config file."""
    # Create a YAML file with different organism
    config_file = tmp_path / "config.yml"
    config_file.write_text("organisms:\n  mpox:\n    schema: {}")

    with pytest.raises(ValueError) as exc_info:
        LoculusSchema.from_file("covid", config_file)
    assert "not found in config" in str(exc_info.value)


def test_from_file_not_exists():
    """Test error when config file doesn't exist."""
    from pathlib import Path
    
    config_file = Path("/nonexistent/config.yml")

    with pytest.raises(FileNotFoundError):
        LoculusSchema.from_file("covid", config_file)


def test_from_file_no_yaml():
    """Test error when PyYAML not available."""
    import sys
    from unittest.mock import patch
    
    # Mock yaml module as None (not installed)
    with patch.dict(sys.modules, {'yaml': None}):
        # Need to reimport to trigger the None check
        import importlib
        from sr2silo.schema import loculus_schema as ls_module
        importlib.reload(ls_module)
        
        # Now test should raise ImportError
        from pathlib import Path
        config_file = Path("/tmp/test.yml")
        
        with pytest.raises(ImportError) as exc_info:
            ls_module.LoculusSchema.from_file("covid", config_file)
        assert "PyYAML is required" in str(exc_info.value)
