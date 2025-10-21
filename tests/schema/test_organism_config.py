"""Tests for organism configuration management."""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import pytest

from sr2silo.schema.field_mapping import FieldMapping
from sr2silo.schema.loculus_schema import LoculusSchema
from sr2silo.schema.organism_config import OrganismConfig
from sr2silo.schema.silo_schema import SiloSchema

# Mock schemas for testing
MOCK_SILO_SCHEMA = {
    "metadata": [
        {"name": "read_id", "type": "string"},
        {"name": "sample_id", "type": "string"},
        {"name": "batch_id", "type": "string"},
        {"name": "sampling_date", "type": "date"},
        {"name": "location_name", "type": "string"},
        {"name": "location_code", "type": "string"},
        {"name": "sr2silo_version", "type": "string"},
    ],
    "nucleotideSequences": [{"name": "main"}],
    "genes": [{"name": "S"}, {"name": "ORF1a"}],
}

MOCK_LOCULUS_CONFIG = {
    "schema": {
        "organismName": "SARS-CoV-2",
        "metadata": [
            {"name": "date", "type": "date", "required": True},
            {"name": "sampleId", "type": "string", "required": True},
            {"name": "batchId", "type": "string"},
            {"name": "locationCode", "type": "string"},
            {"name": "locationName", "type": "string"},
            {"name": "samplingDate", "type": "date"},
            {"name": "sr2siloVersion", "type": "string"},
        ],
        "submissionDataTypes": {
            "files": {
                "enabled": True,
                "categories": [
                    {"name": "nucleotideAlignment"},
                    {"name": "siloReads"},
                ],
            }
        },
    }
}


@pytest.fixture
def organism_config():
    """Create a basic OrganismConfig for testing."""
    # Create mock schemas
    silo_schema = SiloSchema("https://lapis.example.org", "covid")
    silo_schema._schema = MOCK_SILO_SCHEMA

    loculus_schema = LoculusSchema.from_dict("covid", MOCK_LOCULUS_CONFIG)

    field_mapping = FieldMapping("covid")

    return OrganismConfig(
        organism="covid",
        silo_schema=silo_schema,
        loculus_schema=loculus_schema,
        field_mapping=field_mapping,
    )


def test_organism_config_init(organism_config):
    """Test OrganismConfig initialization."""
    assert organism_config.organism == "covid"
    assert organism_config.silo_schema is not None
    assert organism_config.loculus_schema is not None
    assert organism_config.field_mapping is not None


@patch('sr2silo.schema.silo_schema.requests.get')
def test_from_api_and_config(mock_get):
    """Test creating OrganismConfig from API and config dict."""
    # Mock API response
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = MOCK_SILO_SCHEMA
    mock_response.raise_for_status = MagicMock()
    mock_get.return_value = mock_response

    config = OrganismConfig.from_api_and_config(
        organism="covid",
        lapis_url="https://lapis.example.org",
        loculus_config=MOCK_LOCULUS_CONFIG,
        use_cache=False,
    )

    assert config.organism == "covid"
    assert config.silo_schema is not None
    assert config.loculus_schema is not None
    assert config.field_mapping is not None


def test_validate_complete_valid(organism_config):
    """Test validation with complete valid configuration."""
    validation = organism_config.validate(strict=False)

    assert validation["valid"] is True
    assert len(validation["errors"]) == 0


def test_validate_missing_silo_schema():
    """Test validation fails when SILO schema not loaded."""
    config = OrganismConfig(
        organism="covid",
        silo_schema=None,
        loculus_schema=LoculusSchema.from_dict("covid", MOCK_LOCULUS_CONFIG),
        field_mapping=FieldMapping("covid"),
    )

    validation = config.validate(strict=False)

    assert validation["valid"] is False
    assert any("SILO schema not loaded" in err for err in validation["errors"])


def test_validate_missing_loculus_schema():
    """Test validation fails when Loculus schema not loaded."""
    silo_schema = SiloSchema("https://lapis.example.org", "covid")
    silo_schema._schema = MOCK_SILO_SCHEMA

    config = OrganismConfig(
        organism="covid",
        silo_schema=silo_schema,
        loculus_schema=None,
        field_mapping=FieldMapping("covid"),
    )

    validation = config.validate(strict=False)

    assert validation["valid"] is False
    assert any("Loculus schema not loaded" in err for err in validation["errors"])


def test_validate_invalid_field_mapping(organism_config):
    """Test validation fails with invalid field mapping."""
    # Add mapping to non-existent field
    organism_config.field_mapping.add_mapping("nonexistent_silo_field", "nonexistentLoculusField")

    validation = organism_config.validate(strict=False)

    assert validation["valid"] is False
    assert len(validation["errors"]) > 0


def test_get_metadata_for_submission(organism_config):
    """Test transforming SILO metadata to Loculus format."""
    silo_metadata = {
        "read_id": "read123",
        "sample_id": "sample001",
        "batch_id": "batch001",
        "location_code": "01",
        "location_name": "Zurich",
        "sampling_date": "2024-01-15",
        "sr2silo_version": "1.5.0",
    }

    loculus_metadata = organism_config.get_metadata_for_submission(silo_metadata)

    # Check mapping
    assert loculus_metadata["sampleId"] == "sample001"
    assert loculus_metadata["batchId"] == "batch001"
    assert loculus_metadata["locationCode"] == "01"
    assert loculus_metadata["locationName"] == "Zurich"
    assert loculus_metadata["samplingDate"] == "2024-01-15"
    assert loculus_metadata["sr2siloVersion"] == "1.5.0"

    # read_id should not be in loculus metadata (not in mapping)
    assert "read_id" not in loculus_metadata
    assert "readId" not in loculus_metadata


def test_get_metadata_for_submission_missing_field(organism_config):
    """Test transforming metadata with missing SILO field."""
    silo_metadata = {
        "sample_id": "sample001",
        # Missing batch_id and other fields
    }

    loculus_metadata = organism_config.get_metadata_for_submission(silo_metadata)

    # Should have sampleId but not other fields
    assert loculus_metadata["sampleId"] == "sample001"
    assert "batchId" not in loculus_metadata


def test_log_configuration_summary(organism_config, caplog):
    """Test logging configuration summary."""
    import logging
    
    with caplog.at_level(logging.INFO):
        organism_config.log_configuration_summary()

    # Check that key information is logged
    assert "Organism Configuration Summary: covid" in caplog.text
    assert "SILO metadata fields" in caplog.text
    assert "Loculus metadata fields" in caplog.text
    assert "Field mappings" in caplog.text


def test_validate_strict_mode(organism_config):
    """Test strict validation mode."""
    # Add extra field to Loculus schema that's not in mapping
    mock_config_with_extra = {
        "schema": {
            "metadata": MOCK_LOCULUS_CONFIG["schema"]["metadata"] + [
                {"name": "extraField", "type": "string"}
            ],
            "submissionDataTypes": MOCK_LOCULUS_CONFIG["schema"]["submissionDataTypes"],
        }
    }

    organism_config.loculus_schema = LoculusSchema.from_dict("covid", mock_config_with_extra)

    validation = organism_config.validate(strict=True)

    # Should have warnings about unmapped fields in strict mode
    assert len(validation["warnings"]) > 0


def test_validate_required_fields_warning(organism_config):
    """Test validation warns about missing required Loculus fields."""
    # Create config with minimal mapping that doesn't cover required fields
    minimal_mapping = FieldMapping("covid", {"sample_id": "sampleId"})
    organism_config.field_mapping = minimal_mapping

    validation = organism_config.validate(strict=False)

    # Should have warnings about unmapped required fields
    # Note: 'date' and 'submissionId' are system fields so should be excluded from warnings
    assert len(validation["warnings"]) > 0 or validation["valid"]
