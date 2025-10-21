"""Tests for SILO schema fetching and caching."""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

from sr2silo.schema.silo_schema import SiloSchema

# Sample SILO schema for testing
MOCK_SILO_SCHEMA = {
    "metadata": [
        {"name": "read_id", "type": "string"},
        {"name": "sample_id", "type": "string"},
        {"name": "batch_id", "type": "string"},
        {"name": "sampling_date", "type": "date"},
        {"name": "location_name", "type": "string"},
        {"name": "location_code", "type": "string"},
    ],
    "nucleotideSequences": [
        {"name": "main", "sequence": "ACGT..."}
    ],
    "genes": [
        {"name": "S", "sequence": "MFV..."},
        {"name": "ORF1a", "sequence": "MES..."},
        {"name": "N", "sequence": "MPR..."},
    ]
}


@pytest.fixture
def temp_cache_dir(tmp_path):
    """Create a temporary cache directory."""
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    return cache_dir


@pytest.fixture
def silo_schema(temp_cache_dir):
    """Create a SiloSchema instance with temp cache."""
    return SiloSchema(
        lapis_url="https://lapis.example.org",
        organism="covid",
        cache_dir=temp_cache_dir
    )


def test_silo_schema_init(silo_schema, temp_cache_dir):
    """Test SiloSchema initialization."""
    assert silo_schema.lapis_url == "https://lapis.example.org"
    assert silo_schema.organism == "covid"
    assert silo_schema.cache_dir == temp_cache_dir
    assert silo_schema._schema is None


def test_cache_file_path(silo_schema, temp_cache_dir):
    """Test cache file path generation."""
    expected = temp_cache_dir / "silo_schema_covid.json"
    assert silo_schema.cache_file == expected


@patch('requests.get')
def test_fetch_schema_from_api(mock_get, silo_schema):
    """Test fetching schema from API."""
    # Mock successful API response
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = MOCK_SILO_SCHEMA
    mock_response.raise_for_status = MagicMock()
    mock_get.return_value = mock_response

    schema = silo_schema.fetch_schema(use_cache=False)

    assert schema == MOCK_SILO_SCHEMA
    assert silo_schema._schema == MOCK_SILO_SCHEMA
    mock_get.assert_called_once()


@patch('requests.get')
def test_fetch_schema_with_caching(mock_get, silo_schema):
    """Test fetching schema and caching it."""
    # Mock successful API response
    mock_response = MagicMock()
    mock_response.status_code = 200
    mock_response.json.return_value = MOCK_SILO_SCHEMA
    mock_response.raise_for_status = MagicMock()
    mock_get.return_value = mock_response

    # First fetch - should call API and cache
    schema = silo_schema.fetch_schema(use_cache=True)
    assert schema == MOCK_SILO_SCHEMA
    assert silo_schema.cache_file.exists()

    # Second fetch - should load from cache without API call
    silo_schema2 = SiloSchema(
        lapis_url="https://lapis.example.org",
        organism="covid",
        cache_dir=silo_schema.cache_dir
    )
    schema2 = silo_schema2.fetch_schema(use_cache=True)
    assert schema2 == MOCK_SILO_SCHEMA
    # API should only be called once (from first fetch)
    assert mock_get.call_count == 1


@patch('requests.get')
def test_fetch_schema_api_failure_with_cache_fallback(mock_get, silo_schema):
    """Test fallback to cache when API fails."""
    # First, populate cache
    silo_schema.cache_file.write_text(json.dumps(MOCK_SILO_SCHEMA))

    # Mock API failure
    mock_get.side_effect = Exception("Network error")

    # Should fall back to cache
    schema = silo_schema.fetch_schema(use_cache=True)
    assert schema == MOCK_SILO_SCHEMA


@patch('requests.get')
def test_fetch_schema_api_failure_no_cache(mock_get, silo_schema):
    """Test exception when API fails and no cache available."""
    # Mock API failure
    mock_get.side_effect = Exception("Network error")

    # Should raise exception
    with pytest.raises(Exception) as exc_info:
        silo_schema.fetch_schema(use_cache=True)
    assert "Failed to fetch SILO schema" in str(exc_info.value)


def test_get_metadata_fields(silo_schema):
    """Test extracting metadata fields from schema."""
    silo_schema._schema = MOCK_SILO_SCHEMA
    fields = silo_schema.get_metadata_fields()

    expected = ["read_id", "sample_id", "batch_id", "sampling_date", "location_name", "location_code"]
    assert fields == expected


def test_get_metadata_fields_not_loaded(silo_schema):
    """Test error when getting metadata fields before loading schema."""
    with pytest.raises(ValueError) as exc_info:
        silo_schema.get_metadata_fields()
    assert "Schema not loaded" in str(exc_info.value)


def test_get_nucleotide_sequences(silo_schema):
    """Test extracting nucleotide sequences from schema."""
    silo_schema._schema = MOCK_SILO_SCHEMA
    sequences = silo_schema.get_nucleotide_sequences()

    assert sequences == ["main"]


def test_get_genes(silo_schema):
    """Test extracting genes from schema."""
    silo_schema._schema = MOCK_SILO_SCHEMA
    genes = silo_schema.get_genes()

    assert genes == ["S", "ORF1a", "N"]


def test_clear_cache(silo_schema):
    """Test clearing the cache."""
    # Create cache file
    silo_schema.cache_file.write_text(json.dumps(MOCK_SILO_SCHEMA))
    assert silo_schema.cache_file.exists()

    # Clear cache
    silo_schema.clear_cache()
    assert not silo_schema.cache_file.exists()


def test_clear_cache_no_file(silo_schema):
    """Test clearing cache when file doesn't exist."""
    assert not silo_schema.cache_file.exists()
    # Should not raise error
    silo_schema.clear_cache()
    assert not silo_schema.cache_file.exists()
