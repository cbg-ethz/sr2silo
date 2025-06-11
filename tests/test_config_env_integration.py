"""Tests for the updated config module environment variable integration."""

from __future__ import annotations

import os
from unittest.mock import patch

import sr2silo.config as config


def test_timeline_file_env_var():
    """Test get_timeline_file with environment variable and default."""
    # Test with default when env var not set
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_timeline_file("default.tsv")
        assert result == "default.tsv"
    
    # Test with environment variable set
    with patch.dict(os.environ, {"TIMELINE_FILE": "/path/to/timeline.tsv"}):
        result = config.get_timeline_file("default.tsv")
        assert result == "/path/to/timeline.tsv"
    
    # Test with no default and no env var
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_timeline_file()
        assert result is None


def test_primer_file_env_var():
    """Test get_primer_file with environment variable and default."""
    # Test with default when env var not set
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_primer_file("default.yaml")
        assert result == "default.yaml"
    
    # Test with environment variable set
    with patch.dict(os.environ, {"PRIMER_FILE": "/path/to/primer.yaml"}):
        result = config.get_primer_file("default.yaml")
        assert result == "/path/to/primer.yaml"
    
    # Test with no default and no env var
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_primer_file()
        assert result is None


def test_reference_env_var():
    """Test get_reference with environment variable and default."""
    # Test with default when env var not set
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_reference("custom-ref")
        assert result == "custom-ref"
    
    # Test with built-in default
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_reference()
        assert result == "sars-cov-2"
    
    # Test with environment variable set
    with patch.dict(os.environ, {"NEXTCLADE_REFERENCE": "env-ref"}):
        result = config.get_reference("default-ref")
        assert result == "env-ref"


def test_keycloak_token_url_with_defaults():
    """Test get_keycloak_token_url with default parameter."""
    # Test with default when env var not set
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_keycloak_token_url("http://default.keycloak.url")
        assert result == "http://default.keycloak.url"
    
    # Test with environment variable set
    with patch.dict(os.environ, {"KEYCLOAK_TOKEN_URL": "http://env.keycloak.url"}):
        result = config.get_keycloak_token_url("http://default.keycloak.url")
        assert result == "http://env.keycloak.url"
    
    # Test error when neither env var nor default provided
    with patch.dict(os.environ, {}, clear=True):
        try:
            config.get_keycloak_token_url()
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "KEYCLOAK_TOKEN_URL environment variable is not set" in str(e)


def test_submission_url_with_defaults():
    """Test get_submission_url with default parameter."""
    # Test with default when env var not set
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_submission_url("http://default.submission.url")
        assert result == "http://default.submission.url"
    
    # Test with environment variable set
    with patch.dict(os.environ, {"SUBMISSION_URL": "http://env.submission.url"}):
        result = config.get_submission_url("http://default.submission.url")
        assert result == "http://env.submission.url"
    
    # Test error when neither env var nor default provided
    with patch.dict(os.environ, {}, clear=True):
        try:
            config.get_submission_url()
            assert False, "Should have raised ValueError"
        except ValueError as e:
            assert "SUBMISSION_URL environment variable is not set" in str(e)


def test_env_with_default_helper():
    """Test the helper function get_env_with_default."""
    # Test with environment variable set
    with patch.dict(os.environ, {"TEST_VAR": "env_value"}):
        result = config.get_env_with_default("TEST_VAR", "default_value")
        assert result == "env_value"
    
    # Test with default when env var not set
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_env_with_default("TEST_VAR", "default_value")
        assert result == "default_value"
    
    # Test with no default and no env var
    with patch.dict(os.environ, {}, clear=True):
        result = config.get_env_with_default("TEST_VAR")
        assert result is None