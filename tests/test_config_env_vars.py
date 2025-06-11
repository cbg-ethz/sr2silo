"""Tests for environment variable configuration support."""

from __future__ import annotations

import os
from unittest.mock import patch

from sr2silo.config import get_primer_file, get_reference, get_timeline_file


def test_get_timeline_file_with_env_var():
    """Test get_timeline_file with environment variable set."""
    with patch.dict(os.environ, {"TIMELINE_FILE": "/path/to/timeline.tsv"}):
        assert get_timeline_file() == "/path/to/timeline.tsv"


def test_get_timeline_file_without_env_var():
    """Test get_timeline_file without environment variable set."""
    with patch.dict(os.environ, {}, clear=True):
        assert get_timeline_file() is None


def test_get_primer_file_with_env_var():
    """Test get_primer_file with environment variable set."""
    with patch.dict(os.environ, {"PRIMER_FILE": "/path/to/primers.yaml"}):
        assert get_primer_file() == "/path/to/primers.yaml"


def test_get_primer_file_without_env_var():
    """Test get_primer_file without environment variable set."""
    with patch.dict(os.environ, {}, clear=True):
        assert get_primer_file() is None


def test_get_reference_with_env_var():
    """Test get_reference with environment variable set."""
    with patch.dict(os.environ, {"NEXTCLADE_REFERENCE": "custom-reference"}):
        assert get_reference() == "custom-reference"


def test_get_reference_without_env_var():
    """Test get_reference without environment variable set (should return default)."""
    with patch.dict(os.environ, {}, clear=True):
        assert get_reference() == "sars-cov-2"


def test_get_reference_with_empty_env_var():
    """Test get_reference with empty environment variable (should return default)."""
    with patch.dict(os.environ, {"NEXTCLADE_REFERENCE": ""}):
        assert get_reference() == "sars-cov-2"