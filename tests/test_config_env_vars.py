"""Tests for environment variable configuration functionality."""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import patch

import pytest

from sr2silo.config import (
    get_batch_id,
    get_default_input_file,
    get_nextclade_reference,
    get_primer_file,
    get_results_dir,
    get_sample_dir,
    get_sample_id,
    get_timeline_file,
)


def test_get_sample_dir():
    """Test get_sample_dir function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"SAMPLE_DIR": "/path/to/sample"}):
        assert get_sample_dir() == "/path/to/sample"

    # Test with environment variable not set
    with patch.dict(os.environ, {}, clear=True):
        assert get_sample_dir() is None


def test_get_sample_id():
    """Test get_sample_id function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"SAMPLE_ID": "test_sample"}):
        assert get_sample_id() == "test_sample"

    # Test with environment variable not set
    with patch.dict(os.environ, {}, clear=True):
        assert get_sample_id() is None


def test_get_batch_id():
    """Test get_batch_id function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"BATCH_ID": "test_batch"}):
        assert get_batch_id() == "test_batch"

    # Test with environment variable not set
    with patch.dict(os.environ, {}, clear=True):
        assert get_batch_id() is None


def test_get_timeline_file():
    """Test get_timeline_file function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"TIMELINE_FILE": "/path/to/timeline.tsv"}):
        assert get_timeline_file() == "/path/to/timeline.tsv"

    # Test with environment variable not set
    with patch.dict(os.environ, {}, clear=True):
        assert get_timeline_file() is None


def test_get_primer_file():
    """Test get_primer_file function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"PRIMER_FILE": "/path/to/primers.yaml"}):
        assert get_primer_file() == "/path/to/primers.yaml"

    # Test with environment variable not set
    with patch.dict(os.environ, {}, clear=True):
        assert get_primer_file() is None


def test_get_nextclade_reference():
    """Test get_nextclade_reference function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"NEXTCLADE_REFERENCE": "custom-ref"}):
        assert get_nextclade_reference() == "custom-ref"

    # Test with environment variable not set (should return default)
    with patch.dict(os.environ, {}, clear=True):
        assert get_nextclade_reference() == "sars-cov-2"


def test_get_results_dir():
    """Test get_results_dir function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"RESULTS_DIR": "/path/to/results"}):
        assert get_results_dir() == "/path/to/results"

    # Test with environment variable not set
    with patch.dict(os.environ, {}, clear=True):
        assert get_results_dir() is None


def test_get_default_input_file():
    """Test get_default_input_file function."""
    # Test with SAMPLE_DIR set
    with patch.dict(os.environ, {"SAMPLE_DIR": "/path/to/sample"}):
        result = get_default_input_file()
        assert result == Path("/path/to/sample/REF_aln_trim.bam")

    # Test with SAMPLE_DIR not set
    with patch.dict(os.environ, {}, clear=True):
        assert get_default_input_file() is None