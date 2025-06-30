"""Tests for the config module."""

from __future__ import annotations

import importlib
import importlib.metadata
import os
import subprocess
from pathlib import Path
from unittest.mock import patch

from sr2silo.config import (
    get_keycloak_token_url,
    get_primer_file,
    get_submission_url,
    get_timeline_file,
    get_version,
    is_ci_environment,
)


def test_is_ci_environment():
    """Test the is_ci_environment function."""
    # Test with CI=true
    with patch.dict(os.environ, {"CI": "true"}):
        assert is_ci_environment() is True

    # Test with CI=false
    with patch.dict(os.environ, {"CI": "false"}):
        assert is_ci_environment() is False

    # Test with CI not set
    with patch.dict(os.environ, {}, clear=True):
        assert is_ci_environment() is False


def test_get_version_only_package():
    """Test get_version without git info."""
    with patch("importlib.metadata.version", return_value="1.2.3"):
        # Without git info
        version = get_version(add_git_info=False)
        assert version == "1.2.3"


def test_get_version_in_ci():
    """Test get_version in CI environment."""
    with patch("importlib.metadata.version", return_value="1.2.3"):
        with patch("sr2silo.config.is_ci_environment", return_value=True):
            # In CI environment
            version = get_version()
            assert version == "1.2.3"


def test_get_version_with_git_success():
    """Test get_version with successful git command."""
    # Standard format from poetry-dynamic-versioning could be:
    # - Tagged release: 1.2.3
    # - Development build: 1.2.3.dev4+a1b2c3d
    with patch("importlib.metadata.version", return_value="1.2.3.dev4+a1b2c3d"):
        with patch("sr2silo.config.is_ci_environment", return_value=False):
            with patch("subprocess.check_output", return_value="v1.2.3-4-g12345678"):
                version = get_version()
                # Even with dynamic versioning, the git info is still appended
                assert version == "1.2.3.dev4+a1b2c3d (v1.2.3-4-g12345678)"


def test_get_version_with_git_failure():
    """Test get_version when git command fails."""
    with patch("importlib.metadata.version", return_value="1.2.3"):
        with patch("sr2silo.config.is_ci_environment", return_value=False):
            with patch(
                "subprocess.check_output",
                side_effect=subprocess.CalledProcessError(1, "git"),
            ):
                version = get_version()
                assert version == "1.2.3"


def test_get_version_package_not_found():
    """Test get_version when package is not found."""
    with patch(
        "importlib.metadata.version",
        side_effect=importlib.metadata.PackageNotFoundError(),
    ):
        # Without git info
        version = get_version(add_git_info=False)
        assert version == "unknown"

        # With git info but in CI environment
        with patch("sr2silo.config.is_ci_environment", return_value=True):
            version = get_version()
            assert version == "unknown"

        # With git info but git command fails
        with patch("sr2silo.config.is_ci_environment", return_value=False):
            with patch(
                "subprocess.check_output",
                side_effect=subprocess.CalledProcessError(1, "git"),
            ):
                version = get_version()
                assert version == "unknown"


def test_get_timeline_file():
    """Test get_timeline_file function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"TIMELINE_FILE": "/path/to/timeline.tsv"}):
        result = get_timeline_file()
        assert result == Path("/path/to/timeline.tsv")

    # Test with environment variable not set, with default
    with patch.dict(os.environ, {}, clear=True):
        result = get_timeline_file("/default/timeline.tsv")
        assert result == Path("/default/timeline.tsv")

    # Test with environment variable not set, no default
    with patch.dict(os.environ, {}, clear=True):
        result = get_timeline_file()
        assert result is None

    # Test with Path default
    with patch.dict(os.environ, {}, clear=True):
        result = get_timeline_file(Path("/path/default.tsv"))
        assert result == Path("/path/default.tsv")


def test_get_primer_file():
    """Test get_primer_file function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"PRIMER_FILE": "/path/to/primers.yaml"}):
        result = get_primer_file()
        assert result == Path("/path/to/primers.yaml")

    # Test with environment variable not set, with default
    with patch.dict(os.environ, {}, clear=True):
        result = get_primer_file("/default/primers.yaml")
        assert result == Path("/default/primers.yaml")

    # Test with environment variable not set, no default
    with patch.dict(os.environ, {}, clear=True):
        result = get_primer_file()
        assert result is None


def test_get_keycloak_token_url():
    """Test get_keycloak_token_url function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"KEYCLOAK_TOKEN_URL": "https://auth.example.com"}):
        result = get_keycloak_token_url()
        assert result == "https://auth.example.com"

    # Test with environment variable not set, with default
    with patch.dict(os.environ, {}, clear=True):
        result = get_keycloak_token_url("https://default.auth.com")
        assert result == "https://default.auth.com"

    # Test with environment variable not set, no default - should exit
    with patch.dict(os.environ, {}, clear=True):
        with patch("sys.exit") as mock_exit:
            get_keycloak_token_url()
            mock_exit.assert_called_once_with(1)


def test_get_submission_url():
    """Test get_submission_url function."""
    # Test with environment variable set
    with patch.dict(os.environ, {"SUBMISSION_URL": "https://submit.example.com"}):
        result = get_submission_url()
        assert result == "https://submit.example.com"

    # Test with environment variable not set, with default
    with patch.dict(os.environ, {}, clear=True):
        result = get_submission_url("https://default.submit.com")
        assert result == "https://default.submit.com"

    # Test with environment variable not set, no default - should exit
    with patch.dict(os.environ, {}, clear=True):
        with patch("sys.exit") as mock_exit:
            get_submission_url()
            mock_exit.assert_called_once_with(1)
