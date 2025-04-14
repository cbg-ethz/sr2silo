"""Tests for the config module."""

from __future__ import annotations

import importlib
import importlib.metadata
import os
import subprocess
from unittest.mock import patch

from sr2silo.config import get_version, is_ci_environment


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
