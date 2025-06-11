"""Configuration utilities."""

from __future__ import annotations

import importlib.metadata
import os
import subprocess
from pathlib import Path


def is_ci_environment() -> bool:
    """Check if running in a CI environment."""
    return bool(os.getenv("CI", "false").lower() in ("yes", "true", "t", "1"))


def get_version(add_git_info: bool = True) -> str:
    """Get the version information of the sr2silo package.

    Args:
        add_git_info (bool): Whether to include Git version information.
            Will be ignored in CI environments.

    Returns:
        str: The version information (package version and Git info if available)
    """
    try:
        # Get package version from metadata
        package_version = importlib.metadata.version("sr2silo")
    except importlib.metadata.PackageNotFoundError:
        package_version = "unknown"

    # If in CI environment or explicitly disabled, return only package version
    if is_ci_environment() or not add_git_info:
        return package_version

    # Try to get Git version info
    try:
        git_version = subprocess.check_output(
            ["git", "describe", "--tags", "--match", "v*"],
            text=True,
            stderr=subprocess.STDOUT,
            timeout=2,  # Timeout after 2 seconds
        ).strip()
        return f"{package_version} ({git_version})"
    except (
        subprocess.CalledProcessError,
        subprocess.TimeoutExpired,
        FileNotFoundError,
    ):
        # Fall back to package version if Git command fails
        return package_version


def get_keycloak_token_url() -> str:
    """Get the Keycloak token URL from environment, or return default if not set.

    Returns:
        str: The Keycloak token URL
    """
    default_url = (
        "https://authentication-wise-seqs.loculus.org/realms/loculus/"
        "protocol/openid-connect/token"
    )
    return os.getenv("KEYCLOAK_TOKEN_URL", default_url)


def get_submission_url() -> str:
    """Get the LAPIS submission URL from environment, or return default if not set.

    Returns:
        str: The submission URL with group_id placeholder
    """
    default_url = (
        "https://backend-wise-seqs.loculus.org/test/submit?"
        "groupId={group_id}&dataUseTermsType=OPEN"
    )
    return os.getenv("SUBMISSION_URL", default_url)


def get_mock_urls() -> tuple[str, str]:
    """Get mock URLs for CI environment.

    Returns:
        tuple[str, str]: The Keycloak token URL and submission URL
    """
    mock_keycloak_url = (
        "https://authentication-wise-seqs.loculus.org/realms/loculus/"
        "protocol/openid-connect/token"
    )
    mock_submission_url = (
        "https://backend-wise-seqs.loculus.org/test/submit?"
        "groupId={group_id}&dataUseTermsType=OPEN"
    )
    return mock_keycloak_url, mock_submission_url


def get_sample_dir() -> str | None:
    """Get the sample directory from environment.

    Returns:
        str | None: The sample directory path, or None if not set
    """
    return os.getenv("SAMPLE_DIR")


def get_sample_id() -> str | None:
    """Get the sample ID from environment.

    Returns:
        str | None: The sample ID, or None if not set
    """
    return os.getenv("SAMPLE_ID")


def get_batch_id() -> str | None:
    """Get the batch ID from environment.

    Returns:
        str | None: The batch ID, or None if not set
    """
    return os.getenv("BATCH_ID")


def get_timeline_file() -> str | None:
    """Get the timeline file path from environment.

    Returns:
        str | None: The timeline file path, or None if not set
    """
    return os.getenv("TIMELINE_FILE")


def get_primer_file() -> str | None:
    """Get the primer file path from environment.

    Returns:
        str | None: The primer file path, or None if not set
    """
    return os.getenv("PRIMER_FILE")


def get_nextclade_reference() -> str:
    """Get the Nextclade reference from environment, or return default if not set.

    Returns:
        str: The Nextclade reference
    """
    return os.getenv("NEXTCLADE_REFERENCE", "sars-cov-2")


def get_results_dir() -> str | None:
    """Get the results directory from environment.

    Returns:
        str | None: The results directory path, or None if not set
    """
    return os.getenv("RESULTS_DIR")


def get_default_input_file() -> Path | None:
    """Get the default input file path based on SAMPLE_DIR environment variable.
    
    Looks for the standard V-PIPE alignment file 'REF_aln_trim.bam' in the sample directory.

    Returns:
        Path | None: The default input file path, or None if SAMPLE_DIR is not set
    """
    sample_dir = get_sample_dir()
    if sample_dir:
        return Path(sample_dir) / "REF_aln_trim.bam"
    return None
