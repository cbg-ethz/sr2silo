"""Configuration utilities."""

from __future__ import annotations

import importlib.metadata
import logging
import os
import subprocess
import sys
from pathlib import Path

from dotenv import load_dotenv

# Load .env file from the project root
env_path = Path(__file__).parent.parent.parent / ".env"
load_dotenv(dotenv_path=env_path)
logging.info(f"Loaded environment variables from {env_path} using python-dotenv")


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
    url = os.getenv("KEYCLOAK_TOKEN_URL")
    if url is None:
        logging.error("KEYCLOAK_TOKEN_URL environment variable is not set")
        sys.exit(1)
    return url


def get_submission_url() -> str:
    """Get the LAPIS submission URL from environment, or return default if not set.

    Returns:
        str: The submission URL with group_id placeholder
    """
    url = os.getenv("SUBMISSION_URL")
    if url is None:
        logging.error("SUBMISSION_URL environment variable is not set")
        sys.exit(1)
    return url


def get_organism() -> str:
    """Get the organism identifier from environment, or return default if not set.

    Returns:
        str: The organism identifier (e.g., 'sc2', 'sars-cov-2')
    """
    return os.getenv("ORGANISM", "sc2")


def get_frontend_url() -> str:
    """Get the frontend URL by deriving it from the submission URL.

    Returns:
        str: The frontend URL (removes 'backend-' prefix from submission URL)

    Raises:
        ValueError: If the submission URL doesn't contain 'backend-' prefix
    """
    submission_url = get_submission_url()

    # Extract the base URL and remove 'backend-' prefix
    # Example: https://backend-wise-seqs.loculus.org/... -> https://wise-seqs.loculus.org/...
    if "backend-" not in submission_url:
        raise ValueError(
            """
            Cannot derive frontend URL:
            submission URL doesn't contain 'backend-' prefix
            """
        )

    frontend_url = submission_url.replace("backend-", "")
    # Remove the path part and just keep the domain
    if "//" not in frontend_url:
        raise ValueError("Cannot derive frontend URL: invalid URL format")

    protocol_and_domain = frontend_url.split("//")[1].split("/")[0]
    return f"https://{protocol_and_domain}"


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
