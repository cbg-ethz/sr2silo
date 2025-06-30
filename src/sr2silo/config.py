"""Configuration utilities."""

from __future__ import annotations

import importlib.metadata
import logging
import os
import subprocess
import sys
from pathlib import Path

from dotenv import load_dotenv

# Load .env file from the project root if it exists
env_path = Path(__file__).parent.parent.parent / ".env"
if env_path.exists():
    load_dotenv(dotenv_path=env_path)
    logging.info(f"Loaded environment variables from {env_path} using python-dotenv")
else:
    logging.debug(
        f"No .env file found at {env_path}, using only system environment variables"
    )


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


def get_keycloak_token_url(default: str | None = None) -> str:
    """Get the Keycloak token URL from environment, or return default if not set.

    Args:
        default: Optional default value to use if environment variable is not set

    Returns:
        str: The Keycloak token URL

    Raises:
        SystemExit: If no URL is available from environment or default
    """
    url = os.getenv("KEYCLOAK_TOKEN_URL", default)
    if url is None:
        logging.error(
            "KEYCLOAK_TOKEN_URL environment variable is not set and no default provided"
        )
        sys.exit(1)
    return url


def get_submission_url(default: str | None = None) -> str:
    """Get the LAPIS submission URL from environment, or return default if not set.

    Args:
        default: Optional default value to use if environment variable is not set

    Returns:
        str: The submission URL with group_id placeholder

    Raises:
        SystemExit: If no URL is available from environment or default
    """
    url = os.getenv("SUBMISSION_URL", default)
    if url is None:
        logging.error(
            "SUBMISSION_URL environment variable is not set and no default provided"
        )
        sys.exit(1)
    return url


def get_organism() -> str:
    """Get the organism identifier from environment, or return default if not set.

    Returns:
        str: The organism identifier (e.g., 'sc2', 'sars-cov-2')
    """
    return os.getenv("ORGANISM", "sc2")


def get_timeline_file(default: Path | str | None = None) -> Path | None:
    """Get the timeline file path from environment, or return default if not set.

    Args:
        default: Optional default path to use if environment variable is not set

    Returns:
        Path | None: The timeline file path, or None if not available
    """
    timeline_file = os.getenv("TIMELINE_FILE")
    if timeline_file:
        return Path(timeline_file)
    elif default:
        return Path(default) if isinstance(default, str) else default
    return None


def get_primer_file(default: Path | str | None = None) -> Path | None:
    """Get the primer file path from environment, or return default if not set.

    Args:
        default: Optional default path to use if environment variable is not set

    Returns:
        Path | None: The primer file path, or None if not available
    """
    primer_file = os.getenv("PRIMER_FILE")
    if primer_file:
        return Path(primer_file)
    elif default:
        return Path(default) if isinstance(default, str) else default
    return None


def get_nextclade_reference(default: str = "sars-cov-2") -> str:
    """Get the Nextclade reference from environment, or return default if not set.

    Args:
        default: Default reference to use if environment variable is not set

    Returns:
        str: The Nextclade reference identifier
    """
    return os.getenv("NEXTCLADE_REFERENCE", default)


def get_group_id(default: int | None = None) -> int:
    """Get the group ID from environment, or return default if not set.

    Args:
        default: Optional default group ID to use if environment variable is not set

    Returns:
        int: The group ID for submissions

    Raises:
        SystemExit: If no group ID is available from environment or default
    """
    group_id_str = os.getenv("GROUP_ID")
    if group_id_str is None:
        if default is None:
            logging.error(
                "GROUP_ID environment variable is not set and no default provided"
            )
            sys.exit(1)
        return default

    try:
        return int(group_id_str)
    except ValueError:
        logging.error(f"Invalid GROUP_ID value '{group_id_str}', must be an integer")
        sys.exit(1)


def get_username(default: str | None = None) -> str:
    """Get the username from environment, or return default if not set.

    Args:
        default: Optional default username to use if environment variable is not set

    Returns:
        str: The username for authentication

    Raises:
        SystemExit: If no username is available from environment or default
    """
    username = os.getenv("USERNAME", default)
    if username is None:
        logging.error(
            "USERNAME environment variable is not set and no default provided"
        )
        sys.exit(1)
    return username


def get_password(default: str | None = None) -> str:
    """Get the password from environment, or return default if not set.

    Args:
        default: Optional default password to use if environment variable is not set

    Returns:
        str: The password for authentication

    Raises:
        SystemExit: If no password is available from environment or default
    """
    password = os.getenv("PASSWORD", default)
    if password is None:
        logging.error(
            "PASSWORD environment variable is not set and no default provided"
        )
        sys.exit(1)
    return password


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
