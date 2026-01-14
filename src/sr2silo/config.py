"""Configuration utilities."""

from __future__ import annotations

import importlib.metadata
import logging
import os
import sys
from pathlib import Path

import yaml


def is_ci_environment() -> bool:
    """Check if running in a CI environment."""
    return bool(os.getenv("CI", "false").lower() in ("yes", "true", "t", "1"))


def get_version() -> str:
    """Get the version information of the sr2silo package.

    Returns:
        str: The version information (package version only for reliability)
    """
    try:
        # Get package version from metadata
        package_version = importlib.metadata.version("sr2silo")
    except importlib.metadata.PackageNotFoundError:
        package_version = "unknown"

    # For reliable versioning, we use only the package version
    # This ensures consistency between development and production
    # The package version is manually updated at release time
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


def get_backend_url(default: str | None = None) -> str:
    """Get the backend URL from environment, or return default if not set.

    Args:
        default: Optional default value to use if environment variable is not set

    Returns:
        str: The backend URL with group_id placeholder

    Raises:
        SystemExit: If no URL is available from environment or default
    """
    url = os.getenv("BACKEND_URL", default)
    if url is None:
        logging.error(
            "BACKEND_URL environment variable is not set and no default provided"
        )
        sys.exit(1)
    return url


def get_organism() -> str:
    """Get the organism identifier from environment, or return default if not set.

    Returns:
        str: The organism identifier (e.g., 'covid')
    """
    try:
        organism = os.getenv("ORGANISM")
        if organism:
            return organism
        else:
            logging.error("ORGANISM environment variable is not set.")
            sys.exit(1)
    except Exception as e:
        logging.error(f"Error retrieving ORGANISM from environment: {e}")
        sys.exit(1)


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


def get_auto_release(default: bool = False) -> bool:
    """Get the auto-release setting from environment.

    Args:
        default: Default value if environment variable is not set (default: False)

    Returns:
        bool: Whether to automatically release sequences after submission
    """
    auto_release = os.getenv("AUTO_RELEASE", "").lower()
    if auto_release in ("yes", "true", "t", "1"):
        return True
    elif auto_release in ("no", "false", "f", "0", ""):
        return default
    else:
        logging.warning(
            f"Invalid AUTO_RELEASE value '{auto_release}', using default: {default}"
        )
        return default


def get_mock_urls() -> tuple[str, str]:
    """Get mock URLs for CI environment.

    Returns:
        tuple[str, str]: The Keycloak token URL and backend URL
    """
    mock_keycloak_url = (
        "https://authentication-wise-seqs.loculus.org/realms/loculus/"
        "protocol/openid-connect/token"
    )
    mock_backend_url = (
        "https://backend-wise-seqs.loculus.org/test/submit?"
        "groupId={group_id}&dataUseTermsType=OPEN"
    )
    return mock_keycloak_url, mock_backend_url


def get_timeline_column_mappings(organism: str) -> dict[str, str]:
    """Get timeline column name mappings for a specific organism.

    Returns a dictionary mapping internal field names to timeline TSV column names.
    Uses organism-specific configuration from resources/vpipe/timeline_columns.yml.

    Args:
        organism: The organism identifier (e.g., 'covid', 'rsva')

    Returns:
        dict[str, str]: Mapping from internal names to timeline column names.
                       Default mappings are used if organism config not found.

    Examples:
        >>> mappings = get_timeline_column_mappings('covid')
        >>> mappings['sample_id']  # Returns 'sample'
        >>> mappings = get_timeline_column_mappings('rsva')
        >>> mappings['sample_id']  # Returns 'submissionId'
    """
    # Default mappings (backward compatible with COVID timeline format)
    default_mappings = {
        "sample_id": "sample",
        "batch_id": "batch",
        "read_length": "reads",
        "primer_protocol": "proto",
        "location_code": "location_code",
        "sampling_date": "date",
        "location_name": "location",
    }

    # Try to load organism-specific mappings from config file
    try:
        # Find the timeline columns config file
        config_path = (
            Path(__file__).parent.parent.parent
            / "resources"
            / "vpipe"
            / "timeline_columns.yml"
        )

        if not config_path.exists():
            logging.warning(
                f"Timeline columns config file not found at {config_path}. "
                "Using default timeline column mappings."
            )
            return default_mappings

        with open(config_path, "r") as f:
            config = yaml.safe_load(f)

        # Get organism-specific timeline column mappings
        if config and "organisms" in config and organism in config["organisms"]:
            mappings = config["organisms"][organism]
            logging.debug(f"Loaded timeline column mappings for {organism}: {mappings}")
            return mappings
        else:
            logging.warning(
                f"No timeline column mappings found for organism '{organism}'. "
                "Using default mappings."
            )
            return default_mappings

    except Exception as e:
        logging.warning(
            f"Error loading timeline columns config: {e}. "
            "Using default timeline column mappings."
        )
        return default_mappings
