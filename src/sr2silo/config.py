"""Configuration utilities."""

from __future__ import annotations

import os


def is_ci_environment() -> bool:
    """Check if running in a CI environment."""
    return bool(os.getenv("CI", "false").lower() in ("yes", "true", "t", "1"))


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
