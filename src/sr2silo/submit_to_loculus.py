"""Submit processd data to Loculus Silo and upload to S3."""

from __future__ import annotations

import json
import logging
from pathlib import Path

from sr2silo.config import (
    get_group_id,
    get_keycloak_token_url,
    get_mock_urls,
    get_organism,
    get_password,
    get_submission_url,
    get_username,
    is_ci_environment,
)
from sr2silo.loculus import LoculusClient, Submission


def load_config(config_file: Path) -> dict:
    """Load a JSON configuration file."""
    try:
        with config_file.open() as f:
            return json.load(f)
    except FileNotFoundError:
        logging.error(f"Config file not found: {config_file}")
        raise
    except json.JSONDecodeError as e:
        logging.error(f"Error decoding JSON from config file: {config_file} - {e}")
        raise


def submit_to_silo(
    processed_file: Path,
    keycloak_token_url: str | None = None,
    submission_url: str | None = None,
    group_id: int | None = None,
    username: str | None = None,
    password: str | None = None,
) -> bool:
    """Submit data to SILO using the new pre-signed upload approach.

    Args:
        result_dir (Path): Directory where to save submission files.
        processed_file (Path): Path to the processed .ndjson.zst file to upload.
        keycloak_token_url (str | None): Keycloak token URL. If None, uses environment.
        submission_url (str | None): Submission URL. If None, uses environment.
        group_id (int | None): Group ID for submission. If None, uses environment.
        username (str | None): Username for authentication. If None, uses environment.
        password (str | None): Password for authentication. If None, uses environment.

    Returns:
        bool: True if submission was successful, False otherwise.
    """
    logging.info("Submitting to Loculus...")

    # Create the new metadata file format and get the submission ID
    metadata_fp, submission_id = Submission.create_metadata_file(
        processed_file, count_reads=True
    )

    if is_ci_environment():
        logging.info(
            "CI environment active; using mock URLs but executing submission mechanics."
        )
        # Get mock URLs for CI environment, then continue with LAPIS submission
        KEYCLOAK_TOKEN_URL, SUBMISSION_URL = get_mock_urls()
    else:
        # Get URLs from parameters or environment
        KEYCLOAK_TOKEN_URL = keycloak_token_url or get_keycloak_token_url()
        SUBMISSION_URL = submission_url or get_submission_url()
        logging.info(f"Using Keycloak URL: {KEYCLOAK_TOKEN_URL}")
        logging.info(f"Using submission URL: {SUBMISSION_URL}")

    # Resolve authentication parameters with environment fallback
    resolved_group_id = group_id if group_id is not None else get_group_id()
    resolved_username = username or get_username()
    resolved_password = password or get_password()

    # Get organism configuration
    organism = get_organism()
    logging.info(f"Using organism: {organism}")
    logging.info(f"Using group ID: {resolved_group_id}")
    logging.info(f"Using username: {resolved_username}")

    try:
        # Create client with organism parameter
        client = LoculusClient(KEYCLOAK_TOKEN_URL, SUBMISSION_URL, organism)
        client.authenticate(username=resolved_username, password=resolved_password)

        # Submit using new API with both metadata and processed file
        response = client.submit(
            group_id=resolved_group_id,
            metadata_file_path=metadata_fp,
            processed_file_path=processed_file,
            submission_id=submission_id,
        )

        if response["status"] == "success":
            logging.info(response["message"])
            return True
        else:
            logging.error(f"Submission failed: {response}")
            return False

    except Exception as e:
        logging.error(f"Error during submission: {e}")
        return False
